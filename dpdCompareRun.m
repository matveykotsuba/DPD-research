function results = dpdCompareRun(cfg, trainBeforeCompare)

    if nargin < 1 || isempty(cfg)
        cfg = config();
    end

    if nargin < 2 || isempty(trainBeforeCompare)
        trainBeforeCompare = false;
    end

    if ~isfield(cfg, 'pa')
        error('dpdCompareRun:MissingPaConfig', 'cfg.pa is required.');
    end

    cfrCfg = systemCfrConfiguration(cfg);

    if ~isfield(cfg, 'dpd') || isempty(cfg.dpd)
        cfg.dpd = dpdConfig(cfg.pa, true, 'indirect', cfrCfg);
    end

    cfg.pa.enabled = true;
    dpdCfg = cfg.dpd;
    [dpdCfg, ~] = dpdRefreshCfrCompatibility(dpdCfg, cfrCfg);
    cfg.dpd = dpdCfg;
    trainingInfo = struct();

    architecture = learningArchitecture(dpdCfg);
    defaultDpdCfg = dpdConfig(cfg.pa, false, architecture, cfrCfg);
    if isfield(dpdCfg, 'modelFile') && ~isDefaultArchitectureModelFile(dpdCfg.modelFile)
        modelFile = dpdCfg.modelFile;
    else
        modelFile = defaultDpdCfg.modelFile;
    end
    modelFileExists = exist(modelFile, 'file') == 2;
    trainedCoefficientsAvailable = hasTrainedCoefficients(dpdCfg);

    shouldTrain = logical(trainBeforeCompare);

    if ~shouldTrain && ~trainedCoefficientsAvailable && ~modelFileExists
        if exist('dpdTrainRun', 'file') == 2
            shouldTrain = true;
            fprintf(['DPD comparison: no trained %s model was found; ' 'running dpdTrainRun first.\n'], learningArchitecture(dpdCfg));
        else
            error('dpdCompareRun:NoTrainedDpd', ['No trained architecture-specific DPD model was found, ' 'and dpdTrainRun.m is not available.']);
        end
    end

    if shouldTrain
        if exist('dpdTrainRun', 'file') ~= 2
            error('dpdCompareRun:MissingTrainingFunction', 'trainBeforeCompare requires dpdTrainRun.m.');
        end

        fprintf('DPD comparison: training %s learning DPD.\n', learningArchitecture(dpdCfg));
        [dpdCfg, trainingInfo] = dpdTrainRun(cfg);
        cfg.dpd = dpdCfg;
    elseif ~trainedCoefficientsAvailable && modelFileExists
        refreshedDpdCfg = dpdConfig(cfg.pa, true, learningArchitecture(dpdCfg), cfrCfg);

        if hasTrainedCoefficients(refreshedDpdCfg)
            dpdCfg = refreshedDpdCfg;
            cfg.dpd = dpdCfg;
        else
            error('dpdCompareRun:ModelNotLoaded', ['The architecture-specific DPD model exists, but ' 'dpdConfig did not load trained coefficients from it.']);
        end
    end

    dpdCfg.enabled = true;

    comparisonSeed = getNestedOption( dpdCfg, 'comparison', 'seed', 19);
    [originalSystemSignal, txInfo] = generateComparisonSignal(cfg, comparisonSeed);
    originalSystemSignal = originalSystemSignal(:);

    [originalModelSignal, modelInputInfo] = paModelInput(originalSystemSignal, cfg);

    [nominalDpdModelSignal, nominalDpdInfo] = dpdCore(originalModelSignal, dpdCfg);

    validateFiniteSignal( originalModelSignal, 'Original model-rate signal');
    validateFiniteSignal( nominalDpdModelSignal, 'Nominal DPD output');

    paModelSignal = gmpCore( originalModelSignal, cfg.pa);

    equalOutputPower = getNestedOption( dpdCfg, 'comparison', 'equalOutputPower', true);

    if equalOutputPower
        [dpdDriveScale, powerMatchInfo] = matchModelOutputPower( originalModelSignal, nominalDpdModelSignal, cfg.pa, dpdCfg);
    else
        dpdDriveScale = 1;
        powerMatchInfo.enabled = false;
        powerMatchInfo.targetPower = mean(abs(paModelSignal).^2);
        powerMatchInfo.probePower = NaN;
        powerMatchInfo.driveScale = 1;
        powerMatchInfo.driveScaleDb = 0;
    end

    dpdModelSignal = nominalDpdModelSignal * dpdDriveScale;
    cascadeModelSignal = gmpCore(dpdModelSignal, cfg.pa);

    [dpdSystemSignal, dpdSystemInfo] = paSystemOutput( dpdModelSignal, length(originalSystemSignal), modelInputInfo.inputScale, modelInputInfo.rateFilter, cfg);

    [paSystemSignal, paSystemInfo] = paSystemOutput( paModelSignal, length(originalSystemSignal), modelInputInfo.inputScale, modelInputInfo.rateFilter, cfg);

    [cascadeSystemSignal, cascadeSystemInfo] = paSystemOutput( cascadeModelSignal, length(originalSystemSignal), modelInputInfo.inputScale, modelInputInfo.rateFilter, cfg);

    outputPowerToleranceDb = getNestedOption( dpdCfg, 'comparison', 'outputPowerToleranceDb', 0.05);

    outputPowerMismatchDb = powerRatioDb( mean(abs(cascadeSystemSignal).^2), mean(abs(paSystemSignal).^2));

    if equalOutputPower && abs(outputPowerMismatchDb) > outputPowerToleranceDb
        correctedDriveScale = dpdDriveScale * 10^(-outputPowerMismatchDb / 20);

        dpdModelSignal = nominalDpdModelSignal * correctedDriveScale;
        cascadeModelSignal = gmpCore(dpdModelSignal, cfg.pa);

        [dpdSystemSignal, dpdSystemInfo] = paSystemOutput( dpdModelSignal, length(originalSystemSignal), modelInputInfo.inputScale, modelInputInfo.rateFilter, cfg);

        [cascadeSystemSignal, cascadeSystemInfo] = paSystemOutput( cascadeModelSignal, length(originalSystemSignal), modelInputInfo.inputScale, modelInputInfo.rateFilter, cfg);

        dpdDriveScale = correctedDriveScale;
        outputPowerMismatchDb = powerRatioDb( mean(abs(cascadeSystemSignal).^2), mean(abs(paSystemSignal).^2));
    end

    powerMatchInfo.driveScale = dpdDriveScale;
    powerMatchInfo.driveScaleDb = 20 * log10(max(abs(dpdDriveScale), realmin));
    powerMatchInfo.systemOutputMismatchDb = outputPowerMismatchDb;
    powerMatchInfo.toleranceDb = outputPowerToleranceDb;
    powerMatchInfo.achieved = abs(outputPowerMismatchDb) <= outputPowerToleranceDb;

    validateFiniteSignal(dpdModelSignal, 'Actual DPD output');
    validateFiniteSignal(paModelSignal, 'PA output');
    validateFiniteSignal( cascadeModelSignal, 'DPD + PA output');

    dpdInfo = actualDpdInfo( nominalDpdInfo, originalModelSignal, dpdModelSignal, dpdCfg, dpdDriveScale);

    paMaximumDelay = max([ cfg.pa.signalDelays(:); cfg.pa.envelopeDelays(:)]);
    dpdMaximumDelay = max([ dpdCfg.signalDelays(:); dpdCfg.envelopeDelays(:)]);
    cascadeTransientLength = dpdMaximumDelay + paMaximumDelay;

    if length(originalModelSignal) <= cascadeTransientLength
        error('dpdCompareRun:InsufficientModelSamples', ['The model-rate signal must be longer than the combined ' 'DPD and PA memory.']);
    end

    validRows = cascadeTransientLength + 1: length(originalModelSignal);

    validOriginal = originalModelSignal(validRows);
    validPa = paModelSignal(validRows);
    validCascade = cascadeModelSignal(validRows);

    commonPaGain = bestLinearGain(validOriginal, validPa);

    if abs(commonPaGain) <= sqrt(realmin)
        commonPaGain = 1;
    end

    normalizedPaModelSignal = paModelSignal / commonPaGain;
    normalizedCascadeModelSignal = cascadeModelSignal / commonPaGain;

    [paNmseDb, paBestGain] = cascadeNmseDb( validOriginal, validPa);
    [cascadeNmseDbValue, cascadeBestGain] = cascadeNmseDb(validOriginal, validCascade);

    commonGainPaNmseDb = fixedGainNmseDb( validOriginal, validPa, commonPaGain);
    commonGainCascadeNmseDb = fixedGainNmseDb( validOriginal, validCascade, commonPaGain);

    spectrumMaximumSamples = round(getNestedOption( dpdCfg, 'comparison', 'spectrumMaximumSamples', 1048576));
    spectrumSegmentLength = round(getNestedOption( dpdCfg, 'comparison', 'spectrumSegmentLength', 8192));

    if spectrumSegmentLength < 2
        error('dpdCompareRun:InvalidWelchSegmentLength', 'spectrumSegmentLength must be at least two.');
    end

    [spectrumOriginal, spectrumRows] = selectAnalysisBlock( originalModelSignal, cascadeTransientLength, spectrumMaximumSamples);
    spectrumDpd = dpdModelSignal(spectrumRows);
    spectrumPa = normalizedPaModelSignal(spectrumRows);
    spectrumCascade = normalizedCascadeModelSignal(spectrumRows);

    [frequency, originalPsd, welchInfo] = welch( spectrumOriginal, cfg.pa.modelSampleRate, spectrumSegmentLength);
    [~, dpdPsd] = welch( spectrumDpd, cfg.pa.modelSampleRate, spectrumSegmentLength);
    [~, paPsd] = welch( spectrumPa, cfg.pa.modelSampleRate, spectrumSegmentLength);
    [~, cascadePsd] = welch( spectrumCascade, cfg.pa.modelSampleRate, spectrumSegmentLength);

    psdReference = max(originalPsd);

    if ~isfinite(psdReference) || psdReference <= 0
        error('dpdCompareRun:InvalidPsdReference', 'The Original PSD reference must be finite and positive.');
    end

    psdDb = 10 * log10([ originalPsd, dpdPsd, paPsd, cascadePsd] / psdReference + eps);

    originalAclr = aclr( frequency, originalPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    dpdAclr = aclr( frequency, dpdPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    paAclr = aclr( frequency, paPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    cascadeAclr = aclr( frequency, cascadePsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);

    spectrumFigure = plotSpectrumComparison( frequency, psdDb, cfg, dpdCfg);

    amAmMaximumSamples = round(getNestedOption( dpdCfg, 'comparison', 'amAmMaximumSamples', 200000));
    amAmBinCount = round(getNestedOption( dpdCfg, 'comparison', 'amAmBinCount', 180));

    amAmRows = evenlySpacedRows( cascadeTransientLength + 1, length(originalModelSignal), amAmMaximumSamples);
    amAmInput = originalModelSignal(amAmRows);
    amAmOutputs = [ originalModelSignal(amAmRows), dpdModelSignal(amAmRows), normalizedPaModelSignal(amAmRows), normalizedCascadeModelSignal(amAmRows)];

    [amAmInputAmplitude, amAmOutputAmplitude] = averageAmAm(amAmInput, amAmOutputs, amAmBinCount);

    amAmFigure = plotAmAmComparison( amAmInputAmplitude, amAmOutputAmplitude);

    averagePower = zeros(4, 1);
    peakPower = zeros(4, 1);
    [averagePower(1), peakPower(1)] = signalPower( originalModelSignal(validRows));
    [averagePower(2), peakPower(2)] = signalPower( dpdModelSignal(validRows));
    [averagePower(3), peakPower(3)] = signalPower( paModelSignal(validRows));
    [averagePower(4), peakPower(4)] = signalPower( cascadeModelSignal(validRows));
    paprDb = 10 * log10(peakPower ./ averagePower);

    signalName = {'Original'; 'DPD'; 'PA'; 'DPD+PA'};
    aclr1Db = [ originalAclr.worstAclr1Db; dpdAclr.worstAclr1Db; paAclr.worstAclr1Db; cascadeAclr.worstAclr1Db];
    aclr2Db = [ originalAclr.worstAclr2Db; dpdAclr.worstAclr2Db; paAclr.worstAclr2Db; cascadeAclr.worstAclr2Db];

    summary = table( signalName, averagePower, peakPower, paprDb, aclr1Db, aclr2Db, 'VariableNames', { 'Signal', 'AveragePower', 'PeakPower', 'PAPR_dB', 'ACLR1_dB', 'ACLR2_dB'});

    fprintf('\nClean DPD comparison before AWGN:\n');
    disp(summary);
    fprintf('PA cascade NMSE:       %.3f dB.\n', paNmseDb);
    fprintf('DPD + PA cascade NMSE: %.3f dB.\n', cascadeNmseDbValue);
    fprintf('DPD drive scale:       %.3f dB.\n', powerMatchInfo.driveScaleDb);
    fprintf('PA output mismatch:    %.3f dB.\n', outputPowerMismatchDb);

    results.cfg = cfg;
    results.dpdCfg = dpdCfg;
    results.trainingInfo = trainingInfo;
    results.txInfo = txInfo;
    results.dpdInfo = dpdInfo;
    results.modelInputInfo = modelInputInfo;
    results.dpdSystemInfo = dpdSystemInfo;
    results.paSystemInfo = paSystemInfo;
    results.cascadeSystemInfo = cascadeSystemInfo;
    results.powerMatch = powerMatchInfo;

    results.originalModelSignal = originalModelSignal;
    results.dpdModelSignal = dpdModelSignal;
    results.paModelSignal = paModelSignal;
    results.cascadeModelSignal = cascadeModelSignal;
    results.originalSystemSignal = originalSystemSignal;
    results.dpdSystemSignal = dpdSystemSignal;
    results.paSystemSignal = paSystemSignal;
    results.cascadeSystemSignal = cascadeSystemSignal;

    results.commonPaGain = commonPaGain;
    results.paBestLinearGain = paBestGain;
    results.cascadeBestLinearGain = cascadeBestGain;
    results.paNmseDb = paNmseDb;
    results.cascadeNmseDb = cascadeNmseDbValue;
    results.commonGainPaNmseDb = commonGainPaNmseDb;
    results.commonGainCascadeNmseDb = commonGainCascadeNmseDb;

    results.frequency = frequency;
    results.originalPsd = originalPsd;
    results.dpdPsd = dpdPsd;
    results.paPsd = paPsd;
    results.cascadePsd = cascadePsd;
    results.psdDb = psdDb;
    results.psdReference = psdReference;
    results.welchInfo = welchInfo;
    results.spectrumRows = spectrumRows;
    results.originalAclr = originalAclr;
    results.dpdAclr = dpdAclr;
    results.paAclr = paAclr;
    results.cascadeAclr = cascadeAclr;

    results.amAmInputAmplitude = amAmInputAmplitude;
    results.amAmOutputAmplitude = amAmOutputAmplitude;
    results.summary = summary;
    results.spectrumFigure = spectrumFigure;
    results.amAmFigure = amAmFigure;
end

function [txSignal, txInfo] = generateComparisonSignal(cfg, seed)
    [txSignal, txInfo] = tx(cfg, seed);
end

function available = hasTrainedCoefficients(dpdCfg)
    available = false;

    if isfield(dpdCfg, 'modelStale') && logical(dpdCfg.modelStale)
        return;
    end

    if ~isfield(dpdCfg, 'coefficients') || isempty(dpdCfg.coefficients)
        return;
    end

    coefficients = dpdCfg.coefficients;
    if isfield(dpdCfg, 'diagonalCount')
        activeMask = logical(gmpDiagonalMask( dpdCfg, dpdCfg.diagonalCount));
    elseif isfield(dpdCfg, 'coefficientMask') && isequal(size(dpdCfg.coefficientMask), size(coefficients))
        activeMask = logical(dpdCfg.coefficientMask);
    else
        activeMask = true(size(coefficients));
    end

    if ~any(activeMask(:))
        return;
    end

    identityCoefficients = complex(zeros(size(coefficients)));
    if activeMask(1, 1)
        identityCoefficients(1, 1) = 1;
    end

    activeCoefficients = coefficients(activeMask);
    activeIdentity = identityCoefficients(activeMask);
    tolerance = 1e-12 * max(1, norm(activeCoefficients));

    available = norm( activeCoefficients - activeIdentity) > tolerance;
end

function [driveScale, info] = matchModelOutputPower( originalSignal, dpdSignal, paCfg, dpdCfg)
    maximumProbeSamples = round(getNestedOption( dpdCfg, 'comparison', 'powerMatchProbeSamples', 65536));

    if maximumProbeSamples < 2
        error('dpdCompareRun:InvalidPowerMatchProbeLength', 'powerMatchProbeSamples must be at least two.');
    end

    if length(originalSignal) ~= length(dpdSignal)
        error('dpdCompareRun:PowerMatchLengthMismatch', 'Original and DPD model-rate signals must have equal length.');
    end

    maximumDelay = max([ paCfg.signalDelays(:); paCfg.envelopeDelays(:)]);
    signalLength = length(originalSignal);

    if signalLength <= maximumDelay
        error('dpdCompareRun:InsufficientPowerMatchSamples', 'The power-match signal is shorter than the PA memory.');
    end

    probeLength = min(maximumProbeSamples, signalLength);

    if probeLength <= maximumDelay
        error('dpdCompareRun:InsufficientPowerMatchProbe', 'The power-match probe is shorter than the PA memory.');
    end

    firstProbeSample = max(1, floor((signalLength - probeLength) / 2) + 1);
    probeRows = firstProbeSample: firstProbeSample + probeLength - 1;

    originalProbe = originalSignal(probeRows);
    dpdProbe = dpdSignal(probeRows);
    validRows = maximumDelay + 1:probeLength;

    baselineProbeOutput = gmpCore(originalProbe, paCfg);
    targetPower = mean(abs( baselineProbeOutput(validRows)).^2);

    if ~isfinite(targetPower) || targetPower <= 0
        error('dpdCompareRun:InvalidTargetOutputPower', 'The PA target output power must be finite and positive.');
    end

    lowerDb = getNestedOption( dpdCfg, 'comparison', 'minimumDriveScaleDb', -6);
    upperDb = getNestedOption( dpdCfg, 'comparison', 'maximumDriveScaleDb', 6);
    iterationCount = round(getNestedOption( dpdCfg, 'comparison', 'powerMatchIterations', 14));

    lowerScale = 10^(lowerDb / 20);
    upperScale = 10^(upperDb / 20);
    lowerPower = probeOutputPower( dpdProbe, lowerScale, paCfg, validRows);
    upperPower = probeOutputPower( dpdProbe, upperScale, paCfg, validRows);

    if ~isfinite(lowerPower) || ~isfinite(upperPower) || lowerPower <= 0 || upperPower <= 0
        error('dpdCompareRun:InvalidPowerMatchResponse', 'The PA power-match response must be finite and positive.');
    end

    if upperPower < lowerPower
        error('dpdCompareRun:NonmonotonicPowerMatchResponse', ['PA output power is not monotonic over the configured ' 'drive-scale interval.']);
    end

    if targetPower <= lowerPower
        driveScale = lowerScale;
    elseif targetPower >= upperPower
        driveScale = upperScale;
    else
        for iterationIndex = 1:iterationCount
            middleScale = sqrt(lowerScale * upperScale);
            middlePower = probeOutputPower( dpdProbe, middleScale, paCfg, validRows);

            if middlePower < targetPower
                lowerScale = middleScale;
            else
                upperScale = middleScale;
            end
        end

        driveScale = sqrt(lowerScale * upperScale);
    end

    matchedPower = probeOutputPower( dpdProbe, driveScale, paCfg, validRows);

    info.enabled = true;
    info.targetPower = targetPower;
    info.probePower = matchedPower;
    info.probeMismatchDb = powerRatioDb(matchedPower, targetPower);
    info.probeRows = probeRows;
    info.driveScale = driveScale;
    info.driveScaleDb = 20 * log10(driveScale);
end

function outputPower = probeOutputPower( dpdProbe, driveScale, paCfg, validRows)
    probeOutput = gmpCore(dpdProbe * driveScale, paCfg);
    outputPower = mean(abs(probeOutput(validRows)).^2);
end

function [analysisSignal, analysisRows] = selectAnalysisBlock(signal, maximumDelay, maximumSamples)
    if maximumSamples < 2
        error('dpdCompareRun:InvalidSpectrumSampleCount', 'spectrumMaximumSamples must be at least two.');
    end

    firstValidSample = maximumDelay + 1;
    availableLength = length(signal) - maximumDelay;

    if availableLength < 2
        error('dpdCompareRun:InsufficientSpectrumSamples', 'At least two valid samples are required for the spectrum.');
    end

    analysisLength = min(maximumSamples, availableLength);
    firstSample = firstValidSample + floor((availableLength - analysisLength) / 2);
    lastSample = firstSample + analysisLength - 1;
    analysisRows = (firstSample:lastSample).';
    analysisSignal = signal(analysisRows);
end

function rows = evenlySpacedRows( firstRow, lastRow, maximumRowCount)
    if lastRow < firstRow || maximumRowCount < 1
        error('dpdCompareRun:InvalidAmAmRows', 'AM-AM requires a nonempty valid sample interval.');
    end

    rowCount = lastRow - firstRow + 1;
    selectedCount = min(maximumRowCount, rowCount);
    rows = unique(round(linspace( firstRow, lastRow, selectedCount))).';
end

function [inputCurve, outputCurves] = averageAmAm(inputSignal, outputSignals, binCount)
    inputAmplitude = abs(inputSignal(:));
    outputAmplitude = abs(outputSignals);

    if isempty(inputAmplitude) || binCount < 1
        error('dpdCompareRun:InvalidAmAmInput', 'AM-AM requires at least one input sample and one bin.');
    end

    if size(outputAmplitude, 1) ~= length(inputAmplitude)
        error('dpdCompareRun:AmAmLengthMismatch', 'Every AM-AM output must match the input length.');
    end

    [sortedInput, sortIndices] = sort(inputAmplitude);
    sortedOutput = outputAmplitude(sortIndices, :);

    sampleCount = length(sortedInput);
    binCount = min(binCount, sampleCount);
    binEdges = round(linspace(1, sampleCount + 1, binCount + 1));

    inputCurve = zeros(binCount, 1);
    outputCurves = zeros(binCount, size(outputSignals, 2));

    for binIndex = 1:binCount
        firstSample = binEdges(binIndex);
        lastSample = max(firstSample, binEdges(binIndex + 1) - 1);
        rows = firstSample:lastSample;
        inputCurve(binIndex) = mean(sortedInput(rows));
        outputCurves(binIndex, :) = mean(sortedOutput(rows, :), 1);
    end
end

function spectrumFigure = plotSpectrumComparison( frequency, psdDb, cfg, dpdCfg)
    spectrumFigure = figure( 'Name', 'DPD spectrum comparison', 'NumberTitle', 'off', 'Color', 'w');
    hold on;

    colors = comparisonColors();
    lineHandles = gobjects(4, 1);
    lineHandles(1) = plot( frequency / 1e6, psdDb(:, 1), 'Color', colors(1, :), 'LineWidth', 1.4);
    lineHandles(2) = plot( frequency / 1e6, psdDb(:, 2), 'Color', colors(2, :), 'LineWidth', 1.3);
    lineHandles(3) = plot( frequency / 1e6, psdDb(:, 3), 'Color', colors(3, :), 'LineWidth', 1.3);
    lineHandles(4) = plot( frequency / 1e6, psdDb(:, 4), 'Color', colors(4, :), 'LineWidth', 1.5);

    requestedHalfSpanHz = getNestedOption( dpdCfg, 'comparison', 'spectrumHalfSpanHz', 20 * cfg.channelBandwidth);

    if ~isfinite(requestedHalfSpanHz) || requestedHalfSpanHz <= 0
        error('dpdCompareRun:InvalidSpectrumSpan', 'spectrumHalfSpanHz must be finite and positive.');
    end

    spectrumHalfSpanHz = min( cfg.pa.modelSampleRate / 2, requestedHalfSpanHz);
    grid on;
    xlim(spectrumHalfSpanHz / 1e6 * [-1 1]);
    ylim([-90 5]);
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    title('Clean spectra before AWGN');
    legend(lineHandles, {'Original', 'DPD', 'PA', 'DPD+PA'}, 'Location', 'best');
end

function amAmFigure = plotAmAmComparison( inputAmplitude, outputAmplitude)
    amAmFigure = figure( 'Name', 'DPD AM-AM comparison', 'NumberTitle', 'off', 'Color', 'w');
    hold on;

    colors = comparisonColors();
    lineHandles = gobjects(4, 1);
    lineHandles(1) = plot( inputAmplitude, outputAmplitude(:, 1), 'Color', colors(1, :), 'LineStyle', '--', 'LineWidth', 1.2);
    lineHandles(2) = plot( inputAmplitude, outputAmplitude(:, 2), 'Color', colors(2, :), 'LineWidth', 1.3);
    lineHandles(3) = plot( inputAmplitude, outputAmplitude(:, 3), 'Color', colors(3, :), 'LineWidth', 1.3);
    lineHandles(4) = plot( inputAmplitude, outputAmplitude(:, 4), 'Color', colors(4, :), 'LineWidth', 1.5);

    maximumAmplitude = max([ inputAmplitude; outputAmplitude(:)]);

    grid on;
    axis equal;
    xlim([0 maximumAmplitude]);
    ylim([0 maximumAmplitude]);
    xlabel('input amplitude');
    ylabel('output amplitude');
    title('AM-AM');
    legend(lineHandles, {'Original', 'DPD', 'PA', 'DPD+PA'}, 'Location', 'best');
end

function colors = comparisonColors()
    colors = [ 0.1500 0.1500 0.1500; 0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560; 0.1804 0.4902 0.1961];
end

function info = actualDpdInfo( nominalInfo, inputSignal, outputSignal, dpdCfg, driveScale)
    info = nominalInfo;
    info.nominalOutputAveragePower = nominalInfo.outputAveragePower;
    info.nominalOutputRms = nominalInfo.outputRms;
    info.nominalOutputPeakMagnitude = nominalInfo.outputPeakMagnitude;
    info.nominalOutputPaprDb = nominalInfo.outputPaprDb;

    [inputAveragePower, inputPeakPower] = signalPower(inputSignal);
    [outputAveragePower, outputPeakPower] = signalPower(outputSignal);

    info.inputAveragePower = inputAveragePower;
    info.outputAveragePower = outputAveragePower;
    info.inputRms = sqrt(inputAveragePower);
    info.outputRms = sqrt(outputAveragePower);
    info.outputPeakMagnitude = sqrt(outputPeakPower);
    info.inputPaprDb = 10 * log10( max(inputPeakPower / inputAveragePower, realmin));
    info.outputPaprDb = 10 * log10( max(outputPeakPower / outputAveragePower, realmin));
    info.appliedDriveScale = driveScale;
    info.appliedDriveScaleDb = 20 * log10(max(abs(driveScale), realmin));

    if isfield(dpdCfg, 'maximumOutputMagnitude') && isfinite(dpdCfg.maximumOutputMagnitude)
        info.outOfRangeFraction = mean( abs(outputSignal) > dpdCfg.maximumOutputMagnitude);
    else
        info.outOfRangeFraction = 0;
    end
end

function validateFiniteSignal(signal, signalName)
    if isempty(signal)
        error('dpdCompareRun:EmptySignal', '%s is empty.', signalName);
    end

    if any(~isfinite(signal))
        error('dpdCompareRun:NonfiniteSignal', '%s contains NaN or Inf.', signalName);
    end

    if mean(abs(signal).^2) <= 0
        error('dpdCompareRun:ZeroPowerSignal', '%s has zero average power.', signalName);
    end
end

function gain = bestLinearGain(reference, output)
    denominator = sum(abs(reference).^2);

    if denominator <= realmin
        gain = 1;
    else
        gain = sum(conj(reference) .* output) / denominator;
    end
end

function [nmseDb, gain] = cascadeNmseDb(reference, output)
    gain = bestLinearGain(reference, output);
    nmseDb = fixedGainNmseDb(reference, output, gain);
end

function nmseDb = fixedGainNmseDb(reference, output, gain)
    target = gain * reference;
    errorSignal = output - target;
    nmseRatio = sum(abs(errorSignal).^2) / max(sum(abs(target).^2), realmin);
    nmseDb = 10 * log10(max(nmseRatio, realmin));
end

function ratioDb = powerRatioDb(numerator, denominator)
    ratioDb = 10 * log10( max(numerator, realmin) / max(denominator, realmin));
end

function [averagePower, peakPower] = signalPower(signal)
    signalPowerValues = abs(signal).^2;
    averagePower = mean(signalPowerValues);
    peakPower = max(signalPowerValues);
end

function value = getNestedOption( parentStruct, childName, fieldName, defaultValue)
    value = defaultValue;

    if isfield(parentStruct, childName)
        childStruct = parentStruct.(childName);

        if isstruct(childStruct) && isfield(childStruct, fieldName)
            value = childStruct.(fieldName);
        end
    end
end

function architecture = learningArchitecture(dpdCfg)
    if isfield(dpdCfg, 'learningArchitecture') && ischar(dpdCfg.learningArchitecture)
        architecture = lower(strtrim(dpdCfg.learningArchitecture));
    else
        architecture = 'indirect';
    end

    if ~strcmp(architecture, 'indirect') && ~strcmp(architecture, 'direct')
        error('dpdCompareRun:UnknownLearningArchitecture', ['dpdCfg.learningArchitecture must be ''indirect'' ' 'or ''direct''.']);
    end
end

function matched = isDefaultArchitectureModelFile(modelFile)
    matched = false;
    if ~ischar(modelFile)
        return;
    end
    [~, modelName, modelExtension] = fileparts(modelFile);
    defaultNames = {'dpdModel', 'dpdModelIndirect', 'dpdModelDirect', 'dpdModelIndirectCfr', 'dpdModelDirectCfr'};
    matched = strcmpi(modelExtension, '.mat') && any(strcmpi(modelName, defaultNames));
end

function cfrCfg = systemCfrConfiguration(cfg)
    if isfield(cfg, 'cfr') && isstruct(cfg.cfr) && isscalar(cfg.cfr)
        cfrCfg = cfg.cfr;
    else
        cfrCfg = cfrConfig(false);
    end
end
