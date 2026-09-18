function metrics = dpdBlockMetrics(systemCfg, trainingInput, learningInfo, paCfg, dpdCfg, monitoringSignal)
    requiredLearningFields = {'coefficientHistory', 'blockHistory'};
    for fieldIndex = 1:length(requiredLearningFields)
        if ~isfield(learningInfo, requiredLearningFields{fieldIndex})
            error('dpdBlockMetrics:MissingLearningField', 'learningInfo must contain %s.', requiredLearningFields{fieldIndex});
        end
    end

    coefficientHistory = learningInfo.coefficientHistory;
    expectedCoefficientCount = numel(dpdCfg.coefficients);
    if size(coefficientHistory, 1) ~= expectedCoefficientCount || size(coefficientHistory, 2) < 2
        error('dpdBlockMetrics:InvalidCoefficientHistory', 'The coefficient history has an invalid size.');
    end
    if any(~isfinite(coefficientHistory(:)))
        error('dpdBlockMetrics:NonfiniteCoefficientHistory', 'The coefficient history must contain only finite values.');
    end

    trainingInput = trainingInput(:);
    combinedMemory = maximumModelDelay(dpdCfg) + maximumModelDelay(paCfg);
    validTrainingRows = combinedMemory + 1:length(trainingInput);
    if isempty(validTrainingRows)
        error('dpdBlockMetrics:InsufficientTrainingSamples', 'The training signal is shorter than the combined memory.');
    end
    trainingInputPower = mean(abs(trainingInput(validTrainingRows)).^2);
    if ~isfinite(trainingInputPower) || trainingInputPower <= realmin
        error('dpdBlockMetrics:InvalidTrainingPower', 'The training signal power must be finite and positive.');
    end
    chunkSampleCount = round(configurationValue(dpdCfg, 'adaptationMetricChunkSamples', 65536));
    validateattributes(chunkSampleCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    spectrumSegmentLength = round(configurationValue(dpdCfg, 'adaptationMetricSpectrumSegmentLength', 8192));
    validateattributes(spectrumSegmentLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    if spectrumSegmentLength < 2
        error('dpdBlockMetrics:InvalidWelchSegmentLength', 'adaptationMetricSpectrumSegmentLength must be at least two.');
    end

    monitorIsShared = nargin >= 6 && ~isempty(monitoringSignal);
    if monitorIsShared
        requiredMonitorFields = {'cfg', 'systemInput', 'txInfo', 'modelInput', 'inputInfo', 'info'};
        for fieldIndex = 1:length(requiredMonitorFields)
            if ~isfield(monitoringSignal, requiredMonitorFields{fieldIndex})
                error('dpdBlockMetrics:MissingMonitoringField', 'monitoringSignal must contain %s.', requiredMonitorFields{fieldIndex});
            end
        end
        monitorCfg = monitoringSignal.cfg;
        monitorSystemInput = monitoringSignal.systemInput;
        monitorTxInfo = monitoringSignal.txInfo;
        monitorModelInput = monitoringSignal.modelInput;
        monitorInputInfo = monitoringSignal.inputInfo;
        monitorSeed = monitoringSignal.info.seed;
        monitorSymbolCount = monitoringSignal.info.ofdmSymbolCount;
    else
        monitorCfg = systemCfg;
        monitorSymbolCount = round(configurationValue(dpdCfg, 'adaptationMetricOfdmSymbols', systemCfg.symbolsPerSlot));
        monitorSymbolCount = min(monitorSymbolCount, length(systemCfg.cpLengths));
        validateattributes(monitorSymbolCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
        monitorCfg.numOfdmSymbols = monitorSymbolCount;
        monitorCfg.cpLengths = systemCfg.cpLengths(1:monitorSymbolCount);
        monitorSeed = round(configurationValue(dpdCfg, 'adaptationMetricSeed', 14));
        validateattributes(monitorSeed, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
        [monitorSystemInput, monitorTxInfo] = tx(monitorCfg, monitorSeed);
        [monitorModelInput, monitorInputInfo] = paModelInput(monitorSystemInput, monitorCfg);
    end

    [prePaFrequency, prePaPsd, prePaWelchInfo] = welch(monitorSystemInput, monitorCfg.sampleRate, spectrumSegmentLength);
    prePaAclr = aclr(prePaFrequency, prePaPsd, monitorCfg.channelBandwidth, monitorCfg.aclrMeasurementBandwidth);

    snapshotCount = size(coefficientHistory, 2);
    updateCount = snapshotCount - 1;
    diagnosticUpdateIndex = unique([0 floor(updateCount / 2) updateCount]).';
    diagnosticSnapshotColumns = diagnosticUpdateIndex + 1;
    diagnosticSnapshotCount = length(diagnosticSnapshotColumns);
    evmRmsPercent = NaN(snapshotCount, 1);
    evmRmsDb = NaN(snapshotCount, 1);
    worstAclr1Db = NaN(snapshotCount, 1);
    worstAclr2Db = NaN(snapshotCount, 1);
    aclr1LeftDb = NaN(snapshotCount, 1);
    aclr1RightDb = NaN(snapshotCount, 1);
    aclr2LeftDb = NaN(snapshotCount, 1);
    aclr2RightDb = NaN(snapshotCount, 1);
    coefficientScale = NaN(snapshotCount, 1);
    coefficientMask = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);
    diagnosticFrequencyHz = [];
    diagnosticPsdColumns = cell(diagnosticSnapshotCount, 1);
    spectrumWelchInfo = struct();
    diagnosticAmAmReferenceGain = NaN;
    monitorFirstValidSample = combinedMemory + 1;
    if monitorFirstValidSample > length(monitorModelInput)
        error('dpdBlockMetrics:InsufficientMonitorSamples', 'The monitoring signal is shorter than the combined model memory.');
    end
    amAmMaximumSamples = round(nestedConfigurationValue(dpdCfg, 'comparison', 'amAmMaximumSamples', 200000));
    amAmBinCount = round(nestedConfigurationValue(dpdCfg, 'comparison', 'amAmBinCount', 180));
    validateattributes(amAmMaximumSamples, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(amAmBinCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    amAmRows = evenlySpacedRows(monitorFirstValidSample, length(monitorModelInput), amAmMaximumSamples);
    [diagnosticAmAmInputAmplitude, amAmSortIndices, amAmBinEdges] = prepareAmAmInput(monitorModelInput(amAmRows), amAmBinCount);
    diagnosticAmAmOutputAmplitude = NaN(length(diagnosticAmAmInputAmplitude), diagnosticSnapshotCount);
    coefficientHistoryNormalized = isfield(learningInfo, 'coefficientHistoryNormalized') && logical(learningInfo.coefficientHistoryNormalized);

    for snapshotIndex = 1:snapshotCount
        snapshotCfg = dpdCfg;
        snapshotCfg.enabled = true;
        snapshotCfg.coefficients = reshape(coefficientHistory(:, snapshotIndex), size(dpdCfg.coefficients));
        snapshotCfg.coefficients(~coefficientMask) = 0;
        if coefficientHistoryNormalized
            coefficientScale(snapshotIndex) = 1;
        else
            trainingDpdPower = streamedOutputPower(trainingInput, snapshotCfg, validTrainingRows(1), chunkSampleCount);
            if ~isfinite(trainingDpdPower) || trainingDpdPower <= realmin
                error('dpdBlockMetrics:InvalidDpdPower', 'DPD snapshot %d produced invalid training power.', snapshotIndex - 1);
            end
            coefficientScale(snapshotIndex) = sqrt(trainingInputPower / trainingDpdPower);
            snapshotCfg.coefficients = snapshotCfg.coefficients * coefficientScale(snapshotIndex);
        end
        monitorPaOutput = streamedCascade(monitorModelInput, snapshotCfg, paCfg, chunkSampleCount);
        diagnosticColumn = find(diagnosticSnapshotColumns == snapshotIndex, 1);
        if ~isempty(diagnosticColumn)
            diagnosticAmAmOutputAmplitude(:, diagnosticColumn) = averageAmAmOutput(monitorPaOutput(amAmRows), amAmSortIndices, amAmBinEdges);
            if snapshotIndex == 1
                diagnosticAmAmReferenceGain = streamedLinearGain(monitorModelInput, monitorPaOutput, monitorFirstValidSample, chunkSampleCount);
            end
        end
        monitorSystemOutput = paSystemOutput(monitorPaOutput, length(monitorSystemInput), monitorInputInfo.inputScale, monitorInputInfo.rateFilter, monitorCfg);
        clear monitorPaOutput
        [evmRmsPercent(snapshotIndex), evmRmsDb(snapshotIndex)] = qamEvm(monitorSystemOutput, monitorTxInfo.qamSymbols, monitorCfg);
        [frequency, outputPsd, currentWelchInfo] = welch(monitorSystemOutput, monitorCfg.sampleRate, spectrumSegmentLength);
        clear monitorSystemOutput
        if isempty(fieldnames(spectrumWelchInfo))
            spectrumWelchInfo = currentWelchInfo;
        end
        outputAclr = aclr(frequency, outputPsd, monitorCfg.channelBandwidth, monitorCfg.aclrMeasurementBandwidth);
        aclr1LeftDb(snapshotIndex) = outputAclr.aclr1LeftDb;
        aclr1RightDb(snapshotIndex) = outputAclr.aclr1RightDb;
        aclr2LeftDb(snapshotIndex) = outputAclr.aclr2LeftDb;
        aclr2RightDb(snapshotIndex) = outputAclr.aclr2RightDb;
        worstAclr1Db(snapshotIndex) = outputAclr.worstAclr1Db;
        worstAclr2Db(snapshotIndex) = outputAclr.worstAclr2Db;
        if ~isempty(diagnosticColumn)
            if isempty(diagnosticFrequencyHz)
                diagnosticFrequencyHz = frequency;
            elseif length(frequency) ~= length(diagnosticFrequencyHz) || any(frequency ~= diagnosticFrequencyHz)
                error('dpdBlockMetrics:InconsistentSpectrumGrid', 'All diagnostic spectra must use the same frequency grid.');
            end
            diagnosticPsdColumns{diagnosticColumn} = outputPsd;
        end
        clear outputPsd frequency
    end

    diagnosticPsd = horzcat(diagnosticPsdColumns{:});
    if any(~isfinite(diagnosticPsd(:))) || any(diagnosticPsd(:) < 0)
        error('dpdBlockMetrics:InvalidDiagnosticSpectrum', 'Diagnostic PSD values must be finite and nonnegative.');
    end
    diagnosticPsdReference = max(diagnosticPsd(:, 1));
    if ~isfinite(diagnosticPsdReference) || diagnosticPsdReference <= realmin
        error('dpdBlockMetrics:InvalidDiagnosticSpectrumReference', 'The update-zero spectrum reference must be finite and positive.');
    end
    diagnosticPsdDb = 10 * log10(diagnosticPsd / diagnosticPsdReference + eps);
    if ~isfinite(diagnosticAmAmReferenceGain) || abs(diagnosticAmAmReferenceGain) <= realmin
        error('dpdBlockMetrics:InvalidDiagnosticAmAmGain', 'The update-zero AM-AM reference gain must be finite and nonzero.');
    end
    diagnosticAmAmOutputAmplitude = diagnosticAmAmOutputAmplitude / abs(diagnosticAmAmReferenceGain);
    accepted = true(snapshotCount, 1);
    epochIndex = zeros(snapshotCount, 1);
    trainingBlockIndex = (0:updateCount).';
    if isfield(learningInfo.blockHistory, 'updateAccepted')
        if numel(learningInfo.blockHistory.updateAccepted) < updateCount
            error('dpdBlockMetrics:ShortAcceptanceHistory', 'The acceptance history is shorter than the coefficient history.');
        end
        accepted(2:end) = learningInfo.blockHistory.updateAccepted(1:updateCount);
    elseif isfield(learningInfo.blockHistory, 'accepted')
        if numel(learningInfo.blockHistory.accepted) < updateCount
            error('dpdBlockMetrics:ShortAcceptanceHistory', 'The acceptance history is shorter than the coefficient history.');
        end
        accepted(2:end) = learningInfo.blockHistory.accepted(1:updateCount);
    end
    if isfield(learningInfo.blockHistory, 'epochIndex')
        if numel(learningInfo.blockHistory.epochIndex) < updateCount
            error('dpdBlockMetrics:ShortEpochHistory', 'The epoch history is shorter than the coefficient history.');
        end
        epochIndex(2:end) = learningInfo.blockHistory.epochIndex(1:updateCount);
    end
    if isfield(learningInfo.blockHistory, 'blockIndex')
        if numel(learningInfo.blockHistory.blockIndex) < updateCount
            error('dpdBlockMetrics:ShortBlockHistory', 'The block history is shorter than the coefficient history.');
        end
        trainingBlockIndex(2:end) = learningInfo.blockHistory.blockIndex(1:updateCount);
    end
    deploymentHistorySelected = isfield(learningInfo, 'deploymentHistorySelected') && logical(learningInfo.deploymentHistorySelected);
    if deploymentHistorySelected
        verifyDeploymentSeries(evmRmsPercent, worstAclr1Db, worstAclr2Db, dpdCfg);
    end
    if isfield(learningInfo, 'deploymentTargetAclr1Db')
        targetAclr1Db = learningInfo.deploymentTargetAclr1Db;
    else
        targetAclr1Db = min(configurationValue(dpdCfg, 'adaptationTargetAclr1Db', 50), prePaAclr.worstAclr1Db);
    end
    validateattributes(targetAclr1Db, {'numeric'}, {'scalar', 'real', 'finite'});

    metrics.learningArchitecture = dpdCfg.learningArchitecture;
    metrics.learningInitialization = learningInitialization(learningInfo);
    metrics.historyStartsFromIdentity = historyStartsFromIdentity(learningInfo);
    metrics.updateIndex = (0:updateCount).';
    metrics.trainingBlockIndex = trainingBlockIndex;
    metrics.epochIndex = epochIndex;
    metrics.accepted = accepted;
    metrics.evmRmsPercent = evmRmsPercent;
    metrics.evmRmsDb = evmRmsDb;
    metrics.worstAclr1Db = worstAclr1Db;
    metrics.worstAclr2Db = worstAclr2Db;
    metrics.aclr1LeftDb = aclr1LeftDb;
    metrics.aclr1RightDb = aclr1RightDb;
    metrics.aclr2LeftDb = aclr2LeftDb;
    metrics.aclr2RightDb = aclr2RightDb;
    metrics.coefficientScale = coefficientScale;
    metrics.monitorSeed = monitorSeed;
    metrics.monitorOfdmSymbolCount = monitorSymbolCount;
    metrics.monitorSampleRate = monitorCfg.sampleRate;
    metrics.monitorIsSharedWithLearning = monitorIsShared;
    metrics.chunkSampleCount = chunkSampleCount;
    metrics.coefficientHistoryNormalized = coefficientHistoryNormalized;
    metrics.deploymentHistorySelected = deploymentHistorySelected;
    metrics.targetAclr1Db = targetAclr1Db;
    metrics.spectrumEstimator = 'Welch';
    metrics.spectrumWelchInfo = spectrumWelchInfo;
    metrics.prePaSpectrumWelchInfo = prePaWelchInfo;
    metrics.prePaWorstAclr1Db = prePaAclr.worstAclr1Db;
    metrics.prePaWorstAclr2Db = prePaAclr.worstAclr2Db;
    metrics.prePaAclr1LeftDb = prePaAclr.aclr1LeftDb;
    metrics.prePaAclr1RightDb = prePaAclr.aclr1RightDb;
    metrics.prePaAclr2LeftDb = prePaAclr.aclr2LeftDb;
    metrics.prePaAclr2RightDb = prePaAclr.aclr2RightDb;
    metrics.prePaAclr1Db = prePaAclr.worstAclr1Db;
    metrics.prePaAclr2Db = prePaAclr.worstAclr2Db;
    metrics.aclr1ShortfallDb = prePaAclr.worstAclr1Db - worstAclr1Db;
    metrics.aclr2ShortfallDb = prePaAclr.worstAclr2Db - worstAclr2Db;
    metrics.diagnosticUpdateIndex = diagnosticUpdateIndex;
    metrics.diagnosticTrainingBlockIndex = trainingBlockIndex(diagnosticSnapshotColumns);
    metrics.diagnosticFrequencyHz = diagnosticFrequencyHz;
    metrics.diagnosticPsd = diagnosticPsd;
    metrics.diagnosticPsdDb = diagnosticPsdDb;
    metrics.diagnosticPsdReference = diagnosticPsdReference;
    metrics.diagnosticAmAmInputAmplitude = diagnosticAmAmInputAmplitude;
    metrics.diagnosticAmAmOutputAmplitude = diagnosticAmAmOutputAmplitude;
    metrics.diagnosticAmAmReferenceGain = diagnosticAmAmReferenceGain;
    metrics.evmReference = 'QAM symbols after CP removal, FFT and best complex-gain correction';
    metrics.aclrReference = 'PA output before channel and AWGN on the fixed monitoring OFDM signal';
    metrics.spectrumReference = 'Welch PA system-output PSD on the fixed monitoring signal, with one common 0 dB reference equal to the update-zero peak PSD';
    metrics.amAmReference = 'Equal-population bin means of model-input magnitude and DPD-plus-PA output magnitude divided by the update-zero complex linear-gain magnitude';
    if isfield(learningInfo.blockHistory, 'deploymentSourceUpdateIndex')
        sourceUpdateIndex = learningInfo.blockHistory.deploymentSourceUpdateIndex(:);
        if length(sourceUpdateIndex) ~= updateCount
            error('dpdBlockMetrics:InvalidDeploymentSourceHistory', 'The deployment source history must contain one value per update.');
        end
        metrics.deploymentSourceUpdateIndex = [0; sourceUpdateIndex];
    end
end

function initialization = learningInitialization(learningInfo)
    initialization = 'unknown';
    if isfield(learningInfo, 'initialization') && ischar(learningInfo.initialization)
        initialization = lower(strtrim(learningInfo.initialization));
    elseif isfield(learningInfo, 'identityInitialization') && logical(learningInfo.identityInitialization)
        initialization = 'identity';
    end
end

function startsFromIdentity = historyStartsFromIdentity(learningInfo)
    startsFromIdentity = false;
    if ~isfield(learningInfo, 'coefficientHistory') || isempty(learningInfo.coefficientHistory)
        return
    end
    initialCoefficients = learningInfo.coefficientHistory(:, 1);
    identityCoefficients = complex(zeros(size(initialCoefficients)));
    identityCoefficients(1) = 1;
    startsFromIdentity = norm(initialCoefficients - identityCoefficients) <= 1e-10;
end

function verifyDeploymentSeries(evmRmsPercent, aclr1Db, aclr2Db, dpdCfg)
    maximumEvmIncreasePercent = configurationValue(dpdCfg, 'adaptationGateMaximumEvmIncreasePercent', 0);
    maximumAclr1DecreaseDb = configurationValue(dpdCfg, 'adaptationGateMaximumAclr1DecreaseDb', 0);
    maximumAclr2DecreaseDb = configurationValue(dpdCfg, 'adaptationGateMaximumAclr2DecreaseDb', 0.10);
    tolerance = configurationValue(dpdCfg, 'adaptationGateComparisonTolerance', 1e-9);
    if any(diff(evmRmsPercent) > maximumEvmIncreasePercent + tolerance)
        error('dpdBlockMetrics:EvmInvariantViolation', 'The deployment history increases RMS EVM on the fixed monitor.');
    end
    if any(diff(aclr1Db) < -maximumAclr1DecreaseDb - tolerance)
        error('dpdBlockMetrics:Aclr1InvariantViolation', 'The deployment history decreases ACLR1 on the fixed monitor.');
    end
    if any(diff(aclr2Db) < -maximumAclr2DecreaseDb - tolerance)
        error('dpdBlockMetrics:Aclr2InvariantViolation', 'The deployment history exceeds the ACLR2 tolerance on the fixed monitor.');
    end
end

function rows = evenlySpacedRows(firstRow, lastRow, maximumRowCount)
    if lastRow < firstRow || maximumRowCount < 1
        error('dpdBlockMetrics:InvalidAmAmRows', 'AM-AM requires a nonempty valid sample interval.');
    end
    rowCount = lastRow - firstRow + 1;
    selectedCount = min(maximumRowCount, rowCount);
    rows = unique(round(linspace(firstRow, lastRow, selectedCount))).';
end

function [inputCurve, sortIndices, binEdges] = prepareAmAmInput(inputSignal, binCount)
    inputAmplitude = abs(inputSignal(:));
    if isempty(inputAmplitude) || binCount < 1
        error('dpdBlockMetrics:InvalidAmAmInput', 'AM-AM requires at least one input sample and one bin.');
    end
    [sortedInput, sortIndices] = sort(inputAmplitude);
    sampleCount = length(sortedInput);
    binCount = min(binCount, sampleCount);
    binEdges = round(linspace(1, sampleCount + 1, binCount + 1));
    inputCurve = zeros(binCount, 1);
    for binIndex = 1:binCount
        firstSample = binEdges(binIndex);
        lastSample = max(firstSample, binEdges(binIndex + 1) - 1);
        inputCurve(binIndex) = mean(sortedInput(firstSample:lastSample));
    end
end

function outputCurve = averageAmAmOutput(outputSignal, sortIndices, binEdges)
    sortedOutputAmplitude = abs(outputSignal(sortIndices));
    binCount = length(binEdges) - 1;
    outputCurve = zeros(binCount, 1);
    for binIndex = 1:binCount
        firstSample = binEdges(binIndex);
        lastSample = max(firstSample, binEdges(binIndex + 1) - 1);
        outputCurve(binIndex) = mean(sortedOutputAmplitude(firstSample:lastSample));
    end
end

function outputPower = streamedOutputPower(inputSignal, modelCfg, firstValidSample, chunkSampleCount)
    state = [];
    powerSum = 0;
    powerSampleCount = 0;
    for firstSample = 1:chunkSampleCount:length(inputSignal)
        lastSample = min(firstSample + chunkSampleCount - 1, length(inputSignal));
        [outputBlock, state] = gmpCoreBlock(inputSignal(firstSample:lastSample), modelCfg, state);
        firstValidBlockSample = max(firstValidSample, firstSample);
        if firstValidBlockSample <= lastSample
            localFirstSample = firstValidBlockSample - firstSample + 1;
            validOutputBlock = outputBlock(localFirstSample:end);
            powerSum = powerSum + sum(abs(validOutputBlock).^2);
            powerSampleCount = powerSampleCount + length(validOutputBlock);
        end
    end
    outputPower = powerSum / powerSampleCount;
end

function linearGain = streamedLinearGain(inputSignal, outputSignal, firstValidSample, chunkSampleCount)
    numerator = complex(0);
    denominator = 0;
    for firstSample = firstValidSample:chunkSampleCount:length(inputSignal)
        lastSample = min(firstSample + chunkSampleCount - 1, length(inputSignal));
        inputBlock = inputSignal(firstSample:lastSample);
        outputBlock = outputSignal(firstSample:lastSample);
        numerator = numerator + inputBlock' * outputBlock;
        denominator = denominator + real(inputBlock' * inputBlock);
    end
    linearGain = numerator / denominator;
end

function paOutput = streamedCascade(inputSignal, dpdCfg, paCfg, chunkSampleCount)
    dpdState = [];
    paState = [];
    paOutput = complex(zeros(size(inputSignal)));
    for firstSample = 1:chunkSampleCount:length(inputSignal)
        lastSample = min(firstSample + chunkSampleCount - 1, length(inputSignal));
        [dpdOutputBlock, dpdState] = gmpCoreBlock(inputSignal(firstSample:lastSample), dpdCfg, dpdState);
        [paOutputBlock, paState] = gmpCoreBlock(dpdOutputBlock, paCfg, paState);
        paOutput(firstSample:lastSample) = paOutputBlock;
    end
end

function [evmRmsPercent, evmRmsDb] = qamEvm(systemOutput, referenceQamSymbols, cfg)
    symbolLengths = cfg.fftSize + cfg.cpLengths;
    requiredSampleCount = sum(symbolLengths);
    if length(systemOutput) < requiredSampleCount
        error('dpdBlockMetrics:ShortSystemOutput', 'The PA system output is shorter than the monitoring OFDM record.');
    end

    receivedNoCp = complex(zeros(cfg.fftSize, cfg.numOfdmSymbols));
    readIndex = 1;
    for symbolIndex = 1:cfg.numOfdmSymbols
        usefulStart = readIndex + cfg.cpLengths(symbolIndex);
        usefulStop = usefulStart + cfg.fftSize - 1;
        receivedNoCp(:, symbolIndex) = systemOutput(usefulStart:usefulStop);
        readIndex = usefulStop + 1;
    end

    receivedGrid = fftshift(fft(receivedNoCp, [], 1) * sqrt(cfg.numActiveSubcarriers) / cfg.fftSize, 1);
    receivedQamMatrix = receivedGrid(cfg.activeIndices, :);
    receivedQamSymbols = receivedQamMatrix(:);
    referenceQamSymbols = referenceQamSymbols(:);
    referenceEnergy = sum(abs(referenceQamSymbols).^2);
    bestLinearGain = sum(conj(referenceQamSymbols) .* receivedQamSymbols) / referenceEnergy;
    correctedQamSymbols = receivedQamSymbols / bestLinearGain;
    errorEnergy = sum(abs(correctedQamSymbols - referenceQamSymbols).^2);
    evmRms = sqrt(errorEnergy / referenceEnergy);
    evmRmsPercent = 100 * evmRms;
    evmRmsDb = 20 * log10(max(evmRms, realmin));
end

function maximumDelay = maximumModelDelay(modelCfg)
    maximumDelay = max([modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end

function value = nestedConfigurationValue(configuration, structureName, fieldName, defaultValue)
    if isfield(configuration, structureName) && isstruct(configuration.(structureName)) && isfield(configuration.(structureName), fieldName)
        value = configuration.(structureName).(fieldName);
    else
        value = defaultValue;
    end
end
