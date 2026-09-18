function comparison = dpdCurrentRunAnalysis(modelInput, dpdOutput, cascadeOutput, cfg)
    modelInput = modelInput(:);
    dpdOutput = dpdOutput(:);
    cascadeOutput = cascadeOutput(:);

    if length(dpdOutput) ~= length(modelInput) || length(cascadeOutput) ~= length(modelInput)
        error('dpdCurrentRunAnalysis:SignalLengthMismatch', 'All model-rate signals must have equal length.');
    end

    paCfg = cfg.pa;
    dpdCfg = cfg.dpd;
    paDelay = max([paCfg.signalDelays(:); paCfg.envelopeDelays(:)]);
    dpdDelay = max([dpdCfg.signalDelays(:); dpdCfg.envelopeDelays(:)]);
    combinedDelay = paDelay + dpdDelay;
    availableSampleCount = length(modelInput) - combinedDelay;

    if availableSampleCount < 2
        error('dpdCurrentRunAnalysis:InsufficientSamples', 'The current signal is shorter than the combined DPD and PA memory.');
    end

    maximumSpectrumSamples = round(getComparisonOption(dpdCfg, 'spectrumMaximumSamples', 1048576));
    analysisSampleCount = min(maximumSpectrumSamples, availableSampleCount);
    firstAnalysisSample = combinedDelay + 1 + floor((availableSampleCount - analysisSampleCount) / 2);
    lastAnalysisSample = firstAnalysisSample + analysisSampleCount - 1;
    analysisRows = firstAnalysisSample:lastAnalysisSample;

    originalSignal = modelInput(analysisRows);
    dpdSignal = dpdOutput(analysisRows);
    cascadeSignal = cascadeOutput(analysisRows);

    firstPaInputSample = firstAnalysisSample - paDelay;
    paInputBlock = modelInput(firstPaInputSample:lastAnalysisSample);
    paOutputBlock = gmpCore(paInputBlock, paCfg);
    firstPaOutputSample = paDelay + 1;
    paSignal = paOutputBlock(firstPaOutputSample:firstPaOutputSample + analysisSampleCount - 1);

    inputEnergy = sum(abs(originalSignal).^2);
    if inputEnergy <= realmin || ~isfinite(inputEnergy)
        error('dpdCurrentRunAnalysis:InvalidInputEnergy', 'The current signal must have finite positive energy.');
    end

    commonPaGain = sum(conj(originalSignal) .* paSignal) / inputEnergy;
    if ~isfinite(commonPaGain) || abs(commonPaGain) <= sqrt(realmin)
        commonPaGain = 1;
    end

    normalizedPaSignal = paSignal / commonPaGain;
    normalizedCascadeSignal = cascadeSignal / commonPaGain;
    signals = [originalSignal dpdSignal normalizedPaSignal normalizedCascadeSignal];

    spectrumSegmentLength = round(getComparisonOption(dpdCfg, 'spectrumSegmentLength', 8192));
    spectrumSegmentLength = min(spectrumSegmentLength, analysisSampleCount);
    [frequency, originalPsd] = welch(signals(:, 1), paCfg.modelSampleRate, spectrumSegmentLength);
    [~, dpdPsd] = welch(signals(:, 2), paCfg.modelSampleRate, spectrumSegmentLength);
    [~, paPsd] = welch(signals(:, 3), paCfg.modelSampleRate, spectrumSegmentLength);
    [~, cascadePsd] = welch(signals(:, 4), paCfg.modelSampleRate, spectrumSegmentLength);
    psd = [originalPsd dpdPsd paPsd cascadePsd];
    psdReference = max(originalPsd);

    if ~isfinite(psdReference) || psdReference <= 0
        error('dpdCurrentRunAnalysis:InvalidPsdReference', 'The Original PSD reference must be finite and positive.');
    end

    psdDb = 10 * log10(psd / psdReference + eps);
    aclrResults(1) = aclr(frequency, originalPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    aclrResults(2) = aclr(frequency, dpdPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    aclrResults(3) = aclr(frequency, paPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    aclrResults(4) = aclr(frequency, cascadePsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);

    signalNames = {'Original'; 'DPD'; 'PA'; 'DPD+PA'};
    aclr1LeftDb = [aclrResults.aclr1LeftDb].';
    aclr1RightDb = [aclrResults.aclr1RightDb].';
    worstAclr1Db = [aclrResults.worstAclr1Db].';
    aclr2LeftDb = [aclrResults.aclr2LeftDb].';
    aclr2RightDb = [aclrResults.aclr2RightDb].';
    worstAclr2Db = [aclrResults.worstAclr2Db].';
    aclrSummary = table(signalNames, aclr1LeftDb, aclr1RightDb, worstAclr1Db, aclr2LeftDb, aclr2RightDb, worstAclr2Db, 'VariableNames', {'Signal', 'ACLR1_L_dB', 'ACLR1_R_dB', 'ACLR1_Worst_dB', 'ACLR2_L_dB', 'ACLR2_R_dB', 'ACLR2_Worst_dB'});

    maximumAmSamples = round(getComparisonOption(dpdCfg, 'amAmMaximumSamples', 200000));
    selectedAmSampleCount = min(maximumAmSamples, analysisSampleCount);
    selectedAmRows = unique(round(linspace(1, analysisSampleCount, selectedAmSampleCount))).';
    selectedInput = originalSignal(selectedAmRows);
    selectedSignals = signals(selectedAmRows, :);
    [amInputAmplitude, amOutputAmplitude, amOutputPhaseDeg] = averageAmplitudeAndPhase(selectedInput, selectedSignals, round(getComparisonOption(dpdCfg, 'amAmBinCount', 180)));

    paError = normalizedPaSignal - originalSignal;
    cascadeError = normalizedCascadeSignal - originalSignal;
    paNmseDb = 10 * log10(max(sum(abs(paError).^2) / inputEnergy, realmin));
    cascadeNmseDb = 10 * log10(max(sum(abs(cascadeError).^2) / inputEnergy, realmin));
    outputPowerMismatchDb = 10 * log10(mean(abs(cascadeSignal).^2) / mean(abs(paSignal).^2));

    comparison.enabled = true;
    comparison.signalNames = signalNames;
    comparison.referencePlane = 'model rate before channel and AWGN';
    comparison.modelSampleRate = paCfg.modelSampleRate;
    comparison.modelSampleCount = length(modelInput);
    comparison.firstValidSample = combinedDelay + 1;
    comparison.firstAnalysisSample = firstAnalysisSample;
    comparison.lastAnalysisSample = lastAnalysisSample;
    comparison.analysisSampleCount = analysisSampleCount;
    comparison.commonPaGain = commonPaGain;
    comparison.commonPaGainDb = 20 * log10(abs(commonPaGain));
    comparison.paOutputPowerMismatchDb = outputPowerMismatchDb;
    comparison.paNmseDb = paNmseDb;
    comparison.cascadeNmseDb = cascadeNmseDb;
    comparison.frequency = frequency;
    comparison.psd = psd;
    comparison.psdDb = psdDb;
    comparison.psdReference = psdReference;
    comparison.aclr = aclrResults;
    comparison.aclrSummary = aclrSummary;
    comparison.amInputAmplitude = amInputAmplitude;
    comparison.amOutputAmplitude = amOutputAmplitude;
    comparison.amOutputPhaseDeg = amOutputPhaseDeg;
end

function value = getComparisonOption(dpdCfg, fieldName, defaultValue)
    value = defaultValue;
    if isfield(dpdCfg, 'comparison') && isstruct(dpdCfg.comparison) && isfield(dpdCfg.comparison, fieldName)
        value = dpdCfg.comparison.(fieldName);
    end
end

function [inputCurve, amplitudeCurves, phaseCurves] = averageAmplitudeAndPhase(inputSignal, outputSignals, binCount)
    inputSignal = inputSignal(:);
    inputAmplitude = abs(inputSignal);
    [sortedAmplitude, sortIndices] = sort(inputAmplitude);
    sortedInput = inputSignal(sortIndices);
    sortedOutputs = outputSignals(sortIndices, :);
    sampleCount = length(sortedInput);
    binCount = min(max(1, binCount), sampleCount);
    binEdges = round(linspace(1, sampleCount + 1, binCount + 1));
    inputCurve = zeros(binCount, 1);
    amplitudeCurves = zeros(binCount, size(outputSignals, 2));
    phaseCurves = zeros(binCount, size(outputSignals, 2));

    for binIndex = 1:binCount
        firstSample = binEdges(binIndex);
        lastSample = max(firstSample, binEdges(binIndex + 1) - 1);
        rows = firstSample:lastSample;
        inputCurve(binIndex) = mean(sortedAmplitude(rows));
        amplitudeCurves(binIndex, :) = mean(abs(sortedOutputs(rows, :)), 1);
        for signalIndex = 1:size(outputSignals, 2)
            phaseVector = sum(conj(sortedInput(rows)) .* sortedOutputs(rows, signalIndex));
            if abs(phaseVector) > realmin
                phaseCurves(binIndex, signalIndex) = angle(phaseVector) * 180 / pi;
            else
                phaseCurves(binIndex, signalIndex) = NaN;
            end
        end
    end
end
