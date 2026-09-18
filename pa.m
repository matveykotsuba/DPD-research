function [paSignal, paInfo] = pa(inputSignal, cfg)
    inputSignal = inputSignal(:);
    paCfg = cfg.pa;
    if ~isfield(cfg, 'cfr') || isempty(cfg.cfr)
        cfg.cfr = cfrConfig(false);
    end
    if isfield(cfg, 'dpd') && isstruct(cfg.dpd)
        [cfg.dpd, cfrCompatible] = dpdRefreshCfrCompatibility(cfg.dpd, cfg.cfr);
        if isfield(cfg.dpd, 'enabled') && cfg.dpd.enabled && ~cfrCompatible
            error('pa:IncompatibleCfrDpd', 'The active DPD coefficients do not match the current CFR configuration. Retrain DPD before running the PA cascade.');
        end
    end
    dpdEnabled = isfield(cfg, 'dpd') && isfield(cfg.dpd, 'enabled') && cfg.dpd.enabled;

    inputAveragePower = mean(abs(inputSignal).^2);
    inputPeakPower = max(abs(inputSignal).^2);
    inputRms = sqrt(inputAveragePower);

    if ~paCfg.enabled && ~dpdEnabled
        paSignal = inputSignal;
        numberOfPoints = min(24000, length(inputSignal));
        pointIndices = round(linspace(1, length(inputSignal), numberOfPoints)).';
        paInfo.enabled = false;
        paInfo.dpdEnabled = false;
        paInfo.dpd = struct('enabled', false);
        paInfo.inputScale = 1;
        paInfo.inputAveragePower = inputAveragePower;
        paInfo.outputAveragePower = inputAveragePower;
        paInfo.modelOutputAveragePower = inputAveragePower;
        paInfo.averagePowerGainDb = 0;
        paInfo.bestLinearGain = 1;
        paInfo.nonlinearEvmPercent = 0;
        paInfo.inputPaprDb = 10 * log10(inputPeakPower / inputAveragePower);
        paInfo.outputPaprDb = paInfo.inputPaprDb;
        paInfo.amInput = inputSignal(pointIndices);
        paInfo.amOutput = inputSignal(pointIndices);
        paInfo.amReferenceInput = inputSignal(pointIndices);
        paInfo.amDpdOutput = inputSignal(pointIndices);
        paInfo.amPaOutput = inputSignal(pointIndices);
        paInfo.amDpdPaOutput = inputSignal(pointIndices);
        paInfo.paInputRms = inputRms;
        paInfo.systemOutput = struct('outputLength', length(inputSignal), 'inputScale', 1, 'modelOutputAveragePower', inputAveragePower, 'systemOutputAveragePower', inputAveragePower);
        return;
    end

    [modelInput, modelInputInfo] = paModelInput(inputSignal, cfg);
    rateFilter = modelInputInfo.rateFilter;
    rateInfo = modelInputInfo.rateInfo;
    interpolatedInputRms = modelInputInfo.interpolatedInputRms;
    modelInputRms = modelInputInfo.modelInputRms;
    inputScale = modelInputInfo.inputScale;

    if isfield(cfg, 'dpd')
        [paInput, dpdInfo] = dpdCore(modelInput, cfg.dpd);
    else
        paInput = modelInput;
        dpdInfo = struct('enabled', false, 'bypassed', true);
    end

    if paCfg.enabled
        modelOutput = gmpCore(paInput, paCfg);
    else
        modelOutput = paInput;
    end

    modelOutputAveragePower = mean(abs(modelOutput).^2);
    [paSignal, systemOutputInfo] = paSystemOutput(modelOutput, length(inputSignal), inputScale, rateFilter, cfg);
    outputAveragePower = mean(abs(paSignal).^2);
    outputPeakPower = max(abs(paSignal).^2);
    bestLinearGain = sum(conj(inputSignal) .* paSignal) / sum(abs(inputSignal).^2);
    nonlinearError = paSignal - bestLinearGain * inputSignal;
    nonlinearEvmPercent = 100 * sqrt(sum(abs(nonlinearError).^2) / sum(abs(bestLinearGain * inputSignal).^2));

    paMemoryLength = max([paCfg.signalDelays paCfg.envelopeDelays]);
    if dpdEnabled
        dpdMemoryLength = max([cfg.dpd.signalDelays cfg.dpd.envelopeDelays]);
    else
        dpdMemoryLength = 0;
    end
    firstStoredSample = paMemoryLength + dpdMemoryLength + 1;
    numberOfPoints = min(24000, length(modelInput) - firstStoredSample + 1);
    pointIndices = round(linspace(firstStoredSample, length(modelInput), numberOfPoints)).';

    paInfo.enabled = paCfg.enabled;
    paInfo.dpdEnabled = dpdEnabled;
    paInfo.dpd = dpdInfo;
    paInfo.orders = paCfg.orders;
    paInfo.signalDelays = paCfg.signalDelays;
    paInfo.envelopeDelays = paCfg.envelopeDelays;
    paInfo.coefficients = paCfg.coefficients;
    paInfo.inputBackoffDb = paCfg.inputBackoffDb;
    paInfo.referenceInputRms = paCfg.referenceInputRms;
    paInfo.modelSampleRate = paCfg.modelSampleRate;
    paInfo.rateFactor = paCfg.rateFactor;
    paInfo.rateFilter = rateInfo;
    paInfo.inputRms = inputRms;
    paInfo.interpolatedInputRms = interpolatedInputRms;
    paInfo.modelInputRms = modelInputRms;
    paInfo.modelInputPeakMagnitude = max(abs(modelInput));
    paInfo.modelInputPaprDb = 10 * log10(max(abs(modelInput)).^2 / mean(abs(modelInput).^2));
    paInfo.paInputRms = sqrt(mean(abs(paInput).^2));
    paInfo.paInputPeakMagnitude = max(abs(paInput));
    paInfo.paInputPaprDb = 10 * log10(max(abs(paInput)).^2 / mean(abs(paInput).^2));
    paInfo.inputScale = inputScale;
    paInfo.systemOutput = systemOutputInfo;
    paInfo.inputAveragePower = inputAveragePower;
    paInfo.outputAveragePower = outputAveragePower;
    paInfo.modelOutputAveragePower = modelOutputAveragePower;
    paInfo.outputPowerRelativeDb = 10 * log10(modelOutputAveragePower / paCfg.referenceInputRms^2);
    paInfo.averagePowerGainDb = 10 * log10(outputAveragePower / inputAveragePower);
    paInfo.bestLinearGain = bestLinearGain;
    paInfo.bestLinearGainDb = 20 * log10(abs(bestLinearGain));
    paInfo.nonlinearEvmPercent = nonlinearEvmPercent;
    paInfo.inputPaprDb = 10 * log10(inputPeakPower / inputAveragePower);
    paInfo.outputPaprDb = 10 * log10(outputPeakPower / outputAveragePower);
    paInfo.cfrSignature = cfrSignature(cfg.cfr);
    paInfo.amReferenceInput = modelInput(pointIndices);
    paInfo.amDpdOutput = paInput(pointIndices);
    if dpdEnabled
        if paCfg.enabled
            paInfo.amPaOutput = sampledGmpOutput(modelInput, paCfg, pointIndices);
        else
            paInfo.amPaOutput = modelInput(pointIndices);
        end
    else
        paInfo.amPaOutput = modelOutput(pointIndices);
    end
    paInfo.amDpdPaOutput = modelOutput(pointIndices);
    paInfo.amInput = paInput(pointIndices);
    paInfo.amOutput = modelOutput(pointIndices);
end

function outputSamples = sampledGmpOutput(inputSignal, modelCfg, sampleIndices)
    inputSignal = inputSignal(:);
    sampleIndices = sampleIndices(:);
    outputSamples = complex(zeros(size(sampleIndices)));
    signalDelays = modelCfg.signalDelays;
    envelopeDelays = modelCfg.envelopeDelays;
    orders = modelCfg.orders;
    hasDiagonalCount = isfield(modelCfg, 'diagonalCount');

    for signalIndex = 1:length(signalDelays)
        signalDelay = signalDelays(signalIndex);
        validRows = sampleIndices > signalDelay;
        delayedIndices = sampleIndices(validRows) - signalDelay;
        outputSamples(validRows) = outputSamples(validRows) + modelCfg.coefficients(signalIndex, 1) .* inputSignal(delayedIndices);
    end

    for orderIndex = 2:length(orders)
        exponent = orders(orderIndex) - 1;
        for envelopeIndex = 1:length(envelopeDelays)
            envelopeDelay = envelopeDelays(envelopeIndex);
            coefficientColumn = 2 + (orderIndex - 2) * length(envelopeDelays) + envelopeIndex - 1;
            for signalIndex = 1:length(signalDelays)
                signalDelay = signalDelays(signalIndex);
                if hasDiagonalCount && abs(signalDelay - envelopeDelay) > modelCfg.diagonalCount
                    continue
                end
                coefficient = modelCfg.coefficients(signalIndex, coefficientColumn);
                if coefficient == 0
                    continue
                end
                validRows = sampleIndices > max(signalDelay, envelopeDelay);
                signalIndices = sampleIndices(validRows) - signalDelay;
                envelopeIndices = sampleIndices(validRows) - envelopeDelay;
                signalPart = inputSignal(signalIndices);
                envelopePart = abs(inputSignal(envelopeIndices)).^exponent;
                outputSamples(validRows) = outputSamples(validRows) + coefficient .* signalPart .* envelopePart;
            end
        end
    end
end
