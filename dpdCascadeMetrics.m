function metrics = dpdCascadeMetrics( modelInput, dpdCfg, paCfg, systemCfg)

    if nargin < 3
        error('dpdCascadeMetrics:MissingInput', 'modelInput, dpdCfg and paCfg are required.');
    end

    if ~isnumeric(modelInput) || ~isvector(modelInput) || isempty(modelInput)
        error('dpdCascadeMetrics:InvalidInput', 'modelInput must be a nonempty numeric vector.');
    end

    modelInput = modelInput(:);

    [dpdOutput, dpdInfo] = dpdCore(modelInput, dpdCfg);
    paOutput = gmpCore(dpdOutput, paCfg);

    paMemory = maximumModelDelay(paCfg);
    if dpdCfg.enabled
        dpdMemory = maximumModelDelay(dpdCfg);
    else
        dpdMemory = 0;
    end

    firstValidSample = paMemory + dpdMemory + 1;
    if isfield(dpdCfg, 'metricFirstValidSample')
        configuredFirstValidSample = round(dpdCfg.metricFirstValidSample);
        validateattributes(configuredFirstValidSample, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
        firstValidSample = max( firstValidSample, configuredFirstValidSample);
    end
    if length(modelInput) < firstValidSample
        error('dpdCascadeMetrics:InsufficientSamples', 'The signal is shorter than the combined DPD and PA memory.');
    end

    validRows = firstValidSample:length(modelInput);
    validInput = modelInput(validRows);
    validDpdOutput = dpdOutput(validRows);
    validPaOutput = paOutput(validRows);

    inputFinite = all(isfinite(validInput));
    dpdFinite = all(isfinite(validDpdOutput));
    paFinite = all(isfinite(validPaOutput));

    inputPeak = finitePeak(validInput);
    dpdPeak = finitePeak(validDpdOutput);
    paPeak = finitePeak(validPaOutput);

    if inputFinite && paFinite
        inputEnergy = sum(abs(validInput).^2);
        if inputEnergy > realmin && isfinite(inputEnergy)
            bestLinearGain = (validInput' * validPaOutput) / inputEnergy;
            linearTarget = bestLinearGain * validInput;
            targetEnergy = sum(abs(linearTarget).^2);
            errorEnergy = sum(abs(validPaOutput - linearTarget).^2);

            if targetEnergy > realmin && isfinite(targetEnergy) && isfinite(errorEnergy)
                bestGainNmseDb = 10 * log10(max( errorEnergy / targetEnergy, realmin));
            else
                bestGainNmseDb = NaN;
            end
        else
            bestLinearGain = NaN;
            bestGainNmseDb = NaN;
        end
    else
        bestLinearGain = NaN;
        bestGainNmseDb = NaN;
    end

    referenceLinearGain = configurationValue( dpdCfg, 'metricReferenceLinearGain', bestLinearGain);
    if ~isscalar(referenceLinearGain) || ~isnumeric(referenceLinearGain) || ~isfinite(referenceLinearGain)
        error('dpdCascadeMetrics:InvalidReferenceLinearGain', 'metricReferenceLinearGain must be a finite numeric scalar.');
    end
    fixedGainNmseDb = linearNmseDb( validInput, validPaOutput, referenceLinearGain);
    nmseDb = fixedGainNmseDb;

    maximumPeakGainDb = configurationValue( dpdCfg, 'maximumDpdPeakGainDb', 12);
    validateattributes(maximumPeakGainDb, {'numeric'}, {'scalar', 'real'});

    if isfinite(maximumPeakGainDb)
        relativePeakLimit = inputPeak * 10^(maximumPeakGainDb / 20);
    else
        relativePeakLimit = Inf;
    end

    absolutePeakLimit = configurationValue( dpdCfg, 'maximumOutputMagnitude', Inf);
    validateattributes(absolutePeakLimit, {'numeric'}, {'scalar', 'real', 'positive'});

    permittedPeak = min(relativePeakLimit, absolutePeakLimit);
    peakAccepted = isfinite(dpdPeak) && dpdPeak <= permittedPeak;

    maximumCoefficientNorm = configurationValue( dpdCfg, 'maximumCoefficientNorm', Inf);
    coefficientAccepted = isfinite(dpdInfo.coefficientNorm) && dpdInfo.coefficientNorm <= maximumCoefficientNorm;

    nmseAccepted = isfinite(nmseDb) || (isinf(nmseDb) && nmseDb < 0);

    metrics.firstValidSample = firstValidSample;
    metrics.sampleCount = length(validRows);
    metrics.bestLinearGain = bestLinearGain;
    metrics.nmseDb = nmseDb;
    metrics.referenceLinearGain = referenceLinearGain;
    metrics.fixedGainNmseDb = fixedGainNmseDb;
    metrics.bestGainNmseDb = bestGainNmseDb;
    metrics.inputAveragePower = mean(abs(validInput).^2);
    metrics.dpdOutputAveragePower = mean(abs(validDpdOutput).^2);
    metrics.paOutputAveragePower = mean(abs(validPaOutput).^2);
    metrics.inputPeakMagnitude = inputPeak;
    metrics.dpdOutputPeakMagnitude = dpdPeak;
    metrics.paOutputPeakMagnitude = paPeak;
    metrics.dpdPeakGainDb = magnitudeRatioDb(dpdPeak, inputPeak);
    metrics.maximumDpdPeakGainDb = maximumPeakGainDb;
    metrics.maximumOutputMagnitude = absolutePeakLimit;
    metrics.permittedDpdPeakMagnitude = permittedPeak;
    metrics.inputFinite = inputFinite;
    metrics.dpdOutputFinite = dpdFinite;
    metrics.paOutputFinite = paFinite;
    metrics.peakAccepted = peakAccepted;
    metrics.coefficientNorm = dpdInfo.coefficientNorm;
    metrics.coefficientAccepted = coefficientAccepted;
    metrics.outOfRangeFraction = dpdInfo.outOfRangeFraction;
    metrics.isSafe = inputFinite && dpdFinite && paFinite && peakAccepted && coefficientAccepted && nmseAccepted;

    metrics.aclrAvailable = false;
    metrics.worstAclr1Db = NaN;
    metrics.worstAclr2Db = NaN;
    if nargin >= 4 && ~isempty(systemCfg) && paFinite
        requiredAclrFields = {'channelBandwidth', 'aclrMeasurementBandwidth'};
        aclrConfigurationAvailable = true;
        for fieldIndex = 1:length(requiredAclrFields)
            aclrConfigurationAvailable = aclrConfigurationAvailable && isfield(systemCfg, requiredAclrFields{fieldIndex});
        end

        if aclrConfigurationAvailable
            spectrumSegmentLength = configurationValue( dpdCfg, 'validationSpectrumSegmentLength', 8192);
            spectrumSegmentLength = min( round(spectrumSegmentLength), length(validPaOutput));
            validateattributes(spectrumSegmentLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});

            [frequency, outputPsd] = welch( validPaOutput, paCfg.modelSampleRate, spectrumSegmentLength);
            outputAclr = aclr(frequency, outputPsd, systemCfg.channelBandwidth, systemCfg.aclrMeasurementBandwidth);
            metrics.aclrAvailable = true;
            metrics.worstAclr1Db = outputAclr.worstAclr1Db;
            metrics.worstAclr2Db = outputAclr.worstAclr2Db;
        end
    end
end

function nmseDb = linearNmseDb(reference, output, linearGain)
    if any(~isfinite(reference)) || any(~isfinite(output)) || ~isfinite(linearGain)
        nmseDb = NaN;
        return;
    end

    target = linearGain * reference;
    targetEnergy = sum(abs(target).^2);
    errorEnergy = sum(abs(output - target).^2);
    if targetEnergy <= realmin || ~isfinite(targetEnergy) || ~isfinite(errorEnergy)
        nmseDb = NaN;
    else
        nmseDb = 10 * log10(max( errorEnergy / targetEnergy, realmin));
    end
end

function maximumDelay = maximumModelDelay(modelCfg)
    requiredFields = {'signalDelays', 'envelopeDelays'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(modelCfg, requiredFields{fieldIndex})
            error('dpdCascadeMetrics:MissingModelField', 'The model configuration must contain %s.', requiredFields{fieldIndex});
        end
    end

    maximumDelay = max([ modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end

function peakValue = finitePeak(signal)
    if isempty(signal) || any(~isfinite(signal))
        peakValue = Inf;
    else
        peakValue = max(abs(signal));
    end
end

function ratioDb = magnitudeRatioDb(numerator, denominator)
    if ~isfinite(numerator) || ~isfinite(denominator)
        ratioDb = Inf;
    elseif denominator <= realmin
        ratioDb = NaN;
    else
        ratioDb = 20 * log10(max(numerator, realmin) / denominator);
    end
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
