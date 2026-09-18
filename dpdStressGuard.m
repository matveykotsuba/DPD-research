function guardInfo = dpdStressGuard(candidateCfg, paCfg, stressInput)

    if nargin < 3 || ~isstruct(candidateCfg) || ~isstruct(paCfg)
        error('dpdStressGuard:MissingInput', 'candidateCfg, paCfg and stressInput are required.');
    end
    if ~isnumeric(stressInput) || ~isvector(stressInput) || isempty(stressInput) || any(~isfinite(stressInput(:)))
        error('dpdStressGuard:InvalidStressInput', 'stressInput must be a nonempty finite numeric vector.');
    end

    stressInput = stressInput(:);
    firstValidSample = maximumModelDelay(candidateCfg) + maximumModelDelay(paCfg) + 1;
    if firstValidSample > length(stressInput)
        error('dpdStressGuard:StressInputTooShort', 'stressInput is shorter than the combined DPD and PA memory.');
    end
    stressOutput = gmpCore(stressInput, candidateCfg);
    validRows = firstValidSample:length(stressInput);
    inputPeak = max(abs(stressInput(validRows)));
    outputPeak = max(abs(stressOutput(validRows)));
    maximumPeakGainDb = configurationValue(candidateCfg, 'maximumDpdPeakGainDb', 12);
    maximumOutputMagnitude = configurationValue(candidateCfg, 'maximumOutputMagnitude', Inf);
    validateattributes(maximumPeakGainDb, {'numeric'}, {'scalar', 'real'});
    validateattributes(maximumOutputMagnitude, {'numeric'}, {'scalar', 'real', 'positive'});
    relativePeakLimit = inputPeak * 10^(maximumPeakGainDb / 20);
    permittedOutputPeak = min(relativePeakLimit, maximumOutputMagnitude);
    guardInfo.inputPeakMagnitude = inputPeak;
    guardInfo.outputPeakMagnitude = outputPeak;
    guardInfo.maximumDpdPeakGainDb = maximumPeakGainDb;
    guardInfo.maximumOutputMagnitude = maximumOutputMagnitude;
    guardInfo.permittedOutputPeakMagnitude = permittedOutputPeak;
    guardInfo.isSafe = all(isfinite(stressOutput(validRows))) && isfinite(outputPeak) && outputPeak <= permittedOutputPeak;
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
