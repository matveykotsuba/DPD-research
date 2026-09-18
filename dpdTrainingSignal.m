function [modelInput, signalInfo] = dpdTrainingSignal( cfg, seed, requestedSampleCount)

    if nargin < 3
        error('dpdTrainingSignal:MissingInput', 'cfg, seed and requestedSampleCount are required.');
    end

    validateattributes(seed, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
    validateattributes(requestedSampleCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});

    if ~isstruct(cfg) || ~isfield(cfg, 'pa')
        error('dpdTrainingSignal:MissingPaConfiguration', 'cfg.pa is required.');
    end

    rateFactor = round(cfg.pa.rateFactor);
    if rateFactor < 1 || abs(cfg.pa.rateFactor - rateFactor) > eps(max(1, cfg.pa.rateFactor)) * 16
        error('dpdTrainingSignal:NonintegerRateFactor', 'cfg.pa.rateFactor must be a positive integer.');
    end

    [rateFilter, rateInfo] = paRateFilter(cfg);
    edgeTrim = rateInfo.delaySamples;
    requiredSystemSamples = ceil( (requestedSampleCount + 2 * edgeTrim) / rateFactor) + 2;

    [fullTxSignal, txInfo] = tx(cfg, seed, false);
    fullTxSignal = fullTxSignal(:);
    if length(fullTxSignal) < requiredSystemSamples
        error('dpdTrainingSignal:InsufficientTxSamples', ['The TX frame has %d samples, but %d are required for ' 'the requested model-rate DPD record.'], length(fullTxSignal), requiredSystemSamples);
    end

    firstSystemSample = floor( (length(fullTxSignal) - requiredSystemSamples) / 2) + 1;
    lastSystemSample = firstSystemSample + requiredSystemSamples - 1;
    systemSegment = fullTxSignal(firstSystemSample:lastSystemSample);

    [uncroppedModelInput, modelInputInfo] = paModelInput(systemSegment, cfg);
    firstUsableSample = edgeTrim + 1;
    lastUsableSample = length(uncroppedModelInput) - edgeTrim;
    availableSampleCount = lastUsableSample - firstUsableSample + 1;

    if availableSampleCount < requestedSampleCount
        error('dpdTrainingSignal:InsufficientModelSamples', ['Rate conversion produced only %d usable model samples; ' '%d were requested.'], availableSampleCount, requestedSampleCount);
    end

    firstModelSample = firstUsableSample + floor( (availableSampleCount - requestedSampleCount) / 2);
    lastModelSample = firstModelSample + requestedSampleCount - 1;
    modelInput = uncroppedModelInput(firstModelSample:lastModelSample);

    actualRms = sqrt(mean(abs(modelInput).^2));
    targetRms = modelInputInfo.modelInputRms;
    if actualRms <= 0 || ~isfinite(actualRms)
        error('dpdTrainingSignal:InvalidModelInput', 'The generated model-rate signal has invalid RMS.');
    end
    postTrimScale = targetRms / actualRms;
    modelInput = modelInput * postTrimScale;

    signalInfo.referencePoint = 'after interpolation and PA-input scaling, before DPD';
    signalInfo.seed = seed;
    signalInfo.txSampleCount = length(fullTxSignal);
    signalInfo.systemSegmentFirstSample = firstSystemSample;
    signalInfo.systemSegmentLastSample = lastSystemSample;
    signalInfo.systemSegmentSampleCount = length(systemSegment);
    signalInfo.modelFirstSample = firstModelSample;
    signalInfo.modelLastSample = lastModelSample;
    signalInfo.modelSampleCount = length(modelInput);
    signalInfo.edgeTrimSamples = edgeTrim;
    signalInfo.modelSampleRate = cfg.pa.modelSampleRate;
    signalInfo.systemSampleRate = cfg.sampleRate;
    signalInfo.targetRms = targetRms;
    signalInfo.actualRms = sqrt(mean(abs(modelInput).^2));
    signalInfo.peakMagnitude = max(abs(modelInput));
    signalInfo.postTrimScale = postTrimScale;
    signalInfo.txPaprDb = txInfo.paprDb;
    signalInfo.rateFilterOrder = length(rateFilter) - 1;
end
