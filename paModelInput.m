function [modelInput, info] = paModelInput(inputSignal, cfg)
    inputSignal = inputSignal(:);

    [rateFilter, rateInfo] = paRateFilter(cfg);

    modelRateInput = paInterpolate( inputSignal, cfg.pa.rateFactor, rateFilter);

    interpolatedInputRms = sqrt( mean(abs(modelRateInput).^2));

    modelInputRms = cfg.pa.referenceInputRms * 10^(-cfg.pa.inputBackoffDb / 20);

    if interpolatedInputRms > 0
        inputScale = modelInputRms / interpolatedInputRms;
    else
        inputScale = 1;
    end

    modelInput = modelRateInput * inputScale;

    info.rateFilter = rateFilter;
    info.rateInfo = rateInfo;
    info.inputScale = inputScale;
    info.interpolatedInputRms = interpolatedInputRms;
    info.modelInputRms = modelInputRms;
    info.modelSampleRate = cfg.pa.modelSampleRate;
    info.systemSampleRate = cfg.sampleRate;
    info.rateFactor = cfg.pa.rateFactor;
end
