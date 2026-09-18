function [coefficients, info] = paRateFilter(cfg)
    factor = round(cfg.pa.rateFactor);
    neighbors = cfg.pa.resamplerNeighbors;
    attenuation = cfg.pa.resamplerAttenuationDb;

    filterOrder = 2 * neighbors * factor;
    coefficientIndex = (0:filterOrder).';
    centeredIndex = coefficientIndex - filterOrder / 2;

    cutoffFrequency = (cfg.pa.resamplerPassbandFrequency + cfg.pa.resamplerStopbandFrequency) / 2;

    normalizedCutoff = cutoffFrequency / cfg.pa.modelSampleRate;

    coefficients = sin( 2 * pi * normalizedCutoff * centeredIndex) ./ (pi * centeredIndex);

    coefficients(filterOrder / 2 + 1) = 2 * normalizedCutoff;

    if attenuation > 50
        beta = 0.1102 * (attenuation - 8.7);
    elseif attenuation >= 21
        beta = 0.5842 * (attenuation - 21)^0.4 + 0.07886 * (attenuation - 21);
    else
        beta = 0;
    end

    windowArgument = beta * sqrt(max(0, 1 - (2 * coefficientIndex / filterOrder - 1).^2));

    window = ones(filterOrder + 1, 1);
    windowTerm = ones(filterOrder + 1, 1);
    windowDenominator = 1;
    denominatorTerm = 1;

    for seriesIndex = 1:30
        windowTerm = windowTerm .* (windowArgument.^2 / 4) / seriesIndex^2;

        window = window + windowTerm;

        denominatorTerm = denominatorTerm * (beta^2 / 4) / seriesIndex^2;

        windowDenominator = windowDenominator + denominatorTerm;
    end

    coefficients = coefficients .*  window / windowDenominator;

    coefficients = coefficients / sum(coefficients);

    info.order = filterOrder;
    info.delaySamples = filterOrder / 2;
    info.cutoffFrequency = cutoffFrequency;
    info.passbandFrequency = cfg.pa.resamplerPassbandFrequency;
    info.stopbandFrequency = cfg.pa.resamplerStopbandFrequency;
end
