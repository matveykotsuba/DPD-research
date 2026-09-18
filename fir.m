function [filteredSignal, info] = fir(signal, cfg)
    signal = signal(:);

    filterOrder = cfg.filterOrder;

    attenuation = cfg.firStopbandAttenuationDb;
    if attenuation > 50
        beta = 0.1102 * (attenuation - 8.7);
    elseif attenuation >= 21
        beta = 0.5842 * (attenuation - 21)^0.4 + 0.07886 * (attenuation - 21);
    else
        beta = 0;
    end

    cutoffFrequency = (cfg.firPassbandFrequency + cfg.firStopbandFrequency) / 2;
    coefficientIndex = (0:filterOrder).';
    centeredIndex = coefficientIndex - filterOrder / 2;

    coefficients = sin(2 * pi * cutoffFrequency / cfg.sampleRate * centeredIndex) ./ (pi * centeredIndex);
    coefficients(filterOrder / 2 + 1) = 2 * cutoffFrequency / cfg.sampleRate;

    windowArgument = beta * sqrt(1 - (2 * coefficientIndex / filterOrder - 1).^2);
    window = ones(filterOrder + 1, 1);
    windowTerm = ones(filterOrder + 1, 1);
    windowDenominator = 1;
    denominatorTerm = 1;

    for seriesIndex = 1:20
        windowTerm = windowTerm .* (windowArgument.^2 / 4) / seriesIndex^2;
        window = window + windowTerm;
        denominatorTerm = denominatorTerm * (beta^2 / 4) / seriesIndex^2;
        windowDenominator = windowDenominator + denominatorTerm;
    end

    coefficients = coefficients .* window / windowDenominator;
    coefficients = coefficients / sum(coefficients);

    delaySamples = filterOrder / 2;
    extendedSignal = [signal(end-delaySamples+1:end); signal; signal(1:delaySamples)];
    convolutionLength = length(extendedSignal) + length(coefficients) - 1;
    fftLength = 2^nextpow2(convolutionLength);
    filteredExtendedSignal = ifft(fft(extendedSignal, fftLength) .* fft(coefficients, fftLength));
    firstSample = 2 * delaySamples + 1;
    filteredSignal = filteredExtendedSignal(firstSample:firstSample+length(signal)-1);

    info.coefficients = coefficients;
    info.order = filterOrder;
    info.delaySamples = delaySamples;
    info.delaySeconds = delaySamples / cfg.sampleRate;
    info.passbandFrequency = cfg.firPassbandFrequency;
    info.stopbandFrequency = cfg.firStopbandFrequency;
    info.stopbandAttenuationDb = cfg.firStopbandAttenuationDb;
    info.cutoffFrequency = cutoffFrequency;
end
