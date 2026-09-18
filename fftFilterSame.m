function filteredSignal = fftFilterSame(signal, coefficients)
    signal = signal(:);
    coefficients = coefficients(:);

    delaySamples = (length(coefficients) - 1 ) / 2;

    extendedSignal = [signal(end-delaySamples+1:end);
        signal;
        signal(1:delaySamples)
        ];

    convolutionLength = length(extendedSignal) + length(coefficients) - 1;
    ffLength = 2^nextpow2(convolutionLength);

    filteredExtendedSignal = ifft(fft(extendedSignal, ffLength) .* fft(coefficients, ffLength));

    firstSample = 2* delaySamples + 1;

    filteredSignal = filteredExtendedSignal(firstSample: firstSample + length(signal) - 1);
end
