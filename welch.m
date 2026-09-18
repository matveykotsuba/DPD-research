function [frequency, psd, info] = welch( signal, sampleRate, segmentLength)
    signal = signal(:);

    if nargin < 3
        segmentLength = 8192;
    end

    segmentLength = min(segmentLength, length(signal));
    overlapLength = floor(3 * segmentLength / 4);
    step = segmentLength - overlapLength;
    fftLength = 2^nextpow2(2 * segmentLength);
    recordCount = floor((length(signal) - segmentLength) / step) + 1;

    n = (0:segmentLength-1).';
    window = 0.5 - 0.5 * cos(2 * pi * n / (segmentLength - 1));
    windowEnergy = sum(window.^2);

    psdSum = zeros(fftLength, 1);
    firstSample = 1;

    for recordIndex = 1:recordCount
        segment = signal(firstSample:firstSample+segmentLength-1);
        spectrum = fftshift(fft(segment .* window, fftLength));
        psdSum = psdSum + abs(spectrum).^2 / (sampleRate * windowEnergy);
        firstSample = firstSample + step;
    end

    psd = psdSum / recordCount;
    frequency = (-fftLength/2:fftLength/2-1).' * (sampleRate / fftLength);

    info.segmentLength = segmentLength;
    info.overlapLength = overlapLength;
    info.fftLength = fftLength;
    info.frequencyStep = sampleRate / fftLength;
    info.recordCount = recordCount;
end
