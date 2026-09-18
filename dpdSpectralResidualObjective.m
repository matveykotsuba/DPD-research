function [objective, effectiveResidual, info] = dpdSpectralResidualObjective(residual, sampleRate, channelBandwidth, measurementBandwidth, rho1, rho2)
    if ~isnumeric(residual) || ~isvector(residual) || isempty(residual) || any(~isfinite(residual(:)))
        error('dpdSpectralResidualObjective:InvalidResidual', 'residual must be a finite nonempty numeric vector.');
    end
    validateattributes(sampleRate, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(channelBandwidth, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(measurementBandwidth, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(rho1, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(rho2, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});

    residualSize = size(residual);
    residual = residual(:);
    sampleCount = length(residual);
    if sampleCount < 2
        error('dpdSpectralResidualObjective:ResidualTooShort', 'residual must contain at least two samples.');
    end

    timeObjective = sum(abs(residual).^2);
    sampleIndex = (0:sampleCount-1).';
    window = 0.5 - 0.5 * cos(2 * pi * sampleIndex / (sampleCount - 1));
    windowEnergy = sum(window.^2);
    fftLength = 2^nextpow2(2 * sampleCount);
    spectrum = fftshift(fft(residual .* window, fftLength));
    frequencyStep = sampleRate / fftLength;
    frequency = (-fftLength/2:fftLength/2-1).' * frequencyStep;

    bandCenters = channelBandwidth * [-2 -1 0 1 2].';
    halfMeasurementBandwidth = measurementBandwidth / 2;
    bandEdges = [bandCenters-halfMeasurementBandwidth, bandCenters+halfMeasurementBandwidth];
    binLowerEdges = frequency - frequencyStep / 2;
    binUpperEdges = frequency + frequencyStep / 2;
    bandWeights = zeros(fftLength, length(bandCenters));
    for bandIndex = 1:length(bandCenters)
        overlapHz = max(0, min(binUpperEdges, bandEdges(bandIndex, 2)) - max(binLowerEdges, bandEdges(bandIndex, 1)));
        bandWeights(:, bandIndex) = overlapHz / frequencyStep;
    end

    aclr1FrequencyWeight = bandWeights(:, 2) + bandWeights(:, 4);
    aclr2FrequencyWeight = bandWeights(:, 1) + bandWeights(:, 5);
    spectralEnergyScale = sampleCount / (fftLength * windowEnergy);
    spectrumPower = abs(spectrum).^2;
    aclr1Energy = spectralEnergyScale * sum(aclr1FrequencyWeight .* spectrumPower);
    aclr2Energy = spectralEnergyScale * sum(aclr2FrequencyWeight .* spectrumPower);
    aclr1Inverse = ifft(ifftshift(aclr1FrequencyWeight .* spectrum));
    aclr2Inverse = ifft(ifftshift(aclr2FrequencyWeight .* spectrum));
    aclr1Adjoint = spectralEnergyScale * fftLength * window .* aclr1Inverse(1:sampleCount);
    aclr2Adjoint = spectralEnergyScale * fftLength * window .* aclr2Inverse(1:sampleCount);

    objective = real(timeObjective + rho1 * aclr1Energy + rho2 * aclr2Energy);
    effectiveResidual = residual + rho1 * aclr1Adjoint + rho2 * aclr2Adjoint;
    if rho1 == 0 && rho2 == 0
        objective = timeObjective;
        effectiveResidual = residual;
    end
    effectiveResidual = reshape(effectiveResidual, residualSize);

    if nargout >= 3
        info.sampleCount = sampleCount;
        info.sampleRate = sampleRate;
        info.fftLength = fftLength;
        info.frequencyStep = frequencyStep;
        info.windowEnergy = windowEnergy;
        info.timeObjective = timeObjective;
        info.aclr1Energy = aclr1Energy;
        info.aclr2Energy = aclr2Energy;
        info.aclr1WeightedObjective = rho1 * aclr1Energy;
        info.aclr2WeightedObjective = rho2 * aclr2Energy;
        info.totalObjective = objective;
        info.rho1 = rho1;
        info.rho2 = rho2;
        info.frequency = frequency;
        info.bandCenters = bandCenters;
        info.bandEdges = bandEdges;
        info.bandWeights = bandWeights;
        info.aclr1FrequencyWeight = aclr1FrequencyWeight;
        info.aclr2FrequencyWeight = aclr2FrequencyWeight;
    end
end
