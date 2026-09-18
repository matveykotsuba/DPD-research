function [outputSignal, info] = cfrCore(inputSignal, cfrCfg, systemCfg)
    if nargin < 3 || ~isstruct(systemCfg)
        error('cfrCore:MissingSystemConfiguration', 'A full systemCfg is required.');
    end
    if isempty(cfrCfg)
        cfrCfg = cfrConfig(false);
    end
    if ~isnumeric(inputSignal) || ~isvector(inputSignal) || isempty(inputSignal) || any(~isfinite(inputSignal(:)))
        error('cfrCore:InvalidInput', 'Input must be a finite serialized CP/WOLA frame vector.');
    end
    effectiveCfg = cfrConfig(false);
    fields = fieldnames(cfrCfg);
    for k = 1:length(fields)
        effectiveCfg.(fields{k}) = cfrCfg.(fields{k});
    end
    effectiveCfg.fftSize = systemCfg.fftSize;
    effectiveCfg.activeIndices = systemCfg.activeIndices(:);
    effectiveCfg.cpLengths = systemCfg.cpLengths(:);
    effectiveCfg.wolaLength = systemCfg.wolaLength;
    if effectiveCfg.algorithmVersion ~= 2 || ~strcmp(effectiveCfg.processingPoint, 'after-cp-wola') || ...
            ~strcmp(effectiveCfg.boundaryMode, 'periodic-frame') || ...
            ~strcmp(effectiveCfg.pulseShape, 'frame-occupied-band-rectangular-mask')
        error('cfrCore:OldConfiguration', 'Recreate the configuration with cfrConfig for CFR after CP/WOLA.');
    end
    cpLengths = effectiveCfg.cpLengths;
    if isempty(cpLengths) || any(~isfinite(cpLengths)) || any(cpLengths < 0) || ...
            any(cpLengths ~= round(cpLengths)) || any(cpLengths > systemCfg.fftSize)
        error('cfrCore:InvalidCpLengths', 'Invalid integer CP lengths.');
    end
    lengths = systemCfg.fftSize + cpLengths;
    frameLength = sum(lengths);
    if numel(inputSignal) ~= frameLength
        error('cfrCore:FrameLengthMismatch', 'Input length must equal sum(fftSize + cpLengths).');
    end
    signature = cfrSignature(effectiveCfg);
    x = inputSignal(:);
    y = x;
    starts = [1; 1 + cumsum(lengths(1:end-1))];
    symbolCount = length(starts);
    sampleSymbol = zeros(frameLength, 1);
    usefulIndices = zeros(systemCfg.fftSize, symbolCount);
    for s = 1:symbolCount
        sampleSymbol(starts(s):starts(s)+lengths(s)-1) = s;
        usefulIndices(:, s) = starts(s) + cpLengths(s) + (0:systemCfg.fftSize-1).';
    end
    inputPower = mean(abs(x).^2);
    inputRms = sqrt(inputPower);
    [inputPeak, worstPeakIndex] = max(abs(x));
    worstSymbol = sampleSymbol(worstPeakIndex);
    threshold = Inf;
    counts = zeros(symbolCount, 1);
    stalled = false(symbolCount, 1);
    normalizationScale = 1;
    pulseLeakage = 0;
    if signature.enabled && inputPower > 0
        threshold = inputRms * 10^(effectiveCfg.targetPaprDb / 20);
        [pulse, bandMask] = framePulse(frameLength, systemCfg.fftSize, systemCfg.activeIndices);
        pulseSpectrum = fftshift(fft(pulse));
        pulseLeakage = sum(abs(pulseSpectrum(~bandMask)).^2) / sum(abs(pulseSpectrum).^2);
        indexVector = (0:frameLength-1).';
        tolerance = effectiveCfg.relativeTolerance;
        while true
            envelope = abs(y);
            framePeak = max(envelope);
            eligible = ~stalled & counts < effectiveCfg.maximumPulsesPerSymbol;
            envelope(~eligible(sampleSymbol)) = 0;
            [peak, peakIndex] = max(envelope);
            if peak <= threshold * (1 + tolerance)
                break
            end
            s = sampleSymbol(peakIndex);
            amplitude = (peak - threshold) * y(peakIndex) / peak;
            shiftedPulse = pulse(mod(indexVector - (peakIndex - 1), frameLength) + 1);
            accepted = false;
            for factor = effectiveCfg.lineSearchFactors(:).'
                candidate = y - effectiveCfg.cancellationFactor * factor * amplitude .* shiftedPulse;
                if max(abs(candidate)) <= framePeak * (1 + tolerance) && abs(candidate(peakIndex)) < peak
                    y = candidate;
                    accepted = true;
                    counts(s) = counts(s) + 1;
                    break
                end
            end
            if ~accepted
                stalled(s) = true;
            end
        end
        if effectiveCfg.preserveAveragePower && sum(counts) > 0 && mean(abs(y).^2) > 0
            normalizationScale = sqrt(inputPower / mean(abs(y).^2));
            y = y * normalizationScale;
        end
    end
    outputPower = mean(abs(y).^2);
    outputRms = sqrt(outputPower);
    inputSymbolPapr = zeros(symbolCount, 1);
    outputSymbolPapr = zeros(symbolCount, 1);
    inputSymbolPeaks = zeros(symbolCount, 1);
    outputSymbolPeaks = zeros(symbolCount, 1);
    for s = 1:symbolCount
        indices = starts(s):starts(s)+lengths(s)-1;
        inputSymbolPapr(s) = signalPaprDb(x(indices));
        outputSymbolPapr(s) = signalPaprDb(y(indices));
        inputSymbolPeaks(s) = max(abs(x(indices)));
        outputSymbolPeaks(s) = max(abs(y(indices)));
    end
    if inputPower > 0
        gain = (x' * y) / (x' * x);
        evm = 100 * norm(y - gain*x) / max(norm(gain*x), realmin('double'));
    else
        gain = 1;
        evm = 0;
    end
    info.enabled = signature.enabled;
    info.processingPoint = 'after CP and WOLA, before DPD and PA';
    info.signature = signature;
    info.targetPaprDb = effectiveCfg.targetPaprDb;
    info.thresholdMagnitude = threshold;
    info.inputAveragePower = inputPower;
    info.outputAveragePower = outputPower;
    info.inputRms = inputRms;
    info.outputRms = outputRms;
    info.inputPeakMagnitude = inputPeak;
    info.outputPeakMagnitude = max(abs(y));
    info.inputPaprDb = signalPaprDb(x);
    info.outputPaprDb = signalPaprDb(y);
    info.usefulInputPaprDb = signalPaprDb(x(usefulIndices));
    info.usefulOutputPaprDb = signalPaprDb(y(usefulIndices));
    info.paprReductionDb = info.inputPaprDb - info.outputPaprDb;
    info.peakMagnitudeReductionDb = 20*log10(max(inputPeak, realmin('double')) / max(info.outputPeakMagnitude, realmin('double')));
    info.averagePowerChangeDb = 10*log10(max(outputPower, realmin('double')) / max(inputPower, realmin('double')));
    info.normalizationScale = normalizationScale;
    info.inputPaprPerSymbolDb = inputSymbolPapr;
    info.outputPaprPerSymbolDb = outputSymbolPapr;
    info.inputPeakRatioPerSymbolDb = 20*log10(max(inputSymbolPeaks / max(inputRms, realmin('double')), realmin('double')));
    info.outputPeakRatioPerSymbolDb = 20*log10(max(outputSymbolPeaks / max(outputRms, realmin('double')), realmin('double')));
    info.pulsesPerSymbol = counts;
    info.totalPulseCount = sum(counts);
    info.reportingTargetToleranceDb = 0.02;
    converged = info.outputPeakRatioPerSymbolDb <= effectiveCfg.targetPaprDb + info.reportingTargetToleranceDb;
    info.convergedSymbolCount = sum(converged);
    info.stalledSymbolCount = sum(stalled);
    info.remainingPeakSymbolCount = sum(~converged);
    info.correctionPowerDb = 10*log10(max(sum(abs(y-x).^2) / max(sum(abs(x).^2), realmin('double')), realmin('double')));
    info.bestLinearGain = gain;
    info.cfrOnlyEvmPercent = evm;
    info.pulseOutOfBandEnergyRatio = pulseLeakage;
    info.worstSymbolIndex = worstSymbol;
    indices = starts(worstSymbol):starts(worstSymbol)+lengths(worstSymbol)-1;
    info.worstSymbolInputEnvelope = abs(x(indices));
    info.worstSymbolOutputEnvelope = abs(y(indices));
    windowLength = min(512, frameLength);
    first = max(1, min(worstPeakIndex - floor(windowLength/2), frameLength-windowLength+1));
    indices = (first:first+windowLength-1).';
    info.peakWindowSampleIndices = indices - 1;
    info.peakWindowInputEnvelope = abs(x(indices));
    info.peakWindowOutputEnvelope = abs(y(indices));
    outputSignal = reshape(y, size(inputSignal));
end

function [pulse, mask] = framePulse(frameLength, fftSize, activeIndices)
    frequencies = (-floor(frameLength/2):ceil(frameLength/2)-1).' / frameLength;
    nearestBin = mod(floor(frequencies*fftSize + 0.5) + floor(fftSize/2), fftSize) + 1;
    activeMask = false(fftSize, 1);
    activeMask(activeIndices) = true;
    mask = activeMask(nearestBin);
    pulse = ifft(ifftshift(double(mask)));
    pulse = pulse / pulse(1);
end

function value = signalPaprDb(signal)
    power = abs(signal(:)).^2;
    if mean(power) > 0
        value = 10*log10(max(power)/mean(power));
    else
        value = 0;
    end
end
