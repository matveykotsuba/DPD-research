function [systemOutput, info] = paSystemOutput( modelOutput, systemLength, inputScale, rateFilter, cfg)
    modelOutput = modelOutput(:);

    factor = round(cfg.pa.rateFactor);
    if canUsePolyphaseDecimator(modelOutput, factor, rateFilter)
        decimatedOutput = polyphaseCircularDecimate(modelOutput, factor, rateFilter);
    else
        decimatedOutput = paDecimate(modelOutput, factor, rateFilter);
    end

    outputLength = min(systemLength, length(decimatedOutput));
    decimatedOutput = decimatedOutput(1:outputLength);

    if inputScale == 0
        inputScale = 1;
    end

    systemOutput = decimatedOutput / inputScale;

    info.outputLength = outputLength;
    info.inputScale = inputScale;
    info.modelOutputAveragePower = mean(abs(modelOutput).^2);
    info.systemOutputAveragePower = mean(abs(systemOutput).^2);
end

function canUse = canUsePolyphaseDecimator(inputSignal, factor, coefficients)
    coefficientCount = length(coefficients);
    canUse = factor >= 1 && isfinite(factor) && factor == round(factor) && mod(length(inputSignal), factor) == 0 && mod(coefficientCount, 2) == 1;
    if canUse
        outputSampleCount = length(inputSignal) / factor;
        maximumPhaseMemory = ceil((coefficientCount - 1) / factor);
        canUse = outputSampleCount > maximumPhaseMemory;
    end
end

function outputSignal = polyphaseCircularDecimate(inputSignal, factor, coefficients)
    coefficients = coefficients(:);
    inputSampleCount = length(inputSignal);
    outputSampleCount = inputSampleCount / factor;
    delaySamples = (length(coefficients) - 1) / 2;
    outputSignal = complex(zeros(outputSampleCount, 1));
    tapIndices = (0:length(coefficients) - 1).';
    tapPhases = mod(delaySamples - tapIndices, factor);
    tapShifts = (delaySamples - tapIndices - tapPhases) / factor;
    for phaseIndex = 0:factor - 1
        phaseTapRows = find(tapPhases == phaseIndex);
        phaseSignal = inputSignal(phaseIndex + 1:factor:end);
        phaseOutput = complex(zeros(outputSampleCount, 1));
        for phaseTapIndex = 1:length(phaseTapRows)
            tapRow = phaseTapRows(phaseTapIndex);
            lag = -tapShifts(tapRow);
            phaseOutput = accumulateCircularLag(phaseOutput, phaseSignal, coefficients(tapRow), lag);
        end
        outputSignal = outputSignal + phaseOutput;
    end
end

function accumulator = accumulateCircularLag(accumulator, signal, coefficient, lag)
    sampleCount = length(signal);
    if lag == 0
        accumulator = accumulator + coefficient * signal;
    elseif lag > 0
        accumulator(lag + 1:end) = accumulator(lag + 1:end) + coefficient * signal(1:end - lag);
        accumulator(1:lag) = accumulator(1:lag) + coefficient * signal(end - lag + 1:end);
    else
        advance = -lag;
        accumulator(1:sampleCount - advance) = accumulator(1:sampleCount - advance) + coefficient * signal(advance + 1:end);
        accumulator(sampleCount - advance + 1:end) = accumulator(sampleCount - advance + 1:end) + coefficient * signal(1:advance);
    end
end
