function outputSignal = gmpCore(inputSignal, paCfg)

inputSignal = inputSignal(:);

orders = paCfg.orders;
signalDelays = paCfg.signalDelays;
envelopeDelays = paCfg.envelopeDelays;
coefficients = paCfg.coefficients;

outputSignal = complex(zeros(size(inputSignal)));

for signalIndex = 1:length(signalDelays)
    signalDelay = signalDelays(signalIndex);
    outputSignal(signalDelay + 1:end) = outputSignal(signalDelay + 1:end) + coefficients(signalIndex, 1) .* inputSignal(1:end - signalDelay);
end

useHorner = ~isempty(orders) && isequal(orders(:).', 1:2:orders(end));

if useHorner
    outputSignal = evaluateHornerTerms(inputSignal, outputSignal, orders, signalDelays, envelopeDelays, coefficients, paCfg);
else
    outputSignal = evaluateDirectTerms(inputSignal, outputSignal, orders, signalDelays, envelopeDelays, coefficients, paCfg);
end

end

function outputSignal = evaluateHornerTerms(inputSignal, outputSignal, orders, signalDelays, envelopeDelays, coefficients, paCfg)

numberOfNonlinearOrders = length(orders) - 1;
numberOfEnvelopeDelays = length(envelopeDelays);
hasDiagonalCount = isfield(paCfg, 'diagonalCount');

for envelopeIndex = 1:numberOfEnvelopeDelays
    envelopeDelay = envelopeDelays(envelopeIndex);
    coefficientColumns = 1 + envelopeIndex:numberOfEnvelopeDelays:1 + (numberOfNonlinearOrders - 1) * numberOfEnvelopeDelays + envelopeIndex;

    for signalIndex = 1:length(signalDelays)
        signalDelay = signalDelays(signalIndex);

        if hasDiagonalCount && abs(signalDelay - envelopeDelay) > paCfg.diagonalCount
            continue
        end

        nonlinearCoefficients = coefficients(signalIndex, coefficientColumns);
        if isempty(nonlinearCoefficients) || ~any(nonlinearCoefficients ~= 0)
            continue
        end

        firstSample = max(signalDelay, envelopeDelay) + 1;
        signalPart = inputSignal(firstSample - signalDelay:end - signalDelay);
        envelopeMagnitudeSquared = abs(inputSignal(firstSample - envelopeDelay:end - envelopeDelay)).^2;
        polynomialValue = nonlinearCoefficients(end);

        for nonlinearIndex = numberOfNonlinearOrders - 1:-1:1
            polynomialValue = polynomialValue .* envelopeMagnitudeSquared + nonlinearCoefficients(nonlinearIndex);
        end

        outputSignal(firstSample:end) = outputSignal(firstSample:end) + signalPart .* envelopeMagnitudeSquared .* polynomialValue;
    end
end

end

function outputSignal = evaluateDirectTerms(inputSignal, outputSignal, orders, signalDelays, envelopeDelays, coefficients, paCfg)

hasDiagonalCount = isfield(paCfg, 'diagonalCount');

for orderIndex = 2:length(orders)
    exponent = orders(orderIndex) - 1;

    for envelopeIndex = 1:length(envelopeDelays)
        envelopeDelay = envelopeDelays(envelopeIndex);
        coefficientColumn = 2 + (orderIndex - 2) * length(envelopeDelays) + envelopeIndex - 1;

        for signalIndex = 1:length(signalDelays)
            signalDelay = signalDelays(signalIndex);

            if hasDiagonalCount && abs(signalDelay - envelopeDelay) > paCfg.diagonalCount
                continue
            end

            firstSample = max(signalDelay, envelopeDelay) + 1;
            signalPart = inputSignal(firstSample - signalDelay:end - signalDelay);
            envelopePart = abs(inputSignal(firstSample - envelopeDelay:end - envelopeDelay)).^exponent;
            outputSignal(firstSample:end) = outputSignal(firstSample:end) + coefficients(signalIndex, coefficientColumn) .* signalPart .* envelopePart;
        end
    end
end

end
