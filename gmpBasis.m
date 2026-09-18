function [basis, firstModelSample] = gmpBasis(inputSignal, paCfg, coefficientMask)

inputSignal = inputSignal(:);

orders = paCfg.orders;
signalDelays = paCfg.signalDelays;
envelopeDelays = paCfg.envelopeDelays;

maximumDelay = max([signalDelays envelopeDelays]);
firstModelSample = maximumDelay + 1;

numberOfRows = length(inputSignal) - maximumDelay;
numberOfCoefficientColumns = 1 + (length(orders) - 1) * length(envelopeDelays);
if nargin < 3 || isempty(coefficientMask)
    coefficientMask = true(length(signalDelays), numberOfCoefficientColumns);
end
if ~islogical(coefficientMask) || ~isequal(size(coefficientMask), [length(signalDelays), numberOfCoefficientColumns])
    error('gmpBasis:InvalidCoefficientMask', 'coefficientMask has an invalid size or type.');
end
numberOfCoefficients = nnz(coefficientMask);

basis = complex(zeros(numberOfRows, numberOfCoefficients));

columnIndex = 1;
coefficientColumn = 1;

for signalIndex = 1:length(signalDelays)
    if coefficientMask(signalIndex, coefficientColumn)
        signalDelay = signalDelays(signalIndex);
        basis(:, columnIndex) = inputSignal( firstModelSample-signalDelay:end-signalDelay);
        columnIndex = columnIndex + 1;
    end
end

for orderIndex = 2:length(orders)
    exponent = orders(orderIndex) - 1;

    for envelopeIndex = 1:length(envelopeDelays)
        coefficientColumn = 2 + (orderIndex - 2) * length(envelopeDelays) + envelopeIndex - 1;
        envelopeDelay = envelopeDelays(envelopeIndex);

        envelopeTerm = abs(inputSignal( firstModelSample-envelopeDelay:end-envelopeDelay)).^exponent;

        for signalIndex = 1:length(signalDelays)
            if coefficientMask(signalIndex, coefficientColumn)
                signalDelay = signalDelays(signalIndex);
                signalTerm = inputSignal( firstModelSample-signalDelay:end-signalDelay);
                basis(:, columnIndex) = signalTerm .* envelopeTerm;
                columnIndex = columnIndex + 1;
            end
        end
    end
end

end
