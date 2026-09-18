function [basis, firstModelSample] = gmpActiveBasis(inputSignal, modelCfg, coefficientMask)

inputSignal = inputSignal(:);
orders = modelCfg.orders;
signalDelays = modelCfg.signalDelays;
envelopeDelays = modelCfg.envelopeDelays;
maximumDelay = max([signalDelays envelopeDelays]);
firstModelSample = maximumDelay + 1;
numberOfRows = length(inputSignal) - maximumDelay;
numberOfCoefficientColumns = 1 + (length(orders) - 1) * length(envelopeDelays);

if numberOfRows < 1
    error('gmpActiveBasis:InsufficientSamples', 'The input signal is shorter than the configured GMP memory.');
end

if nargin < 3 || isempty(coefficientMask)
    coefficientMask = true(length(signalDelays), numberOfCoefficientColumns);
end

if ~islogical(coefficientMask) || ~isequal(size(coefficientMask), [length(signalDelays), numberOfCoefficientColumns])
    error('gmpActiveBasis:InvalidCoefficientMask', 'coefficientMask has an invalid size or type.');
end

basis = complex(zeros(numberOfRows, nnz(coefficientMask)));
activeColumnIndex = 1;
coefficientColumn = 1;

for signalIndex = 1:length(signalDelays)
    if coefficientMask(signalIndex, coefficientColumn)
        signalDelay = signalDelays(signalIndex);
        basis(:, activeColumnIndex) = inputSignal(firstModelSample - signalDelay:end - signalDelay);
        activeColumnIndex = activeColumnIndex + 1;
    end
end

for orderIndex = 2:length(orders)
    exponent = orders(orderIndex) - 1;
    for envelopeIndex = 1:length(envelopeDelays)
        coefficientColumn = 2 + (orderIndex - 2) * length(envelopeDelays) + envelopeIndex - 1;
        envelopeDelay = envelopeDelays(envelopeIndex);
        envelopeTerm = abs(inputSignal(firstModelSample - envelopeDelay:end - envelopeDelay)).^exponent;
        for signalIndex = 1:length(signalDelays)
            if coefficientMask(signalIndex, coefficientColumn)
                signalDelay = signalDelays(signalIndex);
                signalTerm = inputSignal(firstModelSample - signalDelay:end - signalDelay);
                basis(:, activeColumnIndex) = signalTerm .* envelopeTerm;
                activeColumnIndex = activeColumnIndex + 1;
            end
        end
    end
end

end
