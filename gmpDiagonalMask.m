function [coefficientMask, activeIndices, maskInfo] = gmpDiagonalMask(modelCfg, diagonalCount)

if nargin < 2 || isempty(diagonalCount)
    if ~isfield(modelCfg, 'diagonalCount')
        error('gmpDiagonalMask:MissingDiagonalCount', 'diagonalCount is required.');
    end
    diagonalCount = modelCfg.diagonalCount;
end

requiredFields = {'orders', 'signalDelays', 'envelopeDelays'};
for fieldIndex = 1:length(requiredFields)
    if ~isfield(modelCfg, requiredFields{fieldIndex})
        error('gmpDiagonalMask:MissingConfigurationField', 'The model configuration must contain %s.', requiredFields{fieldIndex});
    end
end

orders = modelCfg.orders(:).';
signalDelays = modelCfg.signalDelays(:).';
envelopeDelays = modelCfg.envelopeDelays(:).';

validateattributes(orders, {'numeric'}, {'vector', 'real', 'finite', 'positive', 'integer', 'nonempty'});
validateattributes(signalDelays, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'integer', 'nonempty'});
validateattributes(envelopeDelays, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'integer', 'nonempty'});
validateattributes(diagonalCount, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative', 'integer'});

if orders(1) ~= 1
    error('gmpDiagonalMask:InvalidOrderLayout', 'The first model order must be the linear order 1.');
end

numberOfSignalDelays = length(signalDelays);
numberOfEnvelopeDelays = length(envelopeDelays);
numberOfColumns = 1 + (length(orders) - 1) * numberOfEnvelopeDelays;

coefficientMask = false(numberOfSignalDelays, numberOfColumns);

coefficientMask(:, 1) = true;

for orderIndex = 2:length(orders)
    for envelopeIndex = 1:numberOfEnvelopeDelays
        envelopeDelay = envelopeDelays(envelopeIndex);
        coefficientColumn = 2 + (orderIndex - 2) * numberOfEnvelopeDelays + envelopeIndex - 1;

        for signalIndex = 1:numberOfSignalDelays
            signalDelay = signalDelays(signalIndex);
            coefficientMask(signalIndex, coefficientColumn) = abs(signalDelay - envelopeDelay) <= diagonalCount;
        end
    end
end

activeIndices = find(coefficientMask(:));

maskInfo.diagonalCount = diagonalCount;
maskInfo.numberOfRows = numberOfSignalDelays;
maskInfo.numberOfColumns = numberOfColumns;
maskInfo.numberOfCoefficients = numel(coefficientMask);
maskInfo.numberOfActiveCoefficients = length(activeIndices);
maskInfo.activeFraction = length(activeIndices) / numel(coefficientMask);

end
