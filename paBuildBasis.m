function [basis, sampleIndices] = paBuildBasis(inputSignal, modelType, orders, memoryOrder, diagonalCount, alignmentMemoryOrder)
    inputSignal = inputSignal(:);
    orders = orders(:).';
    alignmentMemoryOrder = max(alignmentMemoryOrder, memoryOrder);
    sampleIndices = (alignmentMemoryOrder + 1:length(inputSignal)).';
    rowCount = length(sampleIndices);

    if strcmp(modelType, 'memoryless')
        basis = complex(zeros(rowCount, length(orders)));

        for orderIndex = 1:length(orders)
            basis(:, orderIndex) = inputSignal(sampleIndices) .* abs(inputSignal(sampleIndices)).^(orders(orderIndex) - 1);
        end

        return
    end

    if strcmp(modelType, 'memory')
        delayCount = memoryOrder + 1;
        basis = complex(zeros(rowCount, length(orders) * delayCount));
        columnIndex = 1;

        for orderIndex = 1:length(orders)
            exponent = orders(orderIndex) - 1;

            for delay = 0:memoryOrder
                delayedSignal = inputSignal(sampleIndices - delay);
                basis(:, columnIndex) = delayedSignal .* abs(delayedSignal).^exponent;
                columnIndex = columnIndex + 1;
            end
        end

        return
    end

    delayCount = memoryOrder + 1;
    pairCount = 0;

    for envelopeDelay = 0:memoryOrder
        for signalDelay = 0:memoryOrder
            if abs(signalDelay - envelopeDelay) <= diagonalCount
                pairCount = pairCount + 1;
            end
        end
    end

    coefficientCount = delayCount + (length(orders) - 1) * pairCount;
    basis = complex(zeros(rowCount, coefficientCount));
    columnIndex = 1;

    for signalDelay = 0:memoryOrder
        basis(:, columnIndex) = inputSignal(sampleIndices - signalDelay);
        columnIndex = columnIndex + 1;
    end

    for orderIndex = 2:length(orders)
        exponent = orders(orderIndex) - 1;

        for envelopeDelay = 0:memoryOrder
            envelopeTerm = abs(inputSignal(sampleIndices - envelopeDelay)).^exponent;

            for signalDelay = 0:memoryOrder
                if abs(signalDelay - envelopeDelay) <= diagonalCount
                    basis(:, columnIndex) = inputSignal(sampleIndices - signalDelay) .* envelopeTerm;
                    columnIndex = columnIndex + 1;
                end
            end
        end
    end
end
