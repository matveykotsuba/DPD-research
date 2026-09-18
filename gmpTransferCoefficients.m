function [targetCoefficients, info] = gmpTransferCoefficients(sourceCfg, targetCfg)
    validateModel(sourceCfg, 'sourceCfg');
    validateModel(targetCfg, 'targetCfg');
    targetColumnCount = 1 + (length(targetCfg.orders) - 1) * length(targetCfg.envelopeDelays);
    targetCoefficients = complex(zeros(length(targetCfg.signalDelays), targetColumnCount));
    copiedCoefficientCount = 0;
    for targetSignalIndex = 1:length(targetCfg.signalDelays)
        sourceSignalIndex = find(sourceCfg.signalDelays == targetCfg.signalDelays(targetSignalIndex), 1);
        if isempty(sourceSignalIndex)
            continue
        end
        targetCoefficients(targetSignalIndex, 1) = sourceCfg.coefficients(sourceSignalIndex, 1);
        copiedCoefficientCount = copiedCoefficientCount + 1;
        for targetOrderIndex = 2:length(targetCfg.orders)
            sourceOrderIndex = find(sourceCfg.orders == targetCfg.orders(targetOrderIndex), 1);
            if isempty(sourceOrderIndex)
                continue
            end
            for targetEnvelopeIndex = 1:length(targetCfg.envelopeDelays)
                sourceEnvelopeIndex = find(sourceCfg.envelopeDelays == targetCfg.envelopeDelays(targetEnvelopeIndex), 1);
                if isempty(sourceEnvelopeIndex)
                    continue
                end
                sourceColumn = 2 + (sourceOrderIndex - 2) * length(sourceCfg.envelopeDelays) + sourceEnvelopeIndex - 1;
                targetColumn = 2 + (targetOrderIndex - 2) * length(targetCfg.envelopeDelays) + targetEnvelopeIndex - 1;
                targetCoefficients(targetSignalIndex, targetColumn) = sourceCfg.coefficients(sourceSignalIndex, sourceColumn);
                copiedCoefficientCount = copiedCoefficientCount + 1;
            end
        end
    end
    if isfield(targetCfg, 'diagonalCount')
        targetMask = gmpDiagonalMask(targetCfg, targetCfg.diagonalCount);
        targetCoefficients(~targetMask) = 0;
    else
        targetMask = true(size(targetCoefficients));
    end
    info.sourceCoefficientCount = numel(sourceCfg.coefficients);
    info.targetCoefficientCount = numel(targetCoefficients);
    info.copiedCoefficientCount = copiedCoefficientCount;
    info.activeCopiedCoefficientCount = nnz(targetCoefficients(targetMask));
end

function validateModel(modelCfg, fieldName)
    requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'coefficients'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(modelCfg, requiredFields{fieldIndex})
            error('gmpTransferCoefficients:MissingField', '%s must contain %s.', fieldName, requiredFields{fieldIndex});
        end
    end
    expectedColumns = 1 + (length(modelCfg.orders) - 1) * length(modelCfg.envelopeDelays);
    expectedSize = [length(modelCfg.signalDelays), expectedColumns];
    if ~isequal(size(modelCfg.coefficients), expectedSize) || any(~isfinite(modelCfg.coefficients(:)))
        error('gmpTransferCoefficients:InvalidCoefficients', '%s coefficients have an invalid size or contain nonfinite values.', fieldName);
    end
end
