function signature = dpdPaModelSignature(paCfg)
    requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount', 'coefficients', 'modelSampleRate', 'referenceInputRms', 'inputBackoffDb'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(paCfg, requiredFields{fieldIndex})
            error('dpdPaModelSignature:MissingPaField', 'paCfg must contain %s.', requiredFields{fieldIndex});
        end
    end

    signature.orders = paCfg.orders;
    signature.signalDelays = paCfg.signalDelays;
    signature.envelopeDelays = paCfg.envelopeDelays;
    signature.diagonalCount = paCfg.diagonalCount;
    signature.coefficients = paCfg.coefficients;
    signature.modelSampleRate = paCfg.modelSampleRate;
    signature.referenceInputRms = paCfg.referenceInputRms;
    signature.inputBackoffDb = paCfg.inputBackoffDb;
end
