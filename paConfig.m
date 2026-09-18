function paCfg = paConfig(loadModel)
    if nargin < 1
        loadModel = true;
    end

    paCfg.enabled = true;

    paCfg.orders = [1 3 5 7 9 11 13 15 17];
    paCfg.signalDelays = 0:7;
    paCfg.envelopeDelays = 0:7;
    paCfg.diagonalCount = 6;

    paCfg.modelSampleRate = 737.28e6;
    paCfg.inputBackoffDb = 3;
    paCfg.referenceInputRms = 1.0;
    paCfg.maximumInputMagnitude = inf;

    paCfg.resamplerNeighbors = 24;
    paCfg.resamplerAttenuationDb = 80;

    paCfg.sweepInputLevelDb = (-26:2:6).';

    numberOfColumns = 1 +  (length(paCfg.orders) - 1) * length(paCfg.envelopeDelays);

    paCfg.coefficients = complex(zeros( length(paCfg.signalDelays), numberOfColumns));

    if loadModel
        modelFile = fullfile(  fileparts(mfilename('fullpath')),    'paModel.mat');

        model = load(  modelFile, 'coefficientsTikhonov',  'fitInfo');

        paCfg.orders = model.fitInfo.orders;
        paCfg.signalDelays = model.fitInfo.signalDelays;
        paCfg.envelopeDelays = model.fitInfo.envelopeDelays;
        if isfield(model.fitInfo, 'diagonalCount')
            paCfg.diagonalCount = model.fitInfo.diagonalCount;
        else
            paCfg.diagonalCount = max([paCfg.signalDelays paCfg.envelopeDelays]);
        end
        paCfg.modelSampleRate = model.fitInfo.modelSampleRate;
        paCfg.referenceInputRms =  model.fitInfo.reference_input_rms;
        if isfield(model.fitInfo, 'maximumInputMagnitude')
            paCfg.maximumInputMagnitude = model.fitInfo.maximumInputMagnitude;
        end
        paCfg.coefficients =   model.coefficientsTikhonov;
        paCfg.modelFile = modelFile;
    end
end
