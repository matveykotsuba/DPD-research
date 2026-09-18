function [outputBlock, state, blockInfo] = gmpCoreBlock(inputBlock, modelCfg, state)

if nargin < 3
    state = [];
end

requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'coefficients'};
for fieldIndex = 1:length(requiredFields)
    if ~isfield(modelCfg, requiredFields{fieldIndex})
        error('gmpCoreBlock:MissingConfigurationField', 'The model configuration must contain %s.', requiredFields{fieldIndex});
    end
end

if ~isnumeric(inputBlock) || ~isvector(inputBlock)
    error('gmpCoreBlock:InvalidInput', 'inputBlock must be a numeric vector.');
end

inputBlock = inputBlock(:);
maximumDelay = max([modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);

if isempty(state)
    state.inputHistory = complex(zeros(0, 1));
    state.historyLength = maximumDelay;
    state.totalSamples = 0;
    state.orders = modelCfg.orders(:).';
    state.signalDelays = modelCfg.signalDelays(:).';
    state.envelopeDelays = modelCfg.envelopeDelays(:).';
    if isfield(modelCfg, 'diagonalCount')
        state.diagonalCount = modelCfg.diagonalCount;
    else
        state.diagonalCount = inf;
    end
else
    requiredStateFields = {'inputHistory', 'historyLength', 'totalSamples', 'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount'};
    for fieldIndex = 1:length(requiredStateFields)
        if ~isfield(state, requiredStateFields{fieldIndex})
            error('gmpCoreBlock:InvalidState', 'The state must contain %s.', requiredStateFields{fieldIndex});
        end
    end

    if isfield(modelCfg, 'diagonalCount')
        diagonalCount = modelCfg.diagonalCount;
    else
        diagonalCount = inf;
    end

    configurationChanged = state.historyLength ~= maximumDelay || ~isequal(state.orders, modelCfg.orders(:).') || ~isequal(state.signalDelays, modelCfg.signalDelays(:).') || ~isequal(state.envelopeDelays, modelCfg.envelopeDelays(:).') || ~isequal(state.diagonalCount, diagonalCount);

    if configurationChanged
        error('gmpCoreBlock:ConfigurationChanged', ['The GMP delay/order configuration changed while a block ' 'state was active. Start with an empty state.']);
    end

    state.inputHistory = state.inputHistory(:);
    if length(state.inputHistory) > maximumDelay
        state.inputHistory = state.inputHistory(end-maximumDelay+1:end);
    end
end

historySampleCount = length(state.inputHistory);

if isempty(inputBlock)
    outputBlock = complex(zeros(0, 1));
else
    inputWithHistory = [state.inputHistory; inputBlock];
    outputWithHistory = gmpCore(inputWithHistory, modelCfg);
    outputBlock = outputWithHistory(historySampleCount + 1:end);

    if maximumDelay > 0
        retainedSampleCount = min(maximumDelay, length(inputWithHistory));
        state.inputHistory = inputWithHistory( end-retainedSampleCount+1:end);
    else
        state.inputHistory = complex(zeros(0, 1));
    end
end

state.totalSamples = state.totalSamples + length(inputBlock);

blockInfo.blockLength = length(inputBlock);
blockInfo.historySamplesUsed = historySampleCount;
blockInfo.historySamplesRetained = length(state.inputHistory);
blockInfo.totalSamples = state.totalSamples;

end
