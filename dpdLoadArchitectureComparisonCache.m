function [metrics, structure, available, reason, cacheFile] = dpdLoadArchitectureComparisonCache(dpdCfg, architecture)
    metrics = struct();
    structure = struct();
    available = false;
    if ~isstruct(dpdCfg) || ~isscalar(dpdCfg)
        reason = 'architecture comparison cache cannot be checked because dpdCfg is invalid';
        cacheFile = '';
        return
    end
    [architecture, architectureValid] = validatedArchitecture(architecture);
    if ~architectureValid
        reason = 'architecture comparison cache architecture must be direct or indirect';
        cacheFile = '';
        return
    end
    try
        cacheFile = architectureCacheFile(dpdCfg, architecture);
    catch cachePathError
        cacheFile = '';
        reason = ['architecture comparison cache path is invalid: ' cachePathError.message];
        return
    end
    if exist(cacheFile, 'file') ~= 2
        reason = ['architecture comparison cache file not found: ' cacheFile];
        return
    end
    try
        savedData = load(cacheFile, 'comparisonCache');
    catch loadError
        reason = ['architecture comparison cache cannot be loaded: ' loadError.message];
        return
    end
    if ~isfield(savedData, 'comparisonCache') || ~isstruct(savedData.comparisonCache) || ~isscalar(savedData.comparisonCache)
        reason = 'architecture comparison cache record is missing or invalid';
        return
    end
    comparisonCache = savedData.comparisonCache;
    requiredCacheFields = {'cacheContractVersion', 'learningArchitecture', 'paModelSignature', 'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount', 'metrics'};
    for fieldIndex = 1:length(requiredCacheFields)
        if ~isfield(comparisonCache, requiredCacheFields{fieldIndex})
            reason = ['architecture comparison cache is missing ' requiredCacheFields{fieldIndex}];
            return
        end
    end
    if ~isequal(comparisonCache.cacheContractVersion, 2)
        reason = 'architecture comparison cache contract version is unsupported';
        return
    end
    if ~ischar(comparisonCache.learningArchitecture) || ~strcmpi(strtrim(comparisonCache.learningArchitecture), architecture)
        reason = ['architecture comparison cache does not contain ' architecture ' history'];
        return
    end
    if ~isfield(dpdCfg, 'paModelSignature') || isempty(dpdCfg.paModelSignature)
        reason = 'architecture comparison cache cannot be checked because the current PA signature is unavailable';
        return
    end
    if ~isequaln(comparisonCache.paModelSignature, dpdCfg.paModelSignature)
        reason = 'architecture comparison cache was created for a different PA model or operating point';
        return
    end
    if isfield(comparisonCache, 'cfrSignature')
        savedCfrSignature = comparisonCache.cfrSignature;
    else
        savedCfrSignature = cfrSignature(cfrConfig(false));
    end
    if ~isequaln(savedCfrSignature, currentCfrSignature(dpdCfg))
        reason = 'architecture comparison cache was created with a different CFR configuration';
        return
    end
    structureFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount'};
    for fieldIndex = 1:length(structureFields)
        fieldName = structureFields{fieldIndex};
        if ~isfield(dpdCfg, fieldName)
            structure = struct();
            reason = ['architecture comparison cache cannot be checked because dpdCfg is missing ' fieldName];
            return
        end
        if ~isequaln(comparisonCache.(fieldName), dpdCfg.(fieldName))
            structure = struct();
            reason = ['architecture comparison cache uses a different DPD GMP structure: ' fieldName];
            return
        end
        structure.(fieldName) = comparisonCache.(fieldName);
    end
    wrapper.architectureComparisonMetrics = comparisonCache.metrics;
    [candidateMetrics, contractAvailable, contractReason] = dpdMakeArchitectureComparisonMetrics(wrapper);
    if ~contractAvailable
        structure = struct();
        reason = ['architecture comparison cache metrics violate the identity-start contract: ' contractReason];
        return
    end
    [metricsValid, metricsReason] = validateMetrics(candidateMetrics, architecture);
    if ~metricsValid
        structure = struct();
        reason = ['architecture comparison cache metrics are invalid: ' metricsReason];
        return
    end
    [protocolValid, protocolReason] = reportingProtocolMatches(candidateMetrics, dpdCfg);
    if ~protocolValid
        structure = struct();
        reason = ['architecture comparison cache uses a different full-frame reporting signal: ' protocolReason];
        return
    end
    metrics = candidateMetrics;
    available = true;
    reason = 'compatible architecture comparison cache loaded';
end

function [valid, reason] = reportingProtocolMatches(metrics, dpdCfg)
    valid = false;
    expectedSeed = configurationValue(dpdCfg, 'adaptationMetricSeed', 7);
    expectedSymbolCount = configurationValue(dpdCfg, 'adaptationMetricOfdmSymbols', 140);
    expectedSegmentLength = configurationValue(dpdCfg, 'adaptationMetricSpectrumSegmentLength', 8192);
    if metrics.monitorSeed ~= expectedSeed
        reason = 'seed does not match the current main run';
        return
    end
    if metrics.monitorOfdmSymbolCount ~= expectedSymbolCount
        reason = 'OFDM symbol count does not match the current full frame';
        return
    end
    if ~isfield(metrics, 'spectrumWelchInfo') || ~isstruct(metrics.spectrumWelchInfo) || ~isscalar(metrics.spectrumWelchInfo) || ~isfield(metrics.spectrumWelchInfo, 'segmentLength') || metrics.spectrumWelchInfo.segmentLength ~= expectedSegmentLength
        reason = 'Welch segment length does not match the current ACLR measurement';
        return
    end
    valid = true;
    reason = 'full-frame reporting protocol matches';
end

function [valid, reason] = validateMetrics(metrics, architecture)
    valid = false;
    requiredFields = metricFields();
    if ~isstruct(metrics) || ~isscalar(metrics)
        reason = 'metrics must be a scalar structure';
        return
    end
    for fieldIndex = 1:length(requiredFields)
        fieldName = requiredFields{fieldIndex};
        if ~isfield(metrics, fieldName)
            reason = ['metrics do not contain ' fieldName];
            return
        end
    end
    if ~ischar(metrics.learningArchitecture) || ~strcmpi(strtrim(metrics.learningArchitecture), architecture)
        reason = ['metric architecture does not match ' architecture];
        return
    end
    if ~isequal(metrics.comparisonContractVersion, 1) || ~ischar(metrics.historyInitialization) || ~strcmpi(strtrim(metrics.historyInitialization), 'identity') || ~logicalScalar(metrics.historyStartsFromIdentity)
        reason = 'metrics do not use the identity-start comparison contract';
        return
    end
    if ~numericScalar(metrics.initialCoefficientIdentityError) || metrics.initialCoefficientIdentityError < 0 || metrics.initialCoefficientIdentityError > 1e-10
        reason = 'initial coefficients are not identity';
        return
    end
    updateIndex = metrics.updateIndex(:);
    metricCount = length(updateIndex);
    if ~isnumeric(updateIndex) || metricCount < 2 || any(~isfinite(updateIndex)) || any(updateIndex < 0) || any(updateIndex ~= round(updateIndex)) || updateIndex(1) ~= 0 || any(diff(updateIndex) <= 0)
        reason = 'update indices must start at zero and increase strictly';
        return
    end
    numericSeriesFields = {'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db'};
    for fieldIndex = 1:length(numericSeriesFields)
        fieldName = numericSeriesFields{fieldIndex};
        fieldValue = metrics.(fieldName);
        if ~isnumeric(fieldValue) || ~isvector(fieldValue) || length(fieldValue) ~= metricCount || any(~isfinite(fieldValue(:)))
            reason = [fieldName ' must be a finite vector matching updateIndex'];
            return
        end
    end
    if any(metrics.evmRmsPercent(:) < 0)
        reason = 'EVM values must be nonnegative';
        return
    end
    positiveIntegerFields = {'trainingSampleCount', 'trainingBlockLength', 'monitorSeed', 'monitorOfdmSymbolCount'};
    for fieldIndex = 1:length(positiveIntegerFields)
        fieldName = positiveIntegerFields{fieldIndex};
        fieldValue = metrics.(fieldName);
        if ~numericScalar(fieldValue) || fieldValue <= 0 || fieldValue ~= round(fieldValue)
            reason = [fieldName ' must be a positive integer'];
            return
        end
    end
    if ~numericScalar(metrics.trainingSeed) || metrics.trainingSeed < 0 || metrics.trainingSeed ~= round(metrics.trainingSeed)
        reason = 'trainingSeed must be a nonnegative integer';
        return
    end
    if ~numericScalar(metrics.monitorSampleRate) || metrics.monitorSampleRate <= 0
        reason = 'monitorSampleRate must be positive';
        return
    end
    textFields = {'initialPointDefinition', 'evmReference', 'aclrReference'};
    for fieldIndex = 1:length(textFields)
        fieldName = textFields{fieldIndex};
        if ~ischar(metrics.(fieldName)) || isempty(strtrim(metrics.(fieldName)))
            reason = [fieldName ' must be a nonempty character vector'];
            return
        end
    end
    if ~ischar(metrics.spectrumEstimator) || ~strcmpi(strtrim(metrics.spectrumEstimator), 'Welch')
        reason = 'spectrumEstimator must be Welch';
        return
    end
    valid = true;
    reason = 'valid compact identity-start metrics';
end

function fields = metricFields()
    fields = {'learningArchitecture', 'comparisonContractVersion', 'historyInitialization', 'historyStartsFromIdentity', 'initialCoefficientIdentityError', 'initialPointDefinition', 'trainingSeed', 'trainingSampleCount', 'trainingBlockLength', 'updateIndex', 'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db', 'monitorSeed', 'monitorOfdmSymbolCount', 'monitorSampleRate', 'evmReference', 'aclrReference', 'spectrumEstimator'};
end

function [architecture, valid] = validatedArchitecture(candidate)
    architecture = '';
    valid = false;
    if ~ischar(candidate)
        return
    end
    architecture = lower(strtrim(candidate));
    valid = ismember(architecture, {'direct', 'indirect'});
end

function cacheFile = architectureCacheFile(dpdCfg, architecture)
    projectDirectory = fileparts(mfilename('fullpath'));
    if isfield(dpdCfg, 'architectureComparisonCacheDirectory')
        cacheDirectory = dpdCfg.architectureComparisonCacheDirectory;
        if ~ischar(cacheDirectory) || isempty(strtrim(cacheDirectory))
            error('dpdLoadArchitectureComparisonCache:InvalidCacheDirectory', 'architectureComparisonCacheDirectory must be a nonempty character vector.');
        end
    else
        cacheDirectory = projectDirectory;
    end
    [~, cacheName, cacheExtension] = fileparts(dpdArchitectureArtifactFile(architecture, currentCfrSignature(dpdCfg), 'comparison'));
    cacheFile = fullfile(cacheDirectory, [cacheName cacheExtension]);
end

function signature = currentCfrSignature(dpdCfg)
    if isfield(dpdCfg, 'cfrSignature') && ~isempty(dpdCfg.cfrSignature)
        signature = dpdCfg.cfrSignature;
    else
        signature = cfrSignature(cfrConfig(false));
    end
end

function valid = numericScalar(value)
    valid = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end

function value = logicalScalar(candidate)
    value = (islogical(candidate) || isnumeric(candidate)) && isscalar(candidate) && isreal(candidate) && isfinite(candidate) && candidate ~= 0;
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
