function [cacheFile, comparisonCache] = dpdSaveArchitectureComparisonCache(dpdCfg, trainingInfoOrMetrics)
    if ~isstruct(dpdCfg) || ~isscalar(dpdCfg)
        error('dpdSaveArchitectureComparisonCache:InvalidConfiguration', 'dpdCfg must be a scalar structure.');
    end
    architecture = configurationArchitecture(dpdCfg);
    requireConfigurationField(dpdCfg, 'paModelSignature');
    if isempty(dpdCfg.paModelSignature)
        error('dpdSaveArchitectureComparisonCache:MissingPaSignature', 'dpdCfg.paModelSignature must be available.');
    end
    structure = configurationStructure(dpdCfg);
    if isfield(dpdCfg, 'cfrSignature') && ~isempty(dpdCfg.cfrSignature)
        currentCfrSignature = dpdCfg.cfrSignature;
    else
        currentCfrSignature = cfrSignature(cfrConfig(false));
    end
    [metrics, available, reason] = comparisonMetrics(trainingInfoOrMetrics);
    if ~available
        error('dpdSaveArchitectureComparisonCache:InvalidComparisonMetrics', 'Identity-start comparison metrics cannot be saved: %s.', reason);
    end
    compactMetrics = validateAndCompactMetrics(metrics, architecture);
    comparisonCache.cacheContractVersion = 2;
    comparisonCache.learningArchitecture = architecture;
    comparisonCache.paModelSignature = dpdCfg.paModelSignature;
    comparisonCache.cfrSignature = currentCfrSignature;
    comparisonCache.orders = structure.orders;
    comparisonCache.signalDelays = structure.signalDelays;
    comparisonCache.envelopeDelays = structure.envelopeDelays;
    comparisonCache.diagonalCount = structure.diagonalCount;
    comparisonCache.metrics = compactMetrics;
    cacheFile = architectureCacheFile(dpdCfg, architecture);
    cacheDirectory = fileparts(cacheFile);
    if exist(cacheDirectory, 'dir') ~= 7
        error('dpdSaveArchitectureComparisonCache:MissingCacheDirectory', 'The architecture comparison cache directory does not exist: %s', cacheDirectory);
    end
    temporaryCacheFile = [tempname(cacheDirectory) '.mat'];
    temporaryCleanup = onCleanup(@() deleteTemporaryFile(temporaryCacheFile));
    save(temporaryCacheFile, 'comparisonCache', '-v7');
    [moveSucceeded, moveMessage] = movefile(temporaryCacheFile, cacheFile, 'f');
    if ~moveSucceeded
        error('dpdSaveArchitectureComparisonCache:CacheSaveFailed', 'Could not replace the architecture comparison cache: %s', moveMessage);
    end
    clear temporaryCleanup
end

function [metrics, available, reason] = comparisonMetrics(trainingInfoOrMetrics)
    [metrics, available, reason] = dpdMakeArchitectureComparisonMetrics(trainingInfoOrMetrics);
    if available
        return
    end
    wrapper.architectureComparisonMetrics = trainingInfoOrMetrics;
    [metrics, available, reason] = dpdMakeArchitectureComparisonMetrics(wrapper);
end

function compactMetrics = validateAndCompactMetrics(metrics, architecture)
    requiredFields = metricFields();
    if ~isstruct(metrics) || ~isscalar(metrics)
        error('dpdSaveArchitectureComparisonCache:InvalidMetricsType', 'Comparison metrics must be a scalar structure.');
    end
    for fieldIndex = 1:length(requiredFields)
        fieldName = requiredFields{fieldIndex};
        if ~isfield(metrics, fieldName)
            error('dpdSaveArchitectureComparisonCache:MissingMetricField', 'Comparison metrics must contain %s.', fieldName);
        end
    end
    if ~ischar(metrics.learningArchitecture) || ~strcmpi(strtrim(metrics.learningArchitecture), architecture)
        error('dpdSaveArchitectureComparisonCache:ArchitectureMismatch', 'Comparison metric architecture must match dpdCfg.learningArchitecture.');
    end
    if ~isequal(metrics.comparisonContractVersion, 1) || ~ischar(metrics.historyInitialization) || ~strcmpi(strtrim(metrics.historyInitialization), 'identity') || ~logicalScalar(metrics.historyStartsFromIdentity)
        error('dpdSaveArchitectureComparisonCache:NonidentityHistory', 'Comparison metrics must use the identity-start contract.');
    end
    if ~numericScalar(metrics.initialCoefficientIdentityError) || metrics.initialCoefficientIdentityError < 0 || metrics.initialCoefficientIdentityError > 1e-10
        error('dpdSaveArchitectureComparisonCache:NonidentityCoefficients', 'Comparison metrics must start from identity coefficients.');
    end
    updateIndex = metrics.updateIndex(:);
    metricCount = length(updateIndex);
    if metricCount < 2 || any(updateIndex < 0) || any(updateIndex ~= round(updateIndex)) || updateIndex(1) ~= 0 || any(diff(updateIndex) <= 0)
        error('dpdSaveArchitectureComparisonCache:InvalidUpdateIndex', 'Comparison update indices must start at zero and increase strictly.');
    end
    numericSeriesFields = {'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db'};
    for fieldIndex = 1:length(numericSeriesFields)
        fieldName = numericSeriesFields{fieldIndex};
        fieldValue = metrics.(fieldName);
        if ~isnumeric(fieldValue) || ~isvector(fieldValue) || length(fieldValue) ~= metricCount || any(~isfinite(fieldValue(:)))
            error('dpdSaveArchitectureComparisonCache:InvalidMetricSeries', 'Comparison metric %s must be a finite vector matching updateIndex.', fieldName);
        end
    end
    if any(metrics.evmRmsPercent(:) < 0)
        error('dpdSaveArchitectureComparisonCache:NegativeEvm', 'Comparison EVM values must be nonnegative.');
    end
    positiveIntegerFields = {'trainingSampleCount', 'trainingBlockLength', 'monitorSeed', 'monitorOfdmSymbolCount'};
    for fieldIndex = 1:length(positiveIntegerFields)
        fieldName = positiveIntegerFields{fieldIndex};
        fieldValue = metrics.(fieldName);
        if ~numericScalar(fieldValue) || fieldValue <= 0 || fieldValue ~= round(fieldValue)
            error('dpdSaveArchitectureComparisonCache:InvalidProtocolField', 'Comparison protocol field %s must be a positive integer.', fieldName);
        end
    end
    if ~numericScalar(metrics.trainingSeed) || metrics.trainingSeed < 0 || metrics.trainingSeed ~= round(metrics.trainingSeed)
        error('dpdSaveArchitectureComparisonCache:InvalidTrainingSeed', 'Comparison trainingSeed must be a nonnegative integer.');
    end
    if ~numericScalar(metrics.monitorSampleRate) || metrics.monitorSampleRate <= 0
        error('dpdSaveArchitectureComparisonCache:InvalidMonitorSampleRate', 'Comparison monitorSampleRate must be positive.');
    end
    textFields = {'initialPointDefinition', 'evmReference', 'aclrReference'};
    for fieldIndex = 1:length(textFields)
        fieldName = textFields{fieldIndex};
        if ~ischar(metrics.(fieldName)) || isempty(strtrim(metrics.(fieldName)))
            error('dpdSaveArchitectureComparisonCache:InvalidMetricText', 'Comparison metric %s must be a nonempty character vector.', fieldName);
        end
    end
    if ~ischar(metrics.spectrumEstimator) || ~strcmpi(strtrim(metrics.spectrumEstimator), 'Welch')
        error('dpdSaveArchitectureComparisonCache:InvalidSpectrumEstimator', 'Comparison spectra must use Welch.');
    end
    compactMetrics = metrics;
    compactMetrics.learningArchitecture = architecture;
    compactMetrics.updateIndex = updateIndex;
    compactMetrics.evmRmsPercent = metrics.evmRmsPercent(:);
    compactMetrics.worstAclr1Db = metrics.worstAclr1Db(:);
    compactMetrics.worstAclr2Db = metrics.worstAclr2Db(:);
end

function fields = metricFields()
    fields = {'learningArchitecture', 'comparisonContractVersion', 'historyInitialization', 'historyStartsFromIdentity', 'initialCoefficientIdentityError', 'initialPointDefinition', 'trainingSeed', 'trainingSampleCount', 'trainingBlockLength', 'updateIndex', 'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db', 'monitorSeed', 'monitorOfdmSymbolCount', 'monitorSampleRate', 'evmReference', 'aclrReference', 'spectrumEstimator'};
end

function structure = configurationStructure(dpdCfg)
    structureFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount'};
    structure = struct();
    for fieldIndex = 1:length(structureFields)
        fieldName = structureFields{fieldIndex};
        requireConfigurationField(dpdCfg, fieldName);
        fieldValue = dpdCfg.(fieldName);
        if ~isnumeric(fieldValue) || isempty(fieldValue) || any(~isfinite(fieldValue(:)))
            error('dpdSaveArchitectureComparisonCache:InvalidDpdStructure', 'dpdCfg.%s must contain finite numeric values.', fieldName);
        end
        structure.(fieldName) = fieldValue;
    end
    if ~isscalar(structure.diagonalCount) || structure.diagonalCount < 0 || structure.diagonalCount ~= round(structure.diagonalCount)
        error('dpdSaveArchitectureComparisonCache:InvalidDiagonalCount', 'dpdCfg.diagonalCount must be a nonnegative integer.');
    end
end

function architecture = configurationArchitecture(dpdCfg)
    requireConfigurationField(dpdCfg, 'learningArchitecture');
    if ~ischar(dpdCfg.learningArchitecture)
        error('dpdSaveArchitectureComparisonCache:InvalidArchitecture', 'dpdCfg.learningArchitecture must be direct or indirect.');
    end
    architecture = lower(strtrim(dpdCfg.learningArchitecture));
    if ~ismember(architecture, {'direct', 'indirect'})
        error('dpdSaveArchitectureComparisonCache:InvalidArchitecture', 'dpdCfg.learningArchitecture must be direct or indirect.');
    end
end

function requireConfigurationField(dpdCfg, fieldName)
    if ~isfield(dpdCfg, fieldName)
        error('dpdSaveArchitectureComparisonCache:MissingConfigurationField', 'dpdCfg must contain %s.', fieldName);
    end
end

function cacheFile = architectureCacheFile(dpdCfg, architecture)
    projectDirectory = fileparts(mfilename('fullpath'));
    if isfield(dpdCfg, 'architectureComparisonCacheDirectory')
        cacheDirectory = dpdCfg.architectureComparisonCacheDirectory;
        if ~ischar(cacheDirectory) || isempty(strtrim(cacheDirectory))
            error('dpdSaveArchitectureComparisonCache:InvalidCacheDirectory', 'architectureComparisonCacheDirectory must be a nonempty character vector.');
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

function deleteTemporaryFile(fileName)
    if exist(fileName, 'file') == 2
        delete(fileName);
    end
end
