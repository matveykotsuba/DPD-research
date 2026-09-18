function [metrics, available, reason] = dpdMakeArchitectureComparisonMetrics(trainingInfo)
    metrics = struct();
    available = false;
    reason = 'identity-start comparison history is unavailable';
    if ~isstruct(trainingInfo) || ~isscalar(trainingInfo)
        return
    end
    if isfield(trainingInfo, 'architectureComparisonMetrics') && isstruct(trainingInfo.architectureComparisonMetrics) && ~isempty(fieldnames(trainingInfo.architectureComparisonMetrics))
        candidateMetrics = trainingInfo.architectureComparisonMetrics;
        if validComparisonContract(candidateMetrics)
            metrics = candidateMetrics;
            available = true;
            reason = 'stored identity-start comparison history';
            return
        else
            reason = 'stored comparison history does not satisfy the identity-start contract';
        end
    end
    if ~isfield(trainingInfo, 'learning') || ~isstruct(trainingInfo.learning) || ~isfield(trainingInfo, 'blockMetrics') || ~isstruct(trainingInfo.blockMetrics) || isempty(fieldnames(trainingInfo.blockMetrics))
        return
    end
    learningInfo = trainingInfo.learning;
    if ~isfield(learningInfo, 'coefficientHistory') || isempty(learningInfo.coefficientHistory)
        reason = 'learning coefficient history is unavailable';
        return
    end
    coefficientHistory = learningInfo.coefficientHistory;
    if ~isnumeric(coefficientHistory) || any(~isfinite(coefficientHistory(:)))
        reason = 'learning coefficient history is invalid';
        return
    end
    initialCoefficients = coefficientHistory(:, 1);
    identityCoefficients = complex(zeros(size(initialCoefficients)));
    identityCoefficients(1) = 1;
    identityError = norm(initialCoefficients - identityCoefficients);
    identityTolerance = 1e-10;
    declaredIdentity = false;
    if isfield(learningInfo, 'identityInitialization')
        declaredIdentity = logicalScalar(learningInfo.identityInitialization);
    elseif isfield(learningInfo, 'initialization') && ischar(learningInfo.initialization)
        declaredIdentity = strcmpi(strtrim(learningInfo.initialization), 'identity');
    end
    actualIdentity = isfinite(identityError) && identityError <= identityTolerance;
    if ~actualIdentity
        reason = 'learning history starts from configured coefficients instead of identity';
        return
    end
    candidateMetrics = trainingInfo.blockMetrics;
    if ~isfield(candidateMetrics, 'updateIndex') || isempty(candidateMetrics.updateIndex) || candidateMetrics.updateIndex(1) ~= 0
        reason = 'block history does not start at update zero';
        return
    end
    candidateMetrics.comparisonContractVersion = 1;
    candidateMetrics.historyInitialization = 'identity';
    candidateMetrics.historyStartsFromIdentity = true;
    candidateMetrics.identityInitializationDeclared = declaredIdentity;
    candidateMetrics.initialCoefficientIdentityError = identityError;
    candidateMetrics.initialPointDefinition = 'identity DPD, equivalent to PA only on the common monitoring signal';
    candidateMetrics.trainingSeed = optionValue(trainingInfo, 'trainingSeed', NaN);
    candidateMetrics.trainingSampleCount = optionValue(trainingInfo, 'trainingSampleCount', NaN);
    candidateMetrics.trainingBlockLength = optionValue(trainingInfo, 'blockLength', NaN);
    candidateMetrics.spectrumEstimator = 'Welch';
    if ~validComparisonContract(candidateMetrics)
        reason = 'legacy identity history is missing comparison protocol fields';
        return
    end
    metrics = candidateMetrics;
    available = true;
    reason = 'legacy identity-start history promoted to comparison contract';
end

function valid = validComparisonContract(metrics)
    requiredFields = {'learningArchitecture', 'comparisonContractVersion', 'historyInitialization', 'historyStartsFromIdentity', 'initialCoefficientIdentityError', 'initialPointDefinition', 'trainingSeed', 'trainingSampleCount', 'trainingBlockLength', 'updateIndex', 'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db', 'monitorSeed', 'monitorOfdmSymbolCount', 'monitorSampleRate', 'evmReference', 'aclrReference', 'spectrumEstimator'};
    valid = isstruct(metrics) && isscalar(metrics);
    if ~valid
        return
    end
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(metrics, requiredFields{fieldIndex})
            valid = false;
            return
        end
    end
    architectureValid = ischar(metrics.learningArchitecture) && (strcmpi(strtrim(metrics.learningArchitecture), 'direct') || strcmpi(strtrim(metrics.learningArchitecture), 'indirect'));
    valid = isequal(metrics.comparisonContractVersion, 1) && architectureValid && ischar(metrics.historyInitialization) && strcmpi(strtrim(metrics.historyInitialization), 'identity') && logicalScalar(metrics.historyStartsFromIdentity) && numericScalar(metrics.initialCoefficientIdentityError) && metrics.initialCoefficientIdentityError >= 0 && metrics.initialCoefficientIdentityError <= 1e-10;
    if ~valid
        return
    end
    updateIndex = metrics.updateIndex(:);
    metricCount = length(updateIndex);
    if ~isnumeric(updateIndex) || metricCount < 2 || any(~isfinite(updateIndex)) || any(updateIndex < 0) || any(updateIndex ~= round(updateIndex)) || updateIndex(1) ~= 0 || any(diff(updateIndex) <= 0)
        valid = false;
        return
    end
    numericSeriesFields = {'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db'};
    for fieldIndex = 1:length(numericSeriesFields)
        fieldValue = metrics.(numericSeriesFields{fieldIndex});
        if ~isnumeric(fieldValue) || ~isvector(fieldValue) || length(fieldValue) ~= metricCount || any(~isfinite(fieldValue(:)))
            valid = false;
            return
        end
    end
    if any(metrics.evmRmsPercent(:) < 0)
        valid = false;
        return
    end
    positiveIntegerFields = {'trainingSampleCount', 'trainingBlockLength', 'monitorSeed', 'monitorOfdmSymbolCount'};
    for fieldIndex = 1:length(positiveIntegerFields)
        fieldValue = metrics.(positiveIntegerFields{fieldIndex});
        if ~numericScalar(fieldValue) || fieldValue <= 0 || fieldValue ~= round(fieldValue)
            valid = false;
            return
        end
    end
    if ~numericScalar(metrics.trainingSeed) || metrics.trainingSeed < 0 || metrics.trainingSeed ~= round(metrics.trainingSeed) || ~numericScalar(metrics.monitorSampleRate) || metrics.monitorSampleRate <= 0
        valid = false;
        return
    end
    textFields = {'initialPointDefinition', 'evmReference', 'aclrReference'};
    for fieldIndex = 1:length(textFields)
        fieldValue = metrics.(textFields{fieldIndex});
        if ~ischar(fieldValue) || isempty(strtrim(fieldValue))
            valid = false;
            return
        end
    end
    valid = ischar(metrics.spectrumEstimator) && strcmpi(strtrim(metrics.spectrumEstimator), 'Welch');
end

function value = optionValue(trainingInfo, fieldName, defaultValue)
    value = defaultValue;
    if isfield(trainingInfo, 'options') && isstruct(trainingInfo.options) && isfield(trainingInfo.options, fieldName)
        candidateValue = trainingInfo.options.(fieldName);
        if isnumeric(candidateValue) && isscalar(candidateValue) && isreal(candidateValue) && isfinite(candidateValue)
            value = double(candidateValue);
        end
    end
end

function value = logicalScalar(candidate)
    value = (islogical(candidate) || isnumeric(candidate)) && isscalar(candidate) && isreal(candidate) && isfinite(candidate) && candidate ~= 0;
end

function valid = numericScalar(value)
    valid = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end
