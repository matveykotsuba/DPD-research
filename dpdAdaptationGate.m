function [accepted, reason, comparison] = dpdAdaptationGate(currentMetrics, trialMetrics, gateCfg)

    if nargin < 2
        error('dpdAdaptationGate:MissingInput', 'currentMetrics and trialMetrics are required.');
    end
    if nargin < 3 || isempty(gateCfg)
        gateCfg = struct();
    end
    if ~isstruct(currentMetrics) || ~isscalar(currentMetrics)
        error('dpdAdaptationGate:InvalidCurrentMetrics', 'currentMetrics must be a scalar structure.');
    end
    if ~isstruct(trialMetrics) || ~isscalar(trialMetrics)
        error('dpdAdaptationGate:InvalidTrialMetrics', 'trialMetrics must be a scalar structure.');
    end
    if ~isstruct(gateCfg) || ~isscalar(gateCfg)
        error('dpdAdaptationGate:InvalidConfiguration', 'gateCfg must be a scalar structure.');
    end

    current = readMetrics(currentMetrics, 'currentMetrics');
    trial = readMetrics(trialMetrics, 'trialMetrics');
    settings = readSettings(gateCfg, currentMetrics, trialMetrics);

    aclr1ImprovementDb = trial.worstAclr1Db - current.worstAclr1Db;
    aclr2ImprovementDb = trial.worstAclr2Db - current.worstAclr2Db;
    evmReductionPercent = current.evmRmsPercent - trial.evmRmsPercent;

    if isfinite(settings.prePaAclr1Db)
        currentTargetGapDb = max(settings.prePaAclr1Db - current.worstAclr1Db, 0);
        trialTargetGapDb = max(settings.prePaAclr1Db - trial.worstAclr1Db, 0);
        targetGapReductionDb = currentTargetGapDb - trialTargetGapDb;
        currentPrimaryScore = -currentTargetGapDb;
        trialPrimaryScore = -trialTargetGapDb;
    else
        currentTargetGapDb = NaN;
        trialTargetGapDb = NaN;
        targetGapReductionDb = NaN;
        currentPrimaryScore = current.worstAclr1Db;
        trialPrimaryScore = trial.worstAclr1Db;
    end

    currentScore = [currentPrimaryScore, -current.evmRmsPercent, current.worstAclr2Db, current.worstAclr1Db];
    trialScore = [trialPrimaryScore, -trial.evmRmsPercent, trial.worstAclr2Db, trial.worstAclr1Db];
    lexicographicallyBetter = lexicographicGreater(trialScore, currentScore, settings.comparisonTolerance);

    finiteAccepted = current.allFinite && trial.allFinite;
    safetyAccepted = current.isSafe && trial.isSafe;
    evmAccepted = evmReductionPercent + settings.comparisonTolerance >= -settings.maximumEvmIncreasePercent;
    aclr1Accepted = aclr1ImprovementDb + settings.comparisonTolerance >= -settings.maximumAclr1DecreaseDb;
    aclr2Accepted = aclr2ImprovementDb + settings.comparisonTolerance >= -settings.maximumAclr2DecreaseDb;
    if isfinite(settings.prePaAclr1Db)
        aclr1ProgressValueDb = targetGapReductionDb;
    else
        aclr1ProgressValueDb = aclr1ImprovementDb;
    end
    aclr1Progress = aclr1ProgressValueDb > settings.minimumAclr1ProgressDb + settings.comparisonTolerance;
    evmProgress = evmReductionPercent > settings.minimumEvmProgressPercent + settings.comparisonTolerance;
    progressAccepted = aclr1Progress || evmProgress;

    comparison.current = current;
    comparison.trial = trial;
    comparison.settings = settings;
    comparison.aclr1ImprovementDb = aclr1ImprovementDb;
    comparison.aclr2ImprovementDb = aclr2ImprovementDb;
    comparison.evmReductionPercent = evmReductionPercent;
    comparison.currentTargetGapDb = currentTargetGapDb;
    comparison.trialTargetGapDb = trialTargetGapDb;
    comparison.targetGapReductionDb = targetGapReductionDb;
    comparison.aclr1ProgressValueDb = aclr1ProgressValueDb;
    comparison.currentScore = currentScore;
    comparison.trialScore = trialScore;
    comparison.lexicographicallyBetter = lexicographicallyBetter;
    comparison.finiteAccepted = finiteAccepted;
    comparison.safetyAccepted = safetyAccepted;
    comparison.evmAccepted = evmAccepted;
    comparison.aclr1Accepted = aclr1Accepted;
    comparison.aclr2Accepted = aclr2Accepted;
    comparison.aclr1Progress = aclr1Progress;
    comparison.evmProgress = evmProgress;
    comparison.progressAccepted = progressAccepted;

    accepted = false;
    if ~current.allFinite
        reason = 'current_nonfinite';
    elseif ~current.isSafe
        reason = 'current_unsafe';
    elseif ~trial.allFinite
        reason = 'trial_nonfinite';
    elseif ~trial.isSafe
        reason = 'trial_unsafe';
    elseif ~evmAccepted
        reason = 'evm_increased';
    elseif ~aclr1Accepted
        reason = 'aclr1_decreased';
    elseif ~aclr2Accepted
        reason = 'aclr2_degraded';
    elseif ~progressAccepted
        reason = 'no_aclr1_or_evm_progress';
    else
        accepted = true;
        if aclr1Progress && evmProgress
            reason = 'accepted_aclr1_and_evm_progress';
        elseif aclr1Progress
            reason = 'accepted_aclr1_progress';
        else
            reason = 'accepted_evm_progress';
        end
    end
    comparison.accepted = accepted;
    comparison.reason = reason;
end

function metrics = readMetrics(source, sourceName)
    requiredFields = {'worstAclr1Db', 'worstAclr2Db', 'evmRmsPercent', 'isSafe'};
    for fieldIndex = 1:length(requiredFields)
        fieldName = requiredFields{fieldIndex};
        if ~isfield(source, fieldName)
            error('dpdAdaptationGate:MissingMetric', '%s must contain %s.', sourceName, fieldName);
        end
    end

    metrics.worstAclr1Db = numericScalar(source.worstAclr1Db, sourceName, 'worstAclr1Db');
    metrics.worstAclr2Db = numericScalar(source.worstAclr2Db, sourceName, 'worstAclr2Db');
    metrics.evmRmsPercent = numericScalar(source.evmRmsPercent, sourceName, 'evmRmsPercent');
    metrics.isSafe = logicalScalar(source.isSafe, sourceName, 'isSafe');
    metrics.allFinite = isfinite(metrics.worstAclr1Db) && isfinite(metrics.worstAclr2Db) && isfinite(metrics.evmRmsPercent) && metrics.evmRmsPercent >= 0;
end

function settings = readSettings(gateCfg, currentMetrics, trialMetrics)
    settings.maximumEvmIncreasePercent = configurationValue(gateCfg, 'maximumEvmIncreasePercent', 0);
    settings.maximumAclr1DecreaseDb = configurationValue(gateCfg, 'maximumAclr1DecreaseDb', 0);
    settings.maximumAclr2DecreaseDb = configurationValue(gateCfg, 'maximumAclr2DecreaseDb', 0.10);
    settings.minimumAclr1ProgressDb = configurationValue(gateCfg, 'minimumAclr1ProgressDb', 0.02);
    settings.minimumEvmProgressPercent = configurationValue(gateCfg, 'minimumEvmProgressPercent', 0.001);
    settings.comparisonTolerance = configurationValue(gateCfg, 'comparisonTolerance', 1e-9);
    settings.prePaAclr1Db = resolvePrePaAclr1Db(gateCfg, currentMetrics, trialMetrics);

    nonnegativeFields = {'maximumEvmIncreasePercent', 'maximumAclr1DecreaseDb', 'maximumAclr2DecreaseDb', 'minimumAclr1ProgressDb', 'minimumEvmProgressPercent', 'comparisonTolerance'};
    for fieldIndex = 1:length(nonnegativeFields)
        fieldName = nonnegativeFields{fieldIndex};
        validateattributes(settings.(fieldName), {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    end
    if ~isnan(settings.prePaAclr1Db)
        validateattributes(settings.prePaAclr1Db, {'numeric'}, {'scalar', 'real', 'finite'});
    end
end

function prePaAclr1Db = resolvePrePaAclr1Db(gateCfg, currentMetrics, trialMetrics)
    if isfield(gateCfg, 'prePaAclr1Db')
        prePaAclr1Db = gateCfg.prePaAclr1Db;
    elseif isfield(gateCfg, 'targetAclr1Db')
        prePaAclr1Db = gateCfg.targetAclr1Db;
    elseif isfield(currentMetrics, 'prePaAclr1Db')
        prePaAclr1Db = currentMetrics.prePaAclr1Db;
    elseif isfield(currentMetrics, 'prePaWorstAclr1Db')
        prePaAclr1Db = currentMetrics.prePaWorstAclr1Db;
    elseif isfield(trialMetrics, 'prePaAclr1Db')
        prePaAclr1Db = trialMetrics.prePaAclr1Db;
    elseif isfield(trialMetrics, 'prePaWorstAclr1Db')
        prePaAclr1Db = trialMetrics.prePaWorstAclr1Db;
    else
        prePaAclr1Db = NaN;
    end
    if ~isnumeric(prePaAclr1Db) || ~isscalar(prePaAclr1Db) || ~isreal(prePaAclr1Db)
        error('dpdAdaptationGate:InvalidPrePaAclr1', 'prePaAclr1Db must be a real numeric scalar or NaN.');
    end
end

function value = numericScalar(value, sourceName, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
        error('dpdAdaptationGate:InvalidMetric', '%s.%s must be a real numeric scalar.', sourceName, fieldName);
    end
end

function value = logicalScalar(value, sourceName, fieldName)
    if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || (value ~= 0 && value ~= 1)
        error('dpdAdaptationGate:InvalidMetric', '%s.%s must be a logical scalar.', sourceName, fieldName);
    end
    value = logical(value);
end

function greater = lexicographicGreater(leftScore, rightScore, tolerance)
    greater = false;
    for scoreIndex = 1:length(leftScore)
        difference = leftScore(scoreIndex) - rightScore(scoreIndex);
        if difference > tolerance
            greater = true;
            return;
        end
        if difference < -tolerance
            return;
        end
    end
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
