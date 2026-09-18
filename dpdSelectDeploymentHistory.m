function [learningInfo, selectionInfo, selectedCoefficients] = dpdSelectDeploymentHistory(learningInfo, rawMetrics, dpdCfg, paCfg, monitoringSignal)

if nargin < 5
    error('dpdSelectDeploymentHistory:MissingInput', 'learningInfo, rawMetrics, dpdCfg, paCfg and monitoringSignal are required.');
end
if ~isfield(learningInfo, 'coefficientHistory') || ~isfield(learningInfo, 'blockHistory')
    error('dpdSelectDeploymentHistory:MissingLearningHistory', 'learningInfo must contain coefficientHistory and blockHistory.');
end
if ~isstruct(monitoringSignal) || ~isfield(monitoringSignal, 'stressInput')
    error('dpdSelectDeploymentHistory:MissingStressInput', 'monitoringSignal must contain stressInput.');
end
requiredMetricFields = {'updateIndex', 'coefficientScale', 'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db', 'prePaWorstAclr1Db'};
for fieldIndex = 1:length(requiredMetricFields)
    if ~isfield(rawMetrics, requiredMetricFields{fieldIndex})
        error('dpdSelectDeploymentHistory:MissingMetric', 'rawMetrics must contain %s.', requiredMetricFields{fieldIndex});
    end
end

rawHistory = learningInfo.coefficientHistory;
snapshotCount = size(rawHistory, 2);
if snapshotCount < 2 || size(rawHistory, 1) ~= numel(dpdCfg.coefficients)
    error('dpdSelectDeploymentHistory:InvalidCoefficientHistory', 'The coefficient history has an invalid size.');
end
if any(~isfinite(rawHistory(:)))
    error('dpdSelectDeploymentHistory:NonfiniteCoefficientHistory', 'The coefficient history must contain only finite values.');
end
if length(rawMetrics.updateIndex) ~= snapshotCount || length(rawMetrics.coefficientScale) ~= snapshotCount || length(rawMetrics.evmRmsPercent) ~= snapshotCount || length(rawMetrics.worstAclr1Db) ~= snapshotCount || length(rawMetrics.worstAclr2Db) ~= snapshotCount
    error('dpdSelectDeploymentHistory:InvalidMetricHistory', 'The metric history must contain one value per coefficient snapshot.');
end
updateIndex = rawMetrics.updateIndex(:);
if any(~isfinite(updateIndex)) || any(updateIndex < 0) || any(updateIndex ~= round(updateIndex)) || updateIndex(1) ~= 0 || any(diff(updateIndex) <= 0)
    error('dpdSelectDeploymentHistory:InvalidUpdateIndex', 'The update indices must start at zero and increase strictly.');
end
if ~isnumeric(rawMetrics.prePaWorstAclr1Db) || ~isscalar(rawMetrics.prePaWorstAclr1Db) || ~isreal(rawMetrics.prePaWorstAclr1Db) || ~isfinite(rawMetrics.prePaWorstAclr1Db)
    error('dpdSelectDeploymentHistory:InvalidPrePaAclr1', 'prePaWorstAclr1Db must be a finite real scalar.');
end

coefficientScale = rawMetrics.coefficientScale(:).';
if any(~isfinite(coefficientScale)) || any(coefficientScale <= 0)
    error('dpdSelectDeploymentHistory:InvalidCoefficientScale', 'Every coefficient scale must be finite and positive.');
end
normalizedHistory = bsxfun(@times, rawHistory, coefficientScale);
coefficientMask = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);
normalizedHistory(~repmat(coefficientMask(:), 1, snapshotCount)) = 0;

gateCfg.maximumEvmIncreasePercent = configurationValue(dpdCfg, 'adaptationGateMaximumEvmIncreasePercent', 0);
gateCfg.maximumAclr1DecreaseDb = configurationValue(dpdCfg, 'adaptationGateMaximumAclr1DecreaseDb', 0);
gateCfg.maximumAclr2DecreaseDb = configurationValue(dpdCfg, 'adaptationGateMaximumAclr2DecreaseDb', 0.10);
gateCfg.minimumAclr1ProgressDb = configurationValue(dpdCfg, 'adaptationGateMinimumAclr1ProgressDb', 0.02);
gateCfg.minimumEvmProgressPercent = configurationValue(dpdCfg, 'adaptationGateMinimumEvmProgressPercent', 0.001);
gateCfg.comparisonTolerance = configurationValue(dpdCfg, 'adaptationGateComparisonTolerance', 1e-9);
configuredTargetAclr1Db = configurationValue(dpdCfg, 'adaptationTargetAclr1Db', 50);
gateCfg.targetAclr1Db = min(configuredTargetAclr1Db, rawMetrics.prePaWorstAclr1Db);

selectedHistory = complex(zeros(size(normalizedHistory)));
selectedSourceSnapshot = ones(snapshotCount, 1);
deploymentAccepted = false(snapshotCount - 1, 1);
deploymentReason = cell(snapshotCount - 1, 1);
stressSafe = false(snapshotCount, 1);
stressOutputPeak = NaN(snapshotCount, 1);
currentSourceSnapshot = 1;
currentCfg = snapshotConfiguration(dpdCfg, normalizedHistory(:, currentSourceSnapshot));
currentGuard = dpdStressGuard(currentCfg, paCfg, monitoringSignal.stressInput);
currentMetrics = snapshotMetrics(rawMetrics, currentSourceSnapshot, currentGuard.isSafe);
stressSafe(currentSourceSnapshot) = currentGuard.isSafe;
stressOutputPeak(currentSourceSnapshot) = currentGuard.outputPeakMagnitude;
if ~currentMetrics.isSafe
    error('dpdSelectDeploymentHistory:UnsafeInitialSnapshot', 'The initial DPD snapshot is unsafe on the fixed stress monitor.');
end
selectedHistory(:, 1) = normalizedHistory(:, currentSourceSnapshot);

for snapshotIndex = 2:snapshotCount
    trialCfg = snapshotConfiguration(dpdCfg, normalizedHistory(:, snapshotIndex));
    trialGuard = dpdStressGuard(trialCfg, paCfg, monitoringSignal.stressInput);
    trialMetrics = snapshotMetrics(rawMetrics, snapshotIndex, trialGuard.isSafe);
    stressSafe(snapshotIndex) = trialGuard.isSafe;
    stressOutputPeak(snapshotIndex) = trialGuard.outputPeakMagnitude;
    [accepted, reason] = dpdAdaptationGate(currentMetrics, trialMetrics, gateCfg);
    if accepted
        currentSourceSnapshot = snapshotIndex;
        currentMetrics = trialMetrics;
        deploymentAccepted(snapshotIndex - 1) = true;
    end
    deploymentReason{snapshotIndex - 1} = reason;
    selectedSourceSnapshot(snapshotIndex) = currentSourceSnapshot;
    selectedHistory(:, snapshotIndex) = normalizedHistory(:, currentSourceSnapshot);
end

selectedEvmRmsPercent = rawMetrics.evmRmsPercent(selectedSourceSnapshot);
selectedAclr1Db = rawMetrics.worstAclr1Db(selectedSourceSnapshot);
selectedAclr2Db = rawMetrics.worstAclr2Db(selectedSourceSnapshot);
verifySelectedMetrics(selectedEvmRmsPercent, selectedAclr1Db, selectedAclr2Db, gateCfg);

learningInfo.rawCoefficientHistory = rawHistory;
learningInfo.rawBlockHistory = learningInfo.blockHistory;
if isfield(learningInfo, 'acceptedUpdateCount')
    learningInfo.rawAcceptedUpdateCount = learningInfo.acceptedUpdateCount;
end
learningInfo.normalizedRawCoefficientHistory = normalizedHistory;
learningInfo.coefficientHistory = selectedHistory;
learningInfo.coefficientHistoryNormalized = true;
learningInfo.deploymentHistorySelected = true;
learningInfo.deploymentTargetAclr1Db = gateCfg.targetAclr1Db;
learningInfo.deploymentAcceptedUpdateCount = nnz(deploymentAccepted);
learningInfo.acceptedUpdateCount = nnz(deploymentAccepted);
learningInfo.blockHistory.deploymentAccepted = deploymentAccepted;
learningInfo.blockHistory.deploymentSourceUpdateIndex = updateIndex(selectedSourceSnapshot(2:end));
learningInfo.blockHistory.deploymentReason = deploymentReason;
if isfield(learningInfo.blockHistory, 'updateAccepted')
    learningInfo.blockHistory.updateAccepted = deploymentAccepted;
end
if isfield(learningInfo.blockHistory, 'accepted')
    learningInfo.blockHistory.accepted = deploymentAccepted;
end

selectedCoefficients = reshape(selectedHistory(:, end), size(dpdCfg.coefficients));
selectionInfo.enabled = true;
selectionInfo.targetAclr1Db = gateCfg.targetAclr1Db;
selectionInfo.prePaAclr1Db = rawMetrics.prePaWorstAclr1Db;
selectionInfo.selectedSourceSnapshotIndex = selectedSourceSnapshot;
selectionInfo.selectedSourceUpdateIndex = updateIndex(selectedSourceSnapshot);
selectionInfo.deploymentAccepted = deploymentAccepted;
selectionInfo.deploymentReason = deploymentReason;
selectionInfo.stressSafe = stressSafe;
selectionInfo.stressOutputPeak = stressOutputPeak;
selectionInfo.normalizationScale = coefficientScale;
selectionInfo.selectedEvmRmsPercent = selectedEvmRmsPercent;
selectionInfo.selectedAclr1Db = selectedAclr1Db;
selectionInfo.selectedAclr2Db = selectedAclr2Db;
selectionInfo.finalSourceUpdateIndex = updateIndex(selectedSourceSnapshot(end));
selectionInfo.finalEvmRmsPercent = currentMetrics.evmRmsPercent;
selectionInfo.finalAclr1Db = currentMetrics.worstAclr1Db;
selectionInfo.finalAclr2Db = currentMetrics.worstAclr2Db;
selectionInfo.finalAclr1ShortfallDb = rawMetrics.prePaWorstAclr1Db - currentMetrics.worstAclr1Db;

end

function verifySelectedMetrics(evmRmsPercent, aclr1Db, aclr2Db, gateCfg)
tolerance = gateCfg.comparisonTolerance;
if any(diff(evmRmsPercent) > gateCfg.maximumEvmIncreasePercent + tolerance)
    error('dpdSelectDeploymentHistory:EvmInvariantViolation', 'The selected deployment history increases RMS EVM.');
end
if any(diff(aclr1Db) < -gateCfg.maximumAclr1DecreaseDb - tolerance)
    error('dpdSelectDeploymentHistory:Aclr1InvariantViolation', 'The selected deployment history decreases ACLR1.');
end
if any(diff(aclr2Db) < -gateCfg.maximumAclr2DecreaseDb - tolerance)
    error('dpdSelectDeploymentHistory:Aclr2InvariantViolation', 'The selected deployment history exceeds the ACLR2 tolerance.');
end
end

function metrics = snapshotMetrics(rawMetrics, snapshotIndex, isSafe)
metrics.worstAclr1Db = rawMetrics.worstAclr1Db(snapshotIndex);
metrics.worstAclr2Db = rawMetrics.worstAclr2Db(snapshotIndex);
metrics.evmRmsPercent = rawMetrics.evmRmsPercent(snapshotIndex);
metrics.prePaAclr1Db = rawMetrics.prePaWorstAclr1Db;
metrics.isSafe = logical(isSafe) && isfinite(metrics.worstAclr1Db) && isfinite(metrics.worstAclr2Db) && isfinite(metrics.evmRmsPercent);
end

function snapshotCfg = snapshotConfiguration(dpdCfg, coefficientVector)
snapshotCfg = dpdCfg;
snapshotCfg.enabled = true;
snapshotCfg.modelLoaded = false;
snapshotCfg.modelStale = false;
snapshotCfg.modelStaleReason = '';
snapshotCfg.coefficients = reshape(coefficientVector, size(dpdCfg.coefficients));
end

function value = configurationValue(configuration, fieldName, defaultValue)
if isfield(configuration, fieldName)
    value = configuration.(fieldName);
else
    value = defaultValue;
end
end
