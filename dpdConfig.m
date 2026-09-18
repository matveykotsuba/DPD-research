function dpdCfg = dpdConfig(paCfg, loadModel, learningArchitecture, cfrCfg)

configuredEnabled = true;
configuredAdaptationEnabled = true;
configuredLearningArchitecture = 'direct';

if nargin < 1
    paCfg = [];
end

if nargin < 2
    loadModel = true;
end

if nargin < 3 || isempty(learningArchitecture)
    learningArchitecture = configuredLearningArchitecture;
end

if nargin < 4 || isempty(cfrCfg)
    cfrCfg = cfrConfig(false);
end

learningArchitecture = validateLearningArchitecture( learningArchitecture);

dpdCfg.enabled = configuredEnabled;
dpdCfg.adaptationEnabled = configuredAdaptationEnabled;
dpdCfg.learningArchitecture = learningArchitecture;
dpdCfg.cfrSignature = cfrSignature(cfrCfg);

dpdCfg.orders = 1:2:7;
dpdCfg.signalDelays = 0:4;
dpdCfg.envelopeDelays = 0:4;
dpdCfg.diagonalCount = 4;

if ~isempty(paCfg) && isfield(paCfg, 'modelSampleRate')
    dpdCfg.modelSampleRate = paCfg.modelSampleRate;
else
    dpdCfg.modelSampleRate = 737.28e6;
end

dpdCfg.coefficients = identityCoefficients(dpdCfg);

dpdCfg.blockLength = 30000;
dpdCfg.adaptationBlockLength = dpdCfg.blockLength;

dpdCfg.maximumCoefficientNorm = 1e3;
dpdCfg.feedbackAlignmentMaximumDelay = 64;
dpdCfg.feedbackRemoveDc = true;
dpdCfg.minimumFeedbackSamples = 128;

dpdCfg.regularizationLambdaValues = logspace(-12, 2, 200);
dpdCfg.indirectCoefficientBlend = 1;
dpdCfg.indirectLineSearchFactors = [1 0.5 0.25 0.125 0.0625];
dpdCfg.indirectMinimumMonitorImprovementDb = 0;
dpdCfg.indirectMonitorSampleCount = dpdCfg.blockLength;
dpdCfg.indirectRegularizationCenter = 'zero';
dpdCfg.learningInitialization = 'identity';
dpdCfg.directMaximumEpochs = 1;
dpdCfg.directLambdaValues = logspace(-6, 1, 8);
dpdCfg.directLineSearchFactors = [1 0.5 0.25 0.125 0.0625];
dpdCfg.directJacobianStride = 8;
dpdCfg.directJacobianRowsPerChunk = 64;
dpdCfg.directMonitorSampleCount = dpdCfg.blockLength;
dpdCfg.directMinimumMonitorImprovementDb = 0;
dpdCfg.directSelectionObjective = 'nmse';
dpdCfg.directTrainingUsesSpectralObjective = true;
dpdCfg.directSpectralAclr1Weight = 0;
dpdCfg.directSpectralAclr2Weight = 0;
dpdCfg.directMaximumMonitorNmseDegradationDb = 0;
dpdCfg.directMinimumMonitorAclr1ImprovementDb = 0.02;
dpdCfg.directMaximumMonitorAclr2DegradationDb = 0.10;
dpdCfg.directAclrSpectrumSegmentLength = 8192;
dpdCfg.directAclrCandidateCount = 10;
dpdCfg.directAclrMonitorSampleCount = 262144;
dpdCfg.adaptationMetricsEnabled = true;
dpdCfg.adaptationMetricSeed = 7;
dpdCfg.adaptationMetricOfdmSymbols = 140;
dpdCfg.adaptationMetricSpectrumSegmentLength = 8192;
dpdCfg.adaptationMetricChunkSamples = 8192;
dpdCfg.adaptationGateEnabled = true;
dpdCfg.adaptationGateCandidateCount = 3;
dpdCfg.adaptationGateMaximumEvmIncreasePercent = 0;
dpdCfg.adaptationGateMaximumAclr1DecreaseDb = 0;
dpdCfg.adaptationGateMaximumAclr2DecreaseDb = 0.10;
dpdCfg.adaptationGateMinimumAclr1ProgressDb = 0.02;
dpdCfg.adaptationGateMinimumEvmProgressPercent = 0.001;
dpdCfg.adaptationGateComparisonTolerance = 1e-9;
dpdCfg.adaptationTargetAclr1Db = 50;
dpdCfg.adaptationTargetAclr1GapDb = 1;
dpdCfg.architectureComparisonPlotEnabled = true;
dpdCfg.architectureComparisonAclr1ReferenceDb = 52;

dpdCfg.trainingSeed = 11;
dpdCfg.validationSeed = 12;
dpdCfg.testSeed = 13;
dpdCfg.learningMonitorSeed = 14;
dpdCfg.learningMonitorOfdmSymbols = 28;
dpdCfg.trainingSampleCount = 300000;
dpdCfg.validationSampleCount = 131072;
dpdCfg.testSampleCount = 131072;
dpdCfg.learningMonitorSampleCount = 131072;
dpdCfg.learningMonitorStressHeadroomDb = 0.5;
dpdCfg.maximumDpdPeakGainDb = 12;
dpdCfg.minimumValidationImprovementDb = 0;
dpdCfg.minimumTestImprovementDb = 0;
dpdCfg.maximumDpdAveragePowerChangeDb = 0.25;
dpdCfg.maximumPaOutputPowerMismatchDb = 0.25;
dpdCfg.maximumAclr1DegradationDb = 0.10;
dpdCfg.maximumAclr2DegradationDb = 7.00;
dpdCfg.minimumAclr2Db = 45.00;
dpdCfg.incumbentMinimumValidationAclr1ImprovementDb = 0.02;
dpdCfg.incumbentComparisonTolerance = 1e-9;
dpdCfg.saveModel = true;

if ~isempty(paCfg) && isfield(paCfg, 'maximumInputMagnitude')
    dpdCfg.maximumOutputMagnitude = paCfg.maximumInputMagnitude;
else
    dpdCfg.maximumOutputMagnitude = Inf;
end

if ~isempty(paCfg)
    dpdCfg.paModelSignature = dpdPaModelSignature(paCfg);
end

dpdCfg.comparison.seed = 19;
dpdCfg.comparison.equalOutputPower = true;
dpdCfg.comparison.outputPowerToleranceDb = 0.05;
dpdCfg.comparison.minimumDriveScaleDb = -6;
dpdCfg.comparison.maximumDriveScaleDb = 6;
dpdCfg.comparison.powerMatchIterations = 14;
dpdCfg.comparison.powerMatchProbeSamples = 65536;
dpdCfg.comparison.spectrumMaximumSamples = 1048576;
dpdCfg.comparison.spectrumSegmentLength = 8192;
dpdCfg.comparison.amAmMaximumSamples = 200000;
dpdCfg.comparison.amAmBinCount = 180;

dpdCfg.modelFile = dpdArchitectureArtifactFile(dpdCfg.learningArchitecture, dpdCfg.cfrSignature, 'model');
dpdCfg.modelLoaded = false;
dpdCfg.loadedModelArchitecture = '';
dpdCfg.modelStale = false;
dpdCfg.modelStaleReason = '';

if loadModel && exist(dpdCfg.modelFile, 'file') == 2
    savedModel = load(dpdCfg.modelFile);
    [dpdCfg, modelAccepted] = applySavedModel(dpdCfg, savedModel, paCfg);
    dpdCfg.modelLoaded = modelAccepted;
end

if ~isempty(paCfg) && isfield(paCfg, 'signalDelays') && isfield(paCfg, 'envelopeDelays')
    dpdCfg.metricFirstValidSample = 1 + max([dpdCfg.signalDelays(:); dpdCfg.envelopeDelays(:)]) + max([paCfg.signalDelays(:); paCfg.envelopeDelays(:)]);
end

[dpdCfg.coefficientMask, dpdCfg.activeCoefficientIndices, maskInfo] = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);

if ~isequal(size(dpdCfg.coefficients), size(dpdCfg.coefficientMask))
    error('dpdConfig:InvalidSavedCoefficientSize', 'The DPD coefficient matrix is inconsistent with its GMP model.');
end

dpdCfg.coefficients(~dpdCfg.coefficientMask) = 0;
dpdCfg.numberOfActiveCoefficients = maskInfo.numberOfActiveCoefficients;
end

function coefficients = identityCoefficients(modelCfg)
numberOfColumns = 1 + (length(modelCfg.orders) - 1) * length(modelCfg.envelopeDelays);
coefficients = complex(zeros( length(modelCfg.signalDelays), numberOfColumns));
coefficients(1, 1) = 1;
end

function [dpdCfg, modelAccepted] = applySavedModel(dpdCfg, savedModel, paCfg)
modelAccepted = false;
if ~isfield(savedModel, 'coefficients')
    dpdCfg.modelStale = true;
    dpdCfg.modelStaleReason = 'The saved DPD model has no coefficients and must be retrained.';
    return;
end

if ~isfield(savedModel, 'fitInfo') || ~isstruct(savedModel.fitInfo) || ~isfield(savedModel.fitInfo, 'learningArchitecture')
    dpdCfg.modelStale = true;
    dpdCfg.modelStaleReason = 'The saved DPD model has no learning architecture and must be retrained.';
    return;
end

fitInfo = savedModel.fitInfo;
savedLearningArchitecture = validateLearningArchitecture( fitInfo.learningArchitecture);
if ~strcmp(savedLearningArchitecture, dpdCfg.learningArchitecture)
    dpdCfg.modelStale = true;
    dpdCfg.modelStaleReason = 'The saved DPD model belongs to another learning architecture and must be retrained.';
    return;
end

structureFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount'};
for fieldIndex = 1:length(structureFields)
    fieldName = structureFields{fieldIndex};
    if ~isfield(fitInfo, fieldName) || ~isequal(fitInfo.(fieldName), dpdCfg.(fieldName))
        dpdCfg.modelStale = true;
        dpdCfg.modelStaleReason = 'The saved DPD model has a different GMP structure and must be retrained.';
        return;
    end
end

if ~isempty(paCfg)
    if ~isfield(fitInfo, 'paModelSignature')
        dpdCfg.modelStale = true;
        dpdCfg.modelStaleReason = 'The saved DPD model has no PA model signature and must be retrained.';
        return;
    end
    currentPaModelSignature = dpdPaModelSignature(paCfg);
    if ~isequaln(fitInfo.paModelSignature, currentPaModelSignature)
        dpdCfg.modelStale = true;
        dpdCfg.modelStaleReason = 'The saved DPD model was trained for a different PA model or operating point.';
        return;
    end
end

if isfield(fitInfo, 'cfrSignature')
    savedCfrSignature = fitInfo.cfrSignature;
else
    savedCfrSignature = cfrSignature(cfrConfig(false));
end
if ~isequaln(savedCfrSignature, dpdCfg.cfrSignature)
    dpdCfg.modelStale = true;
    dpdCfg.modelStaleReason = 'The saved DPD model was trained with a different CFR configuration.';
    return;
end

expectedColumns = 1 + (length(dpdCfg.orders) - 1) * length(dpdCfg.envelopeDelays);
expectedSize = [length(dpdCfg.signalDelays), expectedColumns];
if ~isnumeric(savedModel.coefficients) || ~isequal(size(savedModel.coefficients), expectedSize) || any(~isfinite(savedModel.coefficients(:)))
    dpdCfg.modelStale = true;
    dpdCfg.modelStaleReason = 'The saved DPD coefficients are invalid and must be retrained.';
    return;
end

modelFields = {'modelSampleRate', 'maximumOutputMagnitude'};

for fieldIndex = 1:length(modelFields)
    fieldName = modelFields{fieldIndex};
    if isfield(fitInfo, fieldName)
        dpdCfg.(fieldName) = fitInfo.(fieldName);
    end
end

if isfield(fitInfo, 'referenceLinearGain')
    dpdCfg.metricReferenceLinearGain = fitInfo.referenceLinearGain;
end

dpdCfg.coefficients = savedModel.coefficients;

if isfield(savedModel, 'trainingInfo')
    dpdCfg.savedTrainingInfo = savedModel.trainingInfo;
end
dpdCfg.modelStale = false;
dpdCfg.modelStaleReason = '';
dpdCfg.loadedModelArchitecture = savedLearningArchitecture;
modelAccepted = true;
end

function learningArchitecture = validateLearningArchitecture( learningArchitecture)
if ~ischar(learningArchitecture)
    error('dpdConfig:InvalidLearningArchitecture', 'learningArchitecture must be a character vector.');
end

learningArchitecture = lower(strtrim(learningArchitecture));
if ~strcmp(learningArchitecture, 'indirect') && ~strcmp(learningArchitecture, 'direct')
    error('dpdConfig:UnknownLearningArchitecture', ['learningArchitecture must be either ''indirect'' or ' '''direct''.']);
end
end
