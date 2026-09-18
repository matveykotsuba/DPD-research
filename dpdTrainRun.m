function [dpdCfg, trainingInfo] = dpdTrainRun(cfg)

    if nargin < 1 || isempty(cfg)
        cfg = config();
    end
    if ~isstruct(cfg) || ~isfield(cfg, 'pa')
        error('dpdTrainRun:MissingPaConfiguration', 'cfg.pa is required.');
    end
    cfrCfg = systemCfrConfiguration(cfg);
    if ~isfield(cfg, 'dpd') || isempty(cfg.dpd)
        cfg.dpd = dpdConfig(cfg.pa, false, 'indirect', cfrCfg);
    end

    paCfg = cfg.pa;
    dpdCfg = prepareDpdConfiguration(cfg.dpd, paCfg, cfrCfg);
    if strcmp(dpdCfg.learningArchitecture, 'direct') && logical(configurationValue(dpdCfg, 'directTrainingUsesSpectralObjective', true))
        dpdCfg.directSelectionObjective = 'aclr';
        dpdCfg.directAclrChannelBandwidth = cfg.channelBandwidth;
        dpdCfg.directAclrMeasurementBandwidth = cfg.aclrMeasurementBandwidth;
    end
    incumbentAvailable = compatibleLoadedIncumbent(dpdCfg);
    incumbentComparisonMetrics = struct();
    incumbentComparisonAvailable = false;
    incumbentComparisonReason = 'compatible incumbent identity-start history is unavailable';
    if incumbentAvailable
        incumbentCfg = dpdCfg;
        if isfield(incumbentCfg, 'savedTrainingInfo')
            [incumbentComparisonMetrics, incumbentComparisonAvailable, incumbentComparisonReason] = dpdMakeArchitectureComparisonMetrics(incumbentCfg.savedTrainingInfo);
        end
    else
        incumbentCfg = struct();
    end
    dpdCfg.modelLoaded = false;
    dpdCfg.loadedModelArchitecture = '';
    dpdCfg.modelStale = false;
    dpdCfg.modelStaleReason = '';
    if isfield(dpdCfg, 'savedTrainingInfo')
        dpdCfg = rmfield(dpdCfg, 'savedTrainingInfo');
    end
    options = trainingOptions(dpdCfg);
    options.modelFile = resolveModelFile( options.modelFile, paCfg, options.saveModel);
    validateTrainingOptions(options, dpdCfg, paCfg);
    [cachedComparisonMetrics, ~, cachedComparisonAvailable, cachedComparisonReason, comparisonCacheFile] = dpdLoadArchitectureComparisonCache(dpdCfg, dpdCfg.learningArchitecture);
    if sameWindowsPath(comparisonCacheFile, options.modelFile)
        error('dpdTrainRun:ComparisonCacheModelCollision', 'The architecture comparison cache must not overwrite the DPD model file.');
    end
    if isfield(paCfg, 'modelFile') && sameWindowsPath(comparisonCacheFile, paCfg.modelFile)
        error('dpdTrainRun:ComparisonCachePaCollision', 'The architecture comparison cache must not overwrite the PA model file.');
    end

    fprintf('\nDPD %s learning at %.2f MHz (before channel/AWGN).\n', dpdCfg.learningArchitecture, paCfg.modelSampleRate / 1e6);
    fprintf(['Training record: %d samples, %d blocks of at most ' '%d samples.\n'], options.trainingSampleCount, ceil(options.trainingSampleCount / options.blockLength), options.blockLength);
    fprintf('Train/validation/test seeds: %d/%d/%d.\n', options.trainingSeed, options.validationSeed, options.testSeed);

    [trainingInput, trainingSignalInfo] = dpdTrainingSignal( cfg, options.trainingSeed, options.trainingSampleCount);
    learningMonitoringSignal = dpdMonitoringSignal(cfg, dpdCfg, 'learning');
    adaptationMonitoringSignal = dpdMonitoringSignal(cfg, dpdCfg, 'adaptation');
    if incumbentComparisonAvailable && ~comparisonProtocolMatchesConfiguration(incumbentComparisonMetrics, options, adaptationMonitoringSignal, dpdCfg.learningArchitecture)
        incumbentComparisonAvailable = false;
        incumbentComparisonReason = 'compatible incumbent identity-start history uses a different comparison protocol';
    end
    if cachedComparisonAvailable && ~comparisonProtocolMatchesConfiguration(cachedComparisonMetrics, options, adaptationMonitoringSignal, dpdCfg.learningArchitecture)
        cachedComparisonAvailable = false;
        cachedComparisonReason = 'architecture comparison cache uses a different comparison protocol';
    end
    identityComparisonRebuildRequired = incumbentAvailable && ~incumbentComparisonAvailable && ~cachedComparisonAvailable;
    if incumbentAvailable
        if identityComparisonRebuildRequired
            dpdCfg.learningInitialization = 'identity';
            fprintf('Identity-start comparison history is missing; training starts from identity while the loaded incumbent remains protected.\n');
        else
            dpdCfg.learningInitialization = 'configured';
        end
    end
    learningMonitorInput = learningMonitoringSignal.learningInput;
    stressMonitorInput = learningMonitoringSignal.stressInput;
    learningMonitorSignalInfo = learningMonitoringSignal.info;
    adaptationMetricSignalInfo = adaptationMonitoringSignal.info;
    trainingPaOutput = gmpCore(trainingInput, paCfg);

    simulatedFeedbackCfg = dpdCfg;
    requestedAlignmentMaximumDelay = configurationValue( dpdCfg, 'feedbackAlignmentMaximumDelay', 64);
    simulatedFeedbackCfg.feedbackAlignmentMaximumDelay = 0;
    [~, ~, alignmentInfo] = dpdFeedbackAlign( trainingInput, trainingPaOutput, simulatedFeedbackCfg);
    alignmentInfo.configuredMaximumDelay = requestedAlignmentMaximumDelay;
    alignmentInfo.usedMaximumDelay = 0;
    alignmentInfo.delayWasFixedForSimulatedPa = true;

    feedbackAlignmentGain = alignmentInfo.linearGain;
    commonMemory = maximumModelDelay(dpdCfg) + maximumModelDelay(paCfg);
    gainRows = commonMemory + 1:length(trainingInput);
    gainInput = trainingInput(gainRows);
    gainOutput = trainingPaOutput(gainRows);
    referenceGain = (gainInput' * gainOutput) / sum(abs(gainInput).^2);
    if ~isfinite(referenceGain) || abs(referenceGain) <= realmin
        error('dpdTrainRun:InvalidReferenceGain', 'The common PA-only training gain is zero or nonfinite.');
    end
    alignmentInfo.feedbackAlignmentGain = feedbackAlignmentGain;
    alignmentInfo.linearGain = referenceGain;
    alignmentInfo.gainReference = 'raw PA-only training samples after the common cascade crop';
    dpdCfg.metricReferenceLinearGain = referenceGain;
    dpdCfg.metricFirstValidSample = 1 + commonMemory;
    if incumbentAvailable
        incumbentCfg.metricReferenceLinearGain = referenceGain;
        incumbentCfg.metricFirstValidSample = 1 + commonMemory;
        incumbentCfg.enabled = true;
    end

    switch dpdCfg.learningArchitecture
        case 'indirect'
            [learnedCoefficients, learningDetails] = dpdIndirectLearn(trainingInput, paCfg, dpdCfg, referenceGain, learningMonitorInput, stressMonitorInput);
        case 'direct'
            [learnedCoefficients, learningDetails] = dpdDirectLearn(trainingInput, paCfg, dpdCfg, referenceGain, learningMonitorInput, stressMonitorInput, learningMonitoringSignal.modelInput);
        otherwise
            error('dpdTrainRun:UnknownLearningArchitecture', 'Unsupported DPD learning architecture: %s.', dpdCfg.learningArchitecture);
    end
    requireFiniteCoefficients(learnedCoefficients);
    deploymentMonitoringSignal = adaptationMonitoringSignal;
    deploymentMonitoringSignal.stressInput = stressMonitorInput;
    clear learningMonitoringSignal learningMonitorInput;

    adaptationMetricsEnabled = logical(configurationValue(dpdCfg, 'adaptationMetricsEnabled', true));
    adaptationGateEnabled = logical(configurationValue(dpdCfg, 'adaptationGateEnabled', true));
    rawBlockMetrics = struct();
    deploymentSelection = struct();
    if adaptationMetricsEnabled || adaptationGateEnabled || identityComparisonRebuildRequired
        rawBlockMetrics = dpdBlockMetrics(cfg, trainingInput, learningDetails, paCfg, dpdCfg, adaptationMonitoringSignal);
        if adaptationGateEnabled
            [learningDetails, deploymentSelection, learnedCoefficients] = dpdSelectDeploymentHistory(learningDetails, rawBlockMetrics, dpdCfg, paCfg, deploymentMonitoringSignal);
            requireFiniteCoefficients(learnedCoefficients);
            blockMetrics = dpdBlockMetrics(cfg, trainingInput, learningDetails, paCfg, dpdCfg, adaptationMonitoringSignal);
            verifyDeploymentBlockMetrics(blockMetrics, deploymentSelection, dpdCfg);
        else
            blockMetrics = rawBlockMetrics;
        end
    else
        blockMetrics = struct();
    end

    attemptedComparisonSource = struct();
    attemptedComparisonSource.options = options;
    attemptedComparisonSource.learning = learningDetails;
    attemptedComparisonSource.blockMetrics = blockMetrics;
    [attemptedComparisonMetrics, attemptedComparisonAvailable, attemptedComparisonReason] = dpdMakeArchitectureComparisonMetrics(attemptedComparisonSource);

    candidateCfg = candidateConfiguration( dpdCfg, learnedCoefficients);
    candidatePreparationAccepted = true;
    candidatePreparationReason = 'candidate_preparation_passed';
    if adaptationGateEnabled
        [normalizedCheckCfg, normalizationInfo] = dpdNormalizeCandidate(trainingInput, candidateCfg, paCfg);
        normalizationCheckScale = normalizationInfo.coefficientScale;
        normalizationToleranceDb = configurationValue(dpdCfg, 'deploymentNormalizationToleranceDb', 1e-8);
        validateattributes(normalizationToleranceDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
        if ~normalizationInfo.succeeded || ~isfinite(normalizationCheckScale) || normalizationCheckScale <= 0 || abs(20 * log10(normalizationCheckScale)) > normalizationToleranceDb
            if incumbentAvailable
                candidatePreparationAccepted = false;
                candidatePreparationReason = 'candidate_normalization_mismatch';
            else
                error('dpdTrainRun:DeploymentNormalizationMismatch', 'The selected deployment coefficients do not satisfy the common training-power normalization.');
            end
        end
        normalizationInfo.verificationCoefficientScale = normalizationCheckScale;
        normalizationInfo.coefficientScale = 1;
        normalizationInfo.appliedCoefficientScale = 1;
        normalizationInfo.reference = 'selected deployment history normalized once on the common training signal';
        requireFiniteCoefficients(normalizedCheckCfg.coefficients);
    else
        [candidateCfg, normalizationInfo] = dpdNormalizeCandidate(trainingInput, candidateCfg, paCfg);
        normalizationInfo.appliedCoefficientScale = normalizationInfo.coefficientScale;
    end
    normalizationInfo.stressGuard = dpdStressGuard(candidateCfg, paCfg, stressMonitorInput);
    if ~normalizationInfo.stressGuard.isSafe
        if incumbentAvailable
            candidatePreparationAccepted = false;
            candidatePreparationReason = 'candidate_stress_monitor_unsafe';
        else
            error('dpdTrainRun:UnsafeStressMonitorPeak', 'The normalized DPD candidate exceeds the permitted output peak on the fixed OFDM stress monitor.');
        end
    end

    [validationInput, validationSignalInfo] = dpdTrainingSignal( cfg, options.validationSeed, options.validationSampleCount);
    baselineCfg = identityConfiguration(dpdCfg);
    validationBaseline = dpdCascadeMetrics( validationInput, baselineCfg, paCfg, cfg);
    validationCandidate = dpdCascadeMetrics( validationInput, candidateCfg, paCfg, cfg);
    validationCandidate = addFairnessChecks( validationCandidate, validationBaseline, options, normalizationInfo.coefficientScale, normalizationInfo.succeeded);
    validationCandidateImprovementDb = validationBaseline.nmseDb - validationCandidate.nmseDb;

    printMetrics('Validation PA baseline', validationBaseline);
    printMetrics(['Validation ' dpdCfg.learningArchitecture ' candidate'], validationCandidate);
    if incumbentAvailable
        validationIncumbent = dpdCascadeMetrics(validationInput, incumbentCfg, paCfg, cfg);
        validationIncumbent = addFairnessChecks(validationIncumbent, validationBaseline, options, 1, true);
        validationIncumbentImprovementDb = validationBaseline.nmseDb - validationIncumbent.nmseDb;
        printMetrics('Validation loaded incumbent', validationIncumbent);
        requireQualifiedResult(validationBaseline, validationIncumbent, validationIncumbentImprovementDb, options.minimumValidationImprovementDb, 'incumbent validation');
    else
        validationIncumbent = struct();
        validationIncumbentImprovementDb = NaN;
        requireQualifiedResult(validationBaseline, validationCandidate, validationCandidateImprovementDb, options.minimumValidationImprovementDb, 'validation');
    end

    [testInput, testSignalInfo] = dpdTrainingSignal( cfg, options.testSeed, options.testSampleCount);
    testBaseline = dpdCascadeMetrics( testInput, baselineCfg, paCfg, cfg);
    testCandidate = dpdCascadeMetrics( testInput, candidateCfg, paCfg, cfg);
    testCandidate = addFairnessChecks( testCandidate, testBaseline, options, normalizationInfo.coefficientScale, normalizationInfo.succeeded);
    testCandidateImprovementDb = testBaseline.nmseDb - testCandidate.nmseDb;
    printMetrics('Test PA baseline', testBaseline);
    printMetrics(['Test ' dpdCfg.learningArchitecture ' candidate'], testCandidate);
    if incumbentAvailable
        testIncumbent = dpdCascadeMetrics(testInput, incumbentCfg, paCfg, cfg);
        testIncumbent = addFairnessChecks(testIncumbent, testBaseline, options, 1, true);
        testIncumbentImprovementDb = testBaseline.nmseDb - testIncumbent.nmseDb;
        printMetrics('Test loaded incumbent', testIncumbent);
        requireQualifiedResult(testBaseline, testIncumbent, testIncumbentImprovementDb, options.minimumTestImprovementDb, 'incumbent test');
        [candidateSelected, incumbentDecision] = incumbentReplacementDecision(validationBaseline, validationCandidate, validationCandidateImprovementDb, validationIncumbent, testBaseline, testCandidate, testCandidateImprovementDb, testIncumbent, options, candidatePreparationAccepted, candidatePreparationReason);
    else
        testIncumbent = struct();
        testIncumbentImprovementDb = NaN;
        requireQualifiedResult(testBaseline, testCandidate, testCandidateImprovementDb, options.minimumTestImprovementDb, 'test');
        candidateSelected = true;
        incumbentDecision.available = false;
        incumbentDecision.candidateSelected = true;
        incumbentDecision.reason = 'no_loaded_incumbent';
    end

    if candidateSelected
        dpdCfg = candidateCfg;
        deployedLearningDetails = learningDetails;
        deployedBlockMetrics = blockMetrics;
        deployedRawBlockMetrics = rawBlockMetrics;
        deployedSelection = deploymentSelection;
        deployedNormalizationInfo = normalizationInfo;
        deployedValidation = validationCandidate;
        deployedValidationImprovementDb = validationCandidateImprovementDb;
        deployedTest = testCandidate;
        deployedTestImprovementDb = testCandidateImprovementDb;
        shouldSaveModel = logical(options.saveModel);
        dpdCfg.modelLoaded = shouldSaveModel;
        dpdCfg.loadedModelArchitecture = dpdCfg.learningArchitecture;
        dpdCfg.trained = true;
        dpdCfg.modelVerified = true;
        dpdCfg.modelFile = options.modelFile;
        dpdCfg.modelStale = false;
        dpdCfg.modelStaleReason = '';
    else
        dpdCfg = incumbentCfg;
        deployedLearningDetails = savedTrainingField(incumbentCfg, 'learning');
        deployedBlockMetrics = savedTrainingField(incumbentCfg, 'blockMetrics');
        deployedRawBlockMetrics = savedTrainingField(incumbentCfg, 'rawBlockMetrics');
        deployedSelection = savedTrainingField(incumbentCfg, 'deploymentSelection');
        deployedNormalizationInfo = savedTrainingField(incumbentCfg, 'trainingNormalization');
        deployedValidation = validationIncumbent;
        deployedValidationImprovementDb = validationIncumbentImprovementDb;
        deployedTest = testIncumbent;
        deployedTestImprovementDb = testIncumbentImprovementDb;
        shouldSaveModel = false;
        dpdCfg.enabled = true;
        dpdCfg.modelLoaded = true;
        dpdCfg.loadedModelArchitecture = dpdCfg.learningArchitecture;
        dpdCfg.modelVerified = true;
        dpdCfg.modelStale = false;
        dpdCfg.modelStaleReason = '';
    end
    [coefficientMask, activeIndices, maskInfo] = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);
    dpdCfg.coefficientMask = coefficientMask;
    dpdCfg.activeCoefficientIndices = activeIndices;
    dpdCfg.numberOfActiveCoefficients = maskInfo.numberOfActiveCoefficients;

    deployedComparisonSource.options = options;
    deployedComparisonSource.learning = deployedLearningDetails;
    deployedComparisonSource.blockMetrics = deployedBlockMetrics;
    [deployedComparisonMetrics, deployedComparisonAvailable, deployedComparisonReason] = dpdMakeArchitectureComparisonMetrics(deployedComparisonSource);
    if attemptedComparisonAvailable
        architectureComparisonMetrics = attemptedComparisonMetrics;
        architectureComparisonAvailable = true;
        architectureComparisonReason = attemptedComparisonReason;
    elseif incumbentComparisonAvailable
        architectureComparisonMetrics = incumbentComparisonMetrics;
        architectureComparisonAvailable = true;
        architectureComparisonReason = incumbentComparisonReason;
    elseif cachedComparisonAvailable
        architectureComparisonMetrics = cachedComparisonMetrics;
        architectureComparisonAvailable = true;
        architectureComparisonReason = cachedComparisonReason;
    else
        architectureComparisonMetrics = deployedComparisonMetrics;
        architectureComparisonAvailable = deployedComparisonAvailable;
        architectureComparisonReason = deployedComparisonReason;
    end

    trainingInfo.learningArchitecture = dpdCfg.learningArchitecture;
    trainingInfo.algorithm = algorithmDescription( dpdCfg.learningArchitecture);
    trainingInfo.createdAt = datestr(now, 30);
    trainingInfo.options = options;
    trainingInfo.trainingSignal = trainingSignalInfo;
    trainingInfo.learningMonitorSignal = learningMonitorSignalInfo;
    trainingInfo.adaptationMetricSignal = adaptationMetricSignalInfo;
    trainingInfo.validationSignal = validationSignalInfo;
    trainingInfo.testSignal = testSignalInfo;
    trainingInfo.feedbackAlignment = alignmentInfo;
    trainingInfo.learning = deployedLearningDetails;
    trainingInfo.blockMetrics = deployedBlockMetrics;
    trainingInfo.rawBlockMetrics = deployedRawBlockMetrics;
    trainingInfo.deploymentSelection = deployedSelection;
    trainingInfo.trainingNormalization = deployedNormalizationInfo;
    trainingInfo.referenceGain = referenceGain;
    trainingInfo.validation.baseline = validationBaseline;
    trainingInfo.validation.candidate = deployedValidation;
    trainingInfo.validation.improvementDb = deployedValidationImprovementDb;
    trainingInfo.test.baseline = testBaseline;
    trainingInfo.test.candidate = deployedTest;
    trainingInfo.test.improvementDb = deployedTestImprovementDb;
    trainingInfo.test.role = 'Deployment qualification holdout; not used for learning.';
    trainingInfo.modelFile = dpdCfg.modelFile;
    trainingInfo.cfrSignature = dpdCfg.cfrSignature;
    trainingInfo.saved = shouldSaveModel;
    trainingInfo.incumbentDecision = incumbentDecision;
    trainingInfo.incumbentAvailable = incumbentAvailable;
    trainingInfo.architectureComparisonMetrics = architectureComparisonMetrics;
    trainingInfo.architectureComparisonAvailable = architectureComparisonAvailable;
    trainingInfo.architectureComparisonReason = architectureComparisonReason;
    trainingInfo.identityComparisonRebuildRequired = identityComparisonRebuildRequired;
    if incumbentAvailable
        trainingInfo.incumbent.validation = validationIncumbent;
        trainingInfo.incumbent.test = testIncumbent;
        trainingInfo.incumbent.validationImprovementDb = validationIncumbentImprovementDb;
        trainingInfo.incumbent.testImprovementDb = testIncumbentImprovementDb;
    end
    if incumbentAvailable && ~candidateSelected
        trainingInfo.attemptedCandidate.coefficients = candidateCfg.coefficients;
        trainingInfo.attemptedCandidate.learning = learningDetails;
        trainingInfo.attemptedCandidate.blockMetrics = blockMetrics;
        trainingInfo.attemptedCandidate.rawBlockMetrics = rawBlockMetrics;
        trainingInfo.attemptedCandidate.deploymentSelection = deploymentSelection;
        trainingInfo.attemptedCandidate.trainingNormalization = normalizationInfo;
        trainingInfo.attemptedCandidate.validation = validationCandidate;
        trainingInfo.attemptedCandidate.test = testCandidate;
        trainingInfo.attemptedCandidate.validationImprovementDb = validationCandidateImprovementDb;
        trainingInfo.attemptedCandidate.testImprovementDb = testCandidateImprovementDb;
        trainingInfo.attemptedCandidate.preparationAccepted = candidatePreparationAccepted;
        trainingInfo.attemptedCandidate.preparationReason = candidatePreparationReason;
        trainingInfo.attemptedCandidate.architectureComparisonMetrics = attemptedComparisonMetrics;
        trainingInfo.attemptedCandidate.architectureComparisonAvailable = attemptedComparisonAvailable;
        trainingInfo.attemptedCandidate.architectureComparisonReason = attemptedComparisonReason;
    end

    coefficients = dpdCfg.coefficients;

    comparisonCacheSaveReason = 'identity-start comparison history was not available for caching';
    if attemptedComparisonAvailable
        try
            [comparisonCacheFile, ~] = dpdSaveArchitectureComparisonCache(dpdCfg, attemptedComparisonSource);
            comparisonCacheSaveReason = 'identity-start comparison history saved';
        catch comparisonCacheSaveError
            comparisonCacheSaveReason = ['identity-start comparison history was not saved: ' comparisonCacheSaveError.message];
            fprintf('DPD comparison cache was not saved: %s\n', comparisonCacheSaveError.message);
        end
    elseif cachedComparisonAvailable
        comparisonCacheSaveReason = 'existing compatible identity-start comparison cache retained';
    end
    trainingInfo.architectureComparisonCacheFile = comparisonCacheFile;
    trainingInfo.architectureComparisonCacheSaveReason = comparisonCacheSaveReason;

    if shouldSaveModel
        fitInfo = modelFitInfo(dpdCfg, paCfg, options, referenceGain, validationBaseline, deployedValidation, deployedValidationImprovementDb, testBaseline, deployedTest, deployedTestImprovementDb, coefficientMask, activeIndices, trainingInfo.createdAt);
        saveVerifiedModel(options.modelFile, coefficients, fitInfo, trainingInfo);
    end

    fprintf('%s learning completed.\n', upper(dpdCfg.learningArchitecture));
    fprintf('Validation fixed-G NMSE: PA %.3f dB, deployed DPD+PA %.3f dB ', validationBaseline.nmseDb, deployedValidation.nmseDb);
    fprintf('(improvement %.3f dB).\n', deployedValidationImprovementDb);
    fprintf('Test fixed-G NMSE:       PA %.3f dB, deployed DPD+PA %.3f dB ', testBaseline.nmseDb, deployedTest.nmseDb);
    fprintf('(improvement %.3f dB).\n', deployedTestImprovementDb);
    if incumbentAvailable && candidateSelected
        fprintf('Incumbent decision: candidate accepted because validation and test are nonworse and validation ACLR1 improved by at least %.3f dB.\n', options.incumbentMinimumValidationAclr1ImprovementDb);
    elseif incumbentAvailable
        fprintf('Incumbent decision: candidate rejected (%s); loaded model retained and no model file was written.\n', incumbentDecision.reason);
    end
    if shouldSaveModel
        fprintf('Saved verified model: %s\n', options.modelFile);
    elseif incumbentAvailable && ~candidateSelected
        fprintf('Loaded incumbent remains active: %s\n', dpdCfg.modelFile);
    else
        fprintf('Model passed verification; saving is disabled.\n');
    end
end

function verifyDeploymentBlockMetrics(blockMetrics, deploymentSelection, dpdCfg)
    requiredFields = {'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db', 'coefficientHistoryNormalized', 'deploymentHistorySelected'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(blockMetrics, requiredFields{fieldIndex})
            error('dpdTrainRun:MissingDeploymentMetric', 'The deployment block metrics must contain %s.', requiredFields{fieldIndex});
        end
    end
    if ~logical(blockMetrics.coefficientHistoryNormalized) || ~logical(blockMetrics.deploymentHistorySelected)
        error('dpdTrainRun:InvalidDeploymentHistory', 'The final block metrics must use the selected normalized deployment history.');
    end
    maximumEvmIncreasePercent = configurationValue(dpdCfg, 'adaptationGateMaximumEvmIncreasePercent', 0);
    maximumAclr1DecreaseDb = configurationValue(dpdCfg, 'adaptationGateMaximumAclr1DecreaseDb', 0);
    maximumAclr2DecreaseDb = configurationValue(dpdCfg, 'adaptationGateMaximumAclr2DecreaseDb', 0.10);
    tolerance = configurationValue(dpdCfg, 'adaptationGateComparisonTolerance', 1e-9);
    if any(diff(blockMetrics.evmRmsPercent) > maximumEvmIncreasePercent + tolerance)
        error('dpdTrainRun:EvmDeploymentRegression', 'RMS EVM increases in the selected deployment history.');
    end
    if any(diff(blockMetrics.worstAclr1Db) < -maximumAclr1DecreaseDb - tolerance)
        error('dpdTrainRun:Aclr1DeploymentRegression', 'ACLR1 decreases in the selected deployment history.');
    end
    if any(diff(blockMetrics.worstAclr2Db) < -maximumAclr2DecreaseDb - tolerance)
        error('dpdTrainRun:Aclr2DeploymentRegression', 'ACLR2 exceeds its tolerance in the selected deployment history.');
    end
    metricTolerance = max(tolerance, 1e-8);
    if abs(blockMetrics.evmRmsPercent(end) - deploymentSelection.finalEvmRmsPercent) > metricTolerance || abs(blockMetrics.worstAclr1Db(end) - deploymentSelection.finalAclr1Db) > metricTolerance || abs(blockMetrics.worstAclr2Db(end) - deploymentSelection.finalAclr2Db) > metricTolerance
        error('dpdTrainRun:DeploymentMetricMismatch', 'The final graph point does not match the selected deployment coefficients.');
    end
end

function dpdCfg = prepareDpdConfiguration(dpdCfg, paCfg, cfrCfg)
    requiredDpdFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount', 'learningArchitecture'};
    for fieldIndex = 1:length(requiredDpdFields)
        if ~isfield(dpdCfg, requiredDpdFields{fieldIndex})
            error('dpdTrainRun:MissingDpdConfigurationField', 'cfg.dpd must contain %s.', requiredDpdFields{fieldIndex});
        end
    end
    requiredPaFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount', 'coefficients', 'modelSampleRate', 'referenceInputRms', 'inputBackoffDb', 'rateFactor'};
    for fieldIndex = 1:length(requiredPaFields)
        if ~isfield(paCfg, requiredPaFields{fieldIndex})
            error('dpdTrainRun:MissingPaConfigurationField', 'cfg.pa must contain %s.', requiredPaFields{fieldIndex});
        end
    end

    learningArchitecture = lower(strtrim( dpdCfg.learningArchitecture));
    if ~strcmp(learningArchitecture, 'indirect') && ~strcmp(learningArchitecture, 'direct')
        error('dpdTrainRun:UnknownLearningArchitecture', ['cfg.dpd.learningArchitecture must be ''indirect'' ' 'or ''direct''.']);
    end
    dpdCfg.learningArchitecture = learningArchitecture;
    [dpdCfg, ~] = dpdRefreshCfrCompatibility(dpdCfg, cfrCfg);
    expectedModelFile = dpdArchitectureArtifactFile(learningArchitecture, dpdCfg.cfrSignature, 'model');
    if ~isfield(dpdCfg, 'modelFile') || isDefaultArchitectureModelFile(dpdCfg.modelFile)
        dpdCfg.modelFile = expectedModelFile;
    end

    if any(~isfinite(paCfg.coefficients(:)))
        error('dpdTrainRun:NonfinitePaCoefficients', 'The PA coefficients must be finite.');
    end
    currentPaModelSignature = dpdPaModelSignature(paCfg);
    modelWasLoaded = false;
    if isfield(dpdCfg, 'modelLoaded')
        loadedState = dpdCfg.modelLoaded;
        modelWasLoaded = (islogical(loadedState) || isnumeric(loadedState)) && isscalar(loadedState) && isreal(loadedState) && isfinite(loadedState) && loadedState ~= 0;
    end
    signatureWasCompatible = isfield(dpdCfg, 'paModelSignature') && ~isempty(dpdCfg.paModelSignature) && isequaln(dpdCfg.paModelSignature, currentPaModelSignature);
    if modelWasLoaded && ~signatureWasCompatible
        dpdCfg.modelLoaded = false;
        dpdCfg.modelStale = true;
        dpdCfg.modelStaleReason = 'loaded DPD model PA signature is missing or incompatible';
    end
    dpdCfg.paModelSignature = currentPaModelSignature;
    if isfield(dpdCfg, 'modelSampleRate') && abs(dpdCfg.modelSampleRate - paCfg.modelSampleRate) > eps(paCfg.modelSampleRate) * 16
        error('dpdTrainRun:SampleRateMismatch', 'DPD and PA model sample rates must be equal.');
    end
    dpdCfg.modelSampleRate = paCfg.modelSampleRate;

    expectedColumns = 1 + (length(dpdCfg.orders) - 1) * length(dpdCfg.envelopeDelays);
    expectedSize = [length(dpdCfg.signalDelays), expectedColumns];
    if ~isfield(dpdCfg, 'coefficients') || ~isequal(size(dpdCfg.coefficients), expectedSize)
        dpdCfg.coefficients = identityCoefficients(dpdCfg);
    end
    [coefficientMask, activeIndices, maskInfo] = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);
    dpdCfg.coefficientMask = coefficientMask;
    dpdCfg.activeCoefficientIndices = activeIndices;
    dpdCfg.numberOfActiveCoefficients = maskInfo.numberOfActiveCoefficients;
    dpdCfg.enabled = true;
end

function options = trainingOptions(dpdCfg)
    options.trainingSeed = configurationValue( dpdCfg, 'trainingSeed', 11);
    options.validationSeed = configurationValue( dpdCfg, 'validationSeed', 12);
    options.testSeed = configurationValue(dpdCfg, 'testSeed', 13);
    options.learningMonitorSeed = configurationValue(dpdCfg, 'learningMonitorSeed', 14);
    options.trainingSampleCount = round(configurationValue( dpdCfg, 'trainingSampleCount', 300000));
    options.validationSampleCount = round(configurationValue( dpdCfg, 'validationSampleCount', 131072));
    options.testSampleCount = round(configurationValue( dpdCfg, 'testSampleCount', 131072));
    options.learningMonitorSampleCount = round(configurationValue(dpdCfg, 'learningMonitorSampleCount', 131072));
    options.blockLength = round(configurationValue( dpdCfg, 'blockLength', configurationValue( dpdCfg, 'adaptationBlockLength', 30000)));
    options.minimumValidationImprovementDb = configurationValue( dpdCfg, 'minimumValidationImprovementDb', 0);
    options.minimumTestImprovementDb = configurationValue( dpdCfg, 'minimumTestImprovementDb', 0);
    options.maximumDpdAveragePowerChangeDb = configurationValue( dpdCfg, 'maximumDpdAveragePowerChangeDb', 0.25);
    options.maximumPaOutputPowerMismatchDb = configurationValue( dpdCfg, 'maximumPaOutputPowerMismatchDb', 0.25);
    options.maximumAclr1DegradationDb = configurationValue( dpdCfg, 'maximumAclr1DegradationDb', 0.1);
    options.maximumAclr2DegradationDb = configurationValue( dpdCfg, 'maximumAclr2DegradationDb', 7.0);
    options.minimumAclr2Db = configurationValue( dpdCfg, 'minimumAclr2Db', 45.0);
    options.incumbentMinimumValidationAclr1ImprovementDb = configurationValue(dpdCfg, 'incumbentMinimumValidationAclr1ImprovementDb', 0.02);
    options.incumbentComparisonTolerance = configurationValue(dpdCfg, 'incumbentComparisonTolerance', 1e-9);
    options.saveModel = logical(configurationValue( dpdCfg, 'saveModel', true));
    options.modelFile = configurationValue(dpdCfg, 'modelFile', dpdArchitectureArtifactFile(dpdCfg.learningArchitecture, dpdCfg.cfrSignature, 'model'));
end

function validateTrainingOptions(options, dpdCfg, paCfg)
    integerFields = {'trainingSeed', 'validationSeed', 'testSeed', 'learningMonitorSeed', 'trainingSampleCount', 'validationSampleCount', 'testSampleCount', 'learningMonitorSampleCount', 'blockLength'};
    for fieldIndex = 1:length(integerFields)
        validateattributes(options.(integerFields{fieldIndex}), {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    end
    if length(unique([options.trainingSeed, options.validationSeed, options.testSeed, options.learningMonitorSeed])) ~= 4
        error('dpdTrainRun:NonindependentSeeds', 'Training, validation, test and learning-monitor seeds must be distinct.');
    end
    combinedMemory = maximumModelDelay(dpdCfg) + maximumModelDelay(paCfg);
    if options.blockLength <= combinedMemory
        error('dpdTrainRun:BlockTooShort', ['blockLength must exceed the combined DPD and PA memory ' 'of %d samples.'], combinedMemory);
    end
    scalarFields = {'minimumValidationImprovementDb', 'minimumTestImprovementDb', 'maximumDpdAveragePowerChangeDb', 'maximumPaOutputPowerMismatchDb', 'maximumAclr1DegradationDb', 'maximumAclr2DegradationDb', 'minimumAclr2Db', 'incumbentMinimumValidationAclr1ImprovementDb', 'incumbentComparisonTolerance'};
    for fieldIndex = 1:length(scalarFields)
        validateattributes(options.(scalarFields{fieldIndex}), {'numeric'}, {'scalar', 'real', 'finite'});
    end
    validateattributes(options.incumbentMinimumValidationAclr1ImprovementDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(options.incumbentComparisonTolerance, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(options.saveModel, {'logical'}, {'scalar'});
end

function metrics = addFairnessChecks(metrics, baselineMetrics, options, coefficientScale, normalizationSucceeded)
    metrics.coefficientScale = coefficientScale;
    metrics.coefficientScaleDb = 20 * log10( max(abs(coefficientScale), realmin));
    metrics.normalizationSucceeded = logical(normalizationSucceeded);
    metrics.dpdAveragePowerChangeDb = powerRatioDb( metrics.dpdOutputAveragePower, metrics.inputAveragePower);
    metrics.paOutputPowerMismatchDb = powerRatioDb( metrics.paOutputAveragePower, baselineMetrics.paOutputAveragePower);
    metrics.powerAccepted = abs(metrics.dpdAveragePowerChangeDb) <= options.maximumDpdAveragePowerChangeDb && abs(metrics.paOutputPowerMismatchDb) <= options.maximumPaOutputPowerMismatchDb;

    metrics.aclr1DegradationDb = baselineMetrics.worstAclr1Db - metrics.worstAclr1Db;
    metrics.aclr2DegradationDb = baselineMetrics.worstAclr2Db - metrics.worstAclr2Db;
    metrics.aclr1Accepted = isfinite(metrics.aclr1DegradationDb) && metrics.aclr1DegradationDb <= options.maximumAclr1DegradationDb;
    metrics.aclr2Accepted = isfinite(metrics.aclr2DegradationDb) && metrics.aclr2DegradationDb <= options.maximumAclr2DegradationDb && metrics.worstAclr2Db >= options.minimumAclr2Db;
    metrics.aclrAccepted = metrics.aclrAvailable && baselineMetrics.aclrAvailable && metrics.aclr1Accepted && metrics.aclr2Accepted;
    metrics.isSafe = metrics.isSafe && metrics.normalizationSucceeded && metrics.powerAccepted && metrics.aclrAccepted;
end

function requireQualifiedResult(baseline, candidate, improvementDb, minimumImprovementDb, stageName)
    if ~baseline.isSafe
        error('dpdTrainRun:InvalidBaseline', 'The PA-only %s baseline is invalid.', stageName);
    end
    if ~candidate.isSafe
        error('dpdTrainRun:UnsafeCandidate', ['The %s DPD candidate failed finite, peak, power, ' 'or ACLR checks.'], stageName);
    end
    if ~isfinite(improvementDb) || improvementDb < minimumImprovementDb
        error('dpdTrainRun:InsufficientImprovement', ['The %s DPD improvement is %.3f dB; at least %.3f dB ' 'is required.'], stageName, improvementDb, minimumImprovementDb);
    end
end

function available = compatibleLoadedIncumbent(dpdCfg)
    available = false;
    if ~isfield(dpdCfg, 'modelLoaded')
        return
    end
    modelLoaded = dpdCfg.modelLoaded;
    if ~(islogical(modelLoaded) || isnumeric(modelLoaded)) || ~isscalar(modelLoaded) || ~isreal(modelLoaded) || ~isfinite(modelLoaded) || modelLoaded == 0
        return
    end
    if isfield(dpdCfg, 'modelStale') && logical(dpdCfg.modelStale)
        return
    end
    if ~isfield(dpdCfg, 'loadedModelArchitecture') || ~ischar(dpdCfg.loadedModelArchitecture) || ~strcmpi(strtrim(dpdCfg.loadedModelArchitecture), dpdCfg.learningArchitecture)
        return
    end
    if ~isfield(dpdCfg, 'coefficients') || ~isnumeric(dpdCfg.coefficients) || any(~isfinite(dpdCfg.coefficients(:)))
        return
    end
    available = true;
end

function matched = comparisonProtocolMatchesConfiguration(metrics, options, monitoringSignal, architecture)
    matched = false;
    requiredMetricFields = {'learningArchitecture', 'trainingSeed', 'trainingSampleCount', 'trainingBlockLength', 'monitorSeed', 'monitorOfdmSymbolCount', 'monitorSampleRate'};
    for fieldIndex = 1:length(requiredMetricFields)
        if ~isfield(metrics, requiredMetricFields{fieldIndex})
            return
        end
    end
    if ~ischar(metrics.learningArchitecture) || ~strcmpi(strtrim(metrics.learningArchitecture), architecture)
        return
    end
    protocolValues = [metrics.trainingSeed, metrics.trainingSampleCount, metrics.trainingBlockLength, metrics.monitorSeed, metrics.monitorOfdmSymbolCount, metrics.monitorSampleRate];
    if ~isnumeric(protocolValues) || any(~isfinite(protocolValues))
        return
    end
    if ~isstruct(monitoringSignal) || ~isfield(monitoringSignal, 'info') || ~isfield(monitoringSignal, 'cfg') || ~isfield(monitoringSignal.info, 'seed') || ~isfield(monitoringSignal.info, 'ofdmSymbolCount') || ~isfield(monitoringSignal.cfg, 'sampleRate')
        return
    end
    sampleRateTolerance = eps(max(abs([metrics.monitorSampleRate, monitoringSignal.cfg.sampleRate]))) * 16;
    matched = metrics.trainingSeed == options.trainingSeed && metrics.trainingSampleCount == options.trainingSampleCount && metrics.trainingBlockLength == options.blockLength && metrics.monitorSeed == monitoringSignal.info.seed && metrics.monitorOfdmSymbolCount == monitoringSignal.info.ofdmSymbolCount && abs(metrics.monitorSampleRate - monitoringSignal.cfg.sampleRate) <= sampleRateTolerance;
end

function value = savedTrainingField(incumbentCfg, fieldName)
    value = struct();
    if ~isfield(incumbentCfg, 'savedTrainingInfo') || ~isstruct(incumbentCfg.savedTrainingInfo) || ~isscalar(incumbentCfg.savedTrainingInfo)
        return
    end
    if ~isfield(incumbentCfg.savedTrainingInfo, fieldName)
        return
    end
    candidateValue = incumbentCfg.savedTrainingInfo.(fieldName);
    if isstruct(candidateValue)
        value = candidateValue;
    end
end

function [candidateSelected, decision] = incumbentReplacementDecision(validationBaseline, validationCandidate, validationCandidateImprovementDb, validationIncumbent, testBaseline, testCandidate, testCandidateImprovementDb, testIncumbent, options, candidatePreparationAccepted, candidatePreparationReason)
    [validationQualified, validationQualificationReason] = qualificationStatus(validationBaseline, validationCandidate, validationCandidateImprovementDb, options.minimumValidationImprovementDb);
    [testQualified, testQualificationReason] = qualificationStatus(testBaseline, testCandidate, testCandidateImprovementDb, options.minimumTestImprovementDb);
    tolerance = options.incumbentComparisonTolerance;
    validationNmseNonWorse = validationCandidate.nmseDb <= validationIncumbent.nmseDb + tolerance;
    validationAclr1NonWorse = validationCandidate.worstAclr1Db >= validationIncumbent.worstAclr1Db - tolerance;
    validationAclr2NonWorse = validationCandidate.worstAclr2Db >= validationIncumbent.worstAclr2Db - tolerance;
    testNmseNonWorse = testCandidate.nmseDb <= testIncumbent.nmseDb + tolerance;
    testAclr1NonWorse = testCandidate.worstAclr1Db >= testIncumbent.worstAclr1Db - tolerance;
    testAclr2NonWorse = testCandidate.worstAclr2Db >= testIncumbent.worstAclr2Db - tolerance;
    validationAclr1ImprovementDb = validationCandidate.worstAclr1Db - validationIncumbent.worstAclr1Db;
    validationAclr1Improved = validationAclr1ImprovementDb + tolerance >= options.incumbentMinimumValidationAclr1ImprovementDb;

    candidateSelected = false;
    if ~candidatePreparationAccepted
        reason = candidatePreparationReason;
    elseif ~validationQualified
        reason = ['candidate_validation_' validationQualificationReason];
    elseif ~testQualified
        reason = ['candidate_test_' testQualificationReason];
    elseif ~validationNmseNonWorse
        reason = 'candidate_validation_nmse_worse';
    elseif ~validationAclr1NonWorse
        reason = 'candidate_validation_aclr1_worse';
    elseif ~validationAclr2NonWorse
        reason = 'candidate_validation_aclr2_worse';
    elseif ~testNmseNonWorse
        reason = 'candidate_test_nmse_worse';
    elseif ~testAclr1NonWorse
        reason = 'candidate_test_aclr1_worse';
    elseif ~testAclr2NonWorse
        reason = 'candidate_test_aclr2_worse';
    elseif ~validationAclr1Improved
        reason = 'candidate_validation_aclr1_gain_below_minimum';
    else
        candidateSelected = true;
        reason = 'candidate_better_than_loaded_incumbent';
    end

    decision.available = true;
    decision.candidateSelected = candidateSelected;
    decision.reason = reason;
    decision.candidatePreparationAccepted = candidatePreparationAccepted;
    decision.candidatePreparationReason = candidatePreparationReason;
    decision.validationQualified = validationQualified;
    decision.validationQualificationReason = validationQualificationReason;
    decision.testQualified = testQualified;
    decision.testQualificationReason = testQualificationReason;
    decision.validationNmseNonWorse = validationNmseNonWorse;
    decision.validationAclr1NonWorse = validationAclr1NonWorse;
    decision.validationAclr2NonWorse = validationAclr2NonWorse;
    decision.testNmseNonWorse = testNmseNonWorse;
    decision.testAclr1NonWorse = testAclr1NonWorse;
    decision.testAclr2NonWorse = testAclr2NonWorse;
    decision.validationAclr1ImprovementDb = validationAclr1ImprovementDb;
    decision.minimumValidationAclr1ImprovementDb = options.incumbentMinimumValidationAclr1ImprovementDb;
    decision.comparisonTolerance = tolerance;
end

function [qualified, reason] = qualificationStatus(baseline, candidate, improvementDb, minimumImprovementDb)
    qualified = false;
    if ~baseline.isSafe
        reason = 'baseline_invalid';
    elseif ~candidate.isSafe
        reason = 'unsafe';
    elseif ~isfinite(improvementDb)
        reason = 'improvement_nonfinite';
    elseif improvementDb < minimumImprovementDb
        reason = 'improvement_insufficient';
    else
        qualified = true;
        reason = 'qualified';
    end
end

function printMetrics(label, metrics)
    fprintf(['%s: fixed-G NMSE %.3f dB, ACLR1/2 ' '%.2f/%.2f dB'], label, metrics.nmseDb, metrics.worstAclr1Db, metrics.worstAclr2Db);
    if isfield(metrics, 'dpdAveragePowerChangeDb')
        fprintf(', DPD power %+0.3f dB, PA power %+0.3f dB', metrics.dpdAveragePowerChangeDb, metrics.paOutputPowerMismatchDb);
    end
    fprintf('.\n');
end

function fitInfo = modelFitInfo(dpdCfg, paCfg, options, referenceGain, validationBaseline, validationCandidate, validationImprovementDb, testBaseline, testCandidate, testImprovementDb, coefficientMask, activeIndices, createdAt)
    fitInfo.learningArchitecture = dpdCfg.learningArchitecture;
    fitInfo.orders = dpdCfg.orders;
    fitInfo.signalDelays = dpdCfg.signalDelays;
    fitInfo.envelopeDelays = dpdCfg.envelopeDelays;
    fitInfo.diagonalCount = dpdCfg.diagonalCount;
    fitInfo.coefficientMask = coefficientMask;
    fitInfo.activeCoefficientIndices = activeIndices;
    fitInfo.numberOfActiveCoefficients = length(activeIndices);
    fitInfo.modelSampleRate = paCfg.modelSampleRate;
    fitInfo.referenceInputRms = paCfg.referenceInputRms * 10^(-paCfg.inputBackoffDb / 20);
    fitInfo.paInputBackoffDb = paCfg.inputBackoffDb;
    fitInfo.maximumOutputMagnitude = configurationValue( dpdCfg, 'maximumOutputMagnitude', Inf);
    fitInfo.referenceLinearGain = referenceGain;
    fitInfo.nmseReference = 'fixed PA-only training gain, common to direct and indirect';
    fitInfo.referenceSignalPoint = 'after interpolation and PA-input scaling, before DPD';
    fitInfo.blockLength = options.blockLength;
    fitInfo.trainingSampleCount = options.trainingSampleCount;
    fitInfo.numberOfTrainingBlocks = ceil( options.trainingSampleCount / options.blockLength);
    fitInfo.validationCascadeNmseDb = validationCandidate.nmseDb;
    fitInfo.validationBaselineNmseDb = validationBaseline.nmseDb;
    fitInfo.validationImprovementDb = validationImprovementDb;
    fitInfo.validationCoefficientScale = validationCandidate.coefficientScale;
    fitInfo.validationDpdAveragePowerChangeDb = validationCandidate.dpdAveragePowerChangeDb;
    fitInfo.validationPaOutputPowerMismatchDb = validationCandidate.paOutputPowerMismatchDb;
    fitInfo.validationAclr1Db = validationCandidate.worstAclr1Db;
    fitInfo.validationAclr2Db = validationCandidate.worstAclr2Db;
    fitInfo.validationBaselineAclr1Db = validationBaseline.worstAclr1Db;
    fitInfo.validationBaselineAclr2Db = validationBaseline.worstAclr2Db;
    fitInfo.testCascadeNmseDb = testCandidate.nmseDb;
    fitInfo.testBaselineNmseDb = testBaseline.nmseDb;
    fitInfo.testImprovementDb = testImprovementDb;
    fitInfo.testDpdAveragePowerChangeDb = testCandidate.dpdAveragePowerChangeDb;
    fitInfo.testPaOutputPowerMismatchDb = testCandidate.paOutputPowerMismatchDb;
    fitInfo.testAclr1Db = testCandidate.worstAclr1Db;
    fitInfo.testAclr2Db = testCandidate.worstAclr2Db;
    fitInfo.testBaselineAclr1Db = testBaseline.worstAclr1Db;
    fitInfo.testBaselineAclr2Db = testBaseline.worstAclr2Db;
    fitInfo.trainingSeed = options.trainingSeed;
    fitInfo.validationSeed = options.validationSeed;
    fitInfo.testSeed = options.testSeed;
    fitInfo.learningMonitorSeed = options.learningMonitorSeed;
    fitInfo.learningMonitorSampleCount = options.learningMonitorSampleCount;
    fitInfo.createdAt = createdAt;
    fitInfo.paModelSignature = dpdPaModelSignature(paCfg);
    fitInfo.cfrSignature = dpdCfg.cfrSignature;
end

function candidateCfg = candidateConfiguration(dpdCfg, coefficients)
    candidateCfg = dpdCfg;
    candidateCfg.enabled = true;
    candidateCfg.coefficients = coefficients;
end

function identityCfg = identityConfiguration(dpdCfg)
    identityCfg = dpdCfg;
    identityCfg.enabled = true;
    identityCfg.coefficients = identityCoefficients(dpdCfg);
end

function coefficients = identityCoefficients(modelCfg)
    numberOfColumns = 1 + (length(modelCfg.orders) - 1) * length(modelCfg.envelopeDelays);
    coefficients = complex(zeros( length(modelCfg.signalDelays), numberOfColumns));
    coefficients(1, 1) = 1;
end

function description = algorithmDescription(learningArchitecture)
    if strcmp(learningArchitecture, 'direct')
        description = [ 'Direct sequential mini-batch learning from the ' 'post-interpolation reference ' 'through a differentiable DPD-PA cascade'];
    else
        description = [ 'Indirect block learning of a postdistorter from normalized ' 'PA output to the signal after DPD'];
    end
end

function matched = isDefaultArchitectureModelFile(modelFile)
    matched = false;
    if ~ischar(modelFile)
        return;
    end
    [~, modelName, modelExtension] = fileparts(modelFile);
    defaultNames = {'dpdModel', 'dpdModelIndirect', 'dpdModelDirect', 'dpdModelIndirectCfr', 'dpdModelDirectCfr'};
    matched = strcmpi(modelExtension, '.mat') && any(strcmpi(modelName, defaultNames));
end

function modelFile = resolveModelFile(configuredFile, paCfg, saveModel)
    if ~ischar(configuredFile) || isempty(strtrim(configuredFile))
        error('dpdTrainRun:InvalidModelFile', 'dpdCfg.modelFile must be a nonempty character vector.');
    end
    [modelDirectory, modelName, modelExtension] = fileparts(configuredFile);
    if ~strcmpi(modelExtension, '.mat')
        error('dpdTrainRun:InvalidModelExtension', 'The DPD model file must use the .mat extension.');
    end
    projectDirectory = fileparts(mfilename('fullpath'));
    if isempty(modelDirectory)
        modelDirectory = projectDirectory;
    elseif isempty(regexp(modelDirectory, '^[A-Za-z]:[\\/]|^\\\\', 'once'))
        modelDirectory = fullfile(projectDirectory, modelDirectory);
    end
    if exist(modelDirectory, 'dir') ~= 7
        error('dpdTrainRun:MissingModelDirectory', 'The DPD model directory does not exist: %s', modelDirectory);
    end
    modelFile = fullfile(modelDirectory, [modelName modelExtension]);
    if saveModel && isfield(paCfg, 'modelFile') && sameWindowsPath(modelFile, paCfg.modelFile)
        error('dpdTrainRun:PaModelCollision', 'The DPD model file must not overwrite the PA model file.');
    end
end

function same = sameWindowsPath(firstPath, secondPath)
    same = strcmpi(strrep(firstPath, '/', '\'), strrep(secondPath, '/', '\'));
end

function saveVerifiedModel(modelFile, coefficients, fitInfo, trainingInfo)
    modelDirectory = fileparts(modelFile);
    temporaryModelFile = [tempname(modelDirectory) '.mat'];
    temporaryCleanup = onCleanup( @() deleteTemporaryFile(temporaryModelFile));
    if ~isnumeric(coefficients) || ~isstruct(fitInfo) || ~isstruct(trainingInfo)
        error('dpdTrainRun:InvalidSavedModelData', 'The verified DPD model data have invalid types.');
    end
    save(temporaryModelFile, 'coefficients', 'fitInfo', 'trainingInfo', '-v7');
    [moveSucceeded, moveMessage] = movefile( temporaryModelFile, modelFile, 'f');
    if ~moveSucceeded
        error('dpdTrainRun:ModelSaveFailed', 'Could not replace the DPD model: %s', moveMessage);
    end
    clear temporaryCleanup;
end

function deleteTemporaryFile(fileName)
    if exist(fileName, 'file') == 2
        delete(fileName);
    end
end

function maximumDelay = maximumModelDelay(modelCfg)
    maximumDelay = max([ modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end

function ratioDb = powerRatioDb(numerator, denominator)
    if ~isfinite(numerator) || ~isfinite(denominator) || numerator < 0 || denominator <= realmin
        ratioDb = Inf;
    else
        ratioDb = 10 * log10(max(numerator, realmin) / denominator);
    end
end

function requireFiniteCoefficients(coefficients)
    if isempty(coefficients) || any(~isfinite(coefficients(:)))
        error('dpdTrainRun:NonfiniteCoefficients', 'The learned DPD coefficients are nonfinite.');
    end
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end

function cfrCfg = systemCfrConfiguration(cfg)
    if isfield(cfg, 'cfr') && isstruct(cfg.cfr) && isscalar(cfg.cfr)
        cfrCfg = cfg.cfr;
    else
        cfrCfg = cfrConfig(false);
    end
end
