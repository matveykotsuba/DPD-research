function [coefficients, info] = dpdIndirectLearn(modelInput, paCfg, dpdCfg, referenceGain, externalMonitorInput, externalStressInput)

    if nargin < 4
        error('dpdIndirectLearn:MissingInput', ['modelInput, paCfg, dpdCfg and referenceGain are ' 'required.']);
    end

    validateModelInput(modelInput);
    modelInput = modelInput(:);
    validateModelConfiguration(paCfg, 'PA');
    validateModelConfiguration(dpdCfg, 'DPD');
    validateReferenceGain(referenceGain);

    blockLength = resolveBlockLength(dpdCfg);
    dpdMemory = maximumModelDelay(dpdCfg);
    paMemory = maximumModelDelay(paCfg);
    combinedMemory = dpdMemory + paMemory;

    if blockLength - combinedMemory < 5
        error('dpdIndirectLearn:BlockTooShort', ['The learning block must contain at least five usable ' 'samples after the combined DPD and PA memory crop.']);
    end

    coefficientMask = gmpDiagonalMask( dpdCfg, dpdCfg.diagonalCount);
    activeIndices = find(coefficientMask(:));
    [currentCoefficients, initialization] = initialCoefficients(dpdCfg);
    currentCoefficients(~coefficientMask) = 0;

    maximumCoefficientNorm = configurationValue( dpdCfg, 'maximumCoefficientNorm', Inf);
    maximumCoefficientUpdateNorm = configurationValue( dpdCfg, 'maximumCoefficientUpdateNorm', Inf);
    maximumOutputMagnitude = configurationValue( dpdCfg, 'maximumOutputMagnitude', Inf);
    lineSearchFactors = configurationValue( dpdCfg, 'indirectLineSearchFactors', [1, 0.5, 0.25, 0.125, 0.0625]);
    minimumMonitorImprovementDb = configurationValue( dpdCfg, 'indirectMinimumMonitorImprovementDb', 0);
    monitorSampleCount = round(configurationValue( dpdCfg, 'indirectMonitorSampleCount', blockLength));
    indirectFitMode = resolveIndirectFitMode(dpdCfg);
    regularizationCenterMode = resolveRegularizationCenterMode(dpdCfg);
    indirectFitForgettingFactor = configurationValue(dpdCfg, 'indirectFitForgettingFactor', 1);
    validateLimit(maximumCoefficientNorm, 'maximumCoefficientNorm', false);
    validateLimit(maximumCoefficientUpdateNorm, 'maximumCoefficientUpdateNorm', true);
    validateLimit(maximumOutputMagnitude, 'maximumOutputMagnitude', false);
    validateLineSearchFactors(lineSearchFactors);
    validateRegularizationConfiguration(dpdCfg);
    validateattributes(minimumMonitorImprovementDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(monitorSampleCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateIndirectFitForgettingFactor(indirectFitForgettingFactor);
    lineSearchFactors = unique(double(lineSearchFactors(:).'), 'stable');

    sampleCount = length(modelInput);
    blockCount = ceil(sampleCount / blockLength);
    coefficientCount = numel(coefficientMask);
    monitorIsIndependent = nargin >= 5 && ~isempty(externalMonitorInput);
    if monitorIsIndependent
        validateModelInput(externalMonitorInput);
        monitorInput = externalMonitorInput(:);
        monitorSampleCount = length(monitorInput);
        monitorFirstSample = NaN;
        monitorLastSample = NaN;
    else
        monitorSampleCount = min(monitorSampleCount, sampleCount);
        monitorFirstSample = sampleCount - monitorSampleCount + 1;
        monitorLastSample = sampleCount;
        monitorInput = modelInput(monitorFirstSample:monitorLastSample);
    end
    if monitorSampleCount - combinedMemory < 5
        error('dpdIndirectLearn:MonitorSegmentTooShort', 'The fixed monitor segment must contain at least five usable samples after the combined DPD and PA memory crop.');
    end
    stressIsIndependent = nargin >= 6 && ~isempty(externalStressInput);
    if stressIsIndependent
        validateModelInput(externalStressInput);
        stressInput = externalStressInput(:);
    else
        stressInput = monitorInput;
    end
    if length(stressInput) - combinedMemory < 5
        error('dpdIndirectLearn:StressSegmentTooShort', 'The stress segment must contain at least five usable samples after the combined DPD and PA memory crop.');
    end
    permittedOutputMagnitude = permittedDpdPeak(dpdCfg, stressInput, combinedMemory);

    blockHistory.blockIndex = (1:blockCount).';
    blockHistory.firstSample = zeros(blockCount, 1);
    blockHistory.lastSample = zeros(blockCount, 1);
    blockHistory.sampleCount = zeros(blockCount, 1);
    blockHistory.usableSampleCount = zeros(blockCount, 1);
    blockHistory.coefficientNormBefore = NaN(blockCount, 1);
    blockHistory.candidateCoefficientNorm = NaN(blockCount, 1);
    blockHistory.coefficientNormAfter = NaN(blockCount, 1);
    blockHistory.coefficientUpdateNorm = NaN(blockCount, 1);
    blockHistory.bestLambda = NaN(blockCount, 1);
    blockHistory.fitValidationNmseDb = NaN(blockCount, 1);
    blockHistory.fitTestNmseDb = NaN(blockCount, 1);
    blockHistory.fitStateUpdateCount = zeros(blockCount, 1);
    blockHistory.fitTrainingSampleCount = zeros(blockCount, 1);
    blockHistory.fitValidationSampleCount = zeros(blockCount, 1);
    blockHistory.fitTestSampleCount = zeros(blockCount, 1);
    blockHistory.forwardCascadeNmseDb = NaN(blockCount, 1);
    blockHistory.dpdOutputRms = NaN(blockCount, 1);
    blockHistory.paOutputRms = NaN(blockCount, 1);
    blockHistory.candidateOutputPeak = NaN(blockCount, 1);
    blockHistory.monitorNmseBeforeDb = NaN(blockCount, 1);
    blockHistory.monitorNmseCandidateDb = NaN(blockCount, 1);
    blockHistory.monitorNmseAfterDb = NaN(blockCount, 1);
    blockHistory.selectedLineSearchFactor = zeros(blockCount, 1);
    blockHistory.updateAccepted = false(blockCount, 1);
    blockHistory.rejectionReason = cell(blockCount, 1);

    coefficientHistory = complex(zeros( coefficientCount, blockCount + 1));
    coefficientHistory(:, 1) = currentCoefficients(:);

    dpdState = [];
    paState = [];
    indirectFitState = [];

    for blockIndex = 1:blockCount
        firstSample = (blockIndex - 1) * blockLength + 1;
        lastSample = min(blockIndex * blockLength, sampleCount);
        inputBlock = modelInput(firstSample:lastSample);

        forwardCfg = dpdCfg;
        forwardCfg.coefficients = currentCoefficients;
        forwardCfg.enabled = true;

        [dpdOutputBlock, dpdState] = gmpCoreBlock( inputBlock, forwardCfg, dpdState);
        [paOutputBlock, paState] = gmpCoreBlock( dpdOutputBlock, paCfg, paState);

        if any(~isfinite(dpdOutputBlock)) || any(~isfinite(paOutputBlock))
            error('dpdIndirectLearn:NonfiniteForwardSignal', 'Block %d produced a nonfinite DPD or PA output.', blockIndex);
        end

        currentNorm = norm(currentCoefficients(activeIndices));
        blockSampleCount = length(inputBlock);
        usableSampleCount = max(0, blockSampleCount - combinedMemory);

        blockHistory.firstSample(blockIndex) = firstSample;
        blockHistory.lastSample(blockIndex) = lastSample;
        blockHistory.sampleCount(blockIndex) = blockSampleCount;
        blockHistory.usableSampleCount(blockIndex) = usableSampleCount;
        blockHistory.coefficientNormBefore(blockIndex) = currentNorm;
        blockHistory.dpdOutputRms(blockIndex) = sqrt( mean(abs(dpdOutputBlock).^2));
        blockHistory.paOutputRms(blockIndex) = sqrt( mean(abs(paOutputBlock).^2));

        if usableSampleCount >= 5
            validRows = combinedMemory + 1:blockSampleCount;
            blockHistory.forwardCascadeNmseDb(blockIndex) = fixedGainNmseDb(inputBlock(validRows), paOutputBlock(validRows), referenceGain);

            fitRows = paMemory + 1:blockSampleCount;
            postdistorterInput = paOutputBlock(fitRows) / referenceGain;
            postdistorterTarget = dpdOutputBlock(fitRows);

            fitSucceeded = true;
            try
                if strcmp(indirectFitMode, 'cumulative')
                    if strcmp(regularizationCenterMode, 'current')
                        [candidateCoefficients, fitInfo, indirectFitState] = dpdBatchFit(postdistorterTarget, postdistorterInput, dpdCfg, indirectFitState, currentCoefficients);
                    else
                        [candidateCoefficients, fitInfo, indirectFitState] = dpdBatchFit(postdistorterTarget, postdistorterInput, dpdCfg, indirectFitState);
                    end
                else
                    [candidateCoefficients, fitInfo] = dpdBatchFit(postdistorterTarget, postdistorterInput, dpdCfg);
                end
            catch fitError
                if ~isRecoverableFitFailure(fitError)
                    rethrow(fitError)
                end
                fitSucceeded = false;
                blockHistory.rejectionReason{blockIndex} = numericalFitFailureReason(fitError);
            end

            if fitSucceeded
                candidateCoefficients(~coefficientMask) = 0;

                candidateNorm = norm( candidateCoefficients(activeIndices));
                updateNorm = norm( candidateCoefficients(activeIndices) - currentCoefficients(activeIndices));

                blockHistory.candidateCoefficientNorm(blockIndex) = candidateNorm;
                blockHistory.coefficientUpdateNorm(blockIndex) = updateNorm;
                blockHistory.bestLambda(blockIndex) = fitInfo.bestLambda;
                blockHistory.fitValidationNmseDb(blockIndex) = fitInfo.bestValidationNmseDb;
                blockHistory.fitTestNmseDb(blockIndex) = fitInfo.testNmseDb;
                if strcmp(indirectFitMode, 'cumulative')
                    blockHistory.fitStateUpdateCount(blockIndex) = fitInfo.fitStateUpdateCount;
                    blockHistory.fitTrainingSampleCount(blockIndex) = fitInfo.accumulatedTrainingSampleCount;
                    blockHistory.fitValidationSampleCount(blockIndex) = fitInfo.accumulatedValidationSampleCount;
                    blockHistory.fitTestSampleCount(blockIndex) = fitInfo.accumulatedTestSampleCount;
                else
                    blockHistory.fitStateUpdateCount(blockIndex) = 1;
                    blockHistory.fitTrainingSampleCount(blockIndex) = length(fitInfo.trainingRows);
                    blockHistory.fitValidationSampleCount(blockIndex) = length(fitInfo.validationRows);
                    blockHistory.fitTestSampleCount(blockIndex) = length(fitInfo.testRows);
                end

                monitorNmseBeforeDb = monitorCascadeNmseDb(currentCoefficients, monitorInput, dpdCfg, paCfg, referenceGain, combinedMemory);
                if ~isfinite(monitorNmseBeforeDb)
                    error('dpdIndirectLearn:InvalidMonitorBaseline', 'The current coefficients produced a nonfinite fixed-monitor NMSE in block %d.', blockIndex);
                end
                [~, candidateReason, candidatePeak, monitorNmseCandidateDb] = validateCandidate(candidateCoefficients, candidateNorm, updateNorm, monitorInput, stressInput, dpdCfg, paCfg, referenceGain, combinedMemory, maximumCoefficientNorm, maximumCoefficientUpdateNorm, permittedOutputMagnitude);
                blockHistory.candidateOutputPeak(blockIndex) = candidatePeak;
                blockHistory.monitorNmseBeforeDb(blockIndex) = monitorNmseBeforeDb;
                blockHistory.monitorNmseCandidateDb(blockIndex) = monitorNmseCandidateDb;

                candidateAccepted = false;
                selectedCoefficients = currentCoefficients;
                selectedMonitorNmseDb = monitorNmseBeforeDb;
                for factorIndex = 1:length(lineSearchFactors)
                    lineSearchFactor = lineSearchFactors(factorIndex);
                    trialCoefficients = currentCoefficients + lineSearchFactor * (candidateCoefficients - currentCoefficients);
                    trialCoefficients(~coefficientMask) = 0;
                    trialNorm = norm(trialCoefficients(activeIndices));
                    trialUpdateNorm = norm(trialCoefficients(activeIndices) - currentCoefficients(activeIndices));
                    [trialSafe, ~, ~, trialMonitorNmseDb] = validateCandidate(trialCoefficients, trialNorm, trialUpdateNorm, monitorInput, stressInput, dpdCfg, paCfg, referenceGain, combinedMemory, maximumCoefficientNorm, maximumCoefficientUpdateNorm, permittedOutputMagnitude);
                    trialImprovementDb = monitorNmseBeforeDb - trialMonitorNmseDb;
                    if trialSafe && isfinite(trialMonitorNmseDb) && trialImprovementDb > minimumMonitorImprovementDb && trialMonitorNmseDb < selectedMonitorNmseDb
                        candidateAccepted = true;
                        selectedCoefficients = trialCoefficients;
                        selectedMonitorNmseDb = trialMonitorNmseDb;
                        blockHistory.selectedLineSearchFactor(blockIndex) = lineSearchFactor;
                    end
                end

                if candidateAccepted
                    currentCoefficients = selectedCoefficients;
                    blockHistory.updateAccepted(blockIndex) = true;
                    rejectionReason = 'accepted';
                elseif strcmp(candidateReason, 'accepted')
                    rejectionReason = 'no line-search factor improved fixed-monitor NMSE';
                else
                    rejectionReason = strjoin({candidateReason, 'no line-search factor improved fixed-monitor NMSE'}, '; ');
                end
                blockHistory.monitorNmseAfterDb(blockIndex) = selectedMonitorNmseDb;
                blockHistory.rejectionReason{blockIndex} = rejectionReason;
            end
        else
            blockHistory.rejectionReason{blockIndex} = 'insufficient usable samples';
        end

        blockHistory.coefficientNormAfter(blockIndex) = norm( currentCoefficients(activeIndices));
        coefficientHistory(:, blockIndex + 1) = currentCoefficients(:);
    end

    coefficients = currentCoefficients;

    info.method = 'indirect learning architecture';
    info.trainingPair = 'PA output / referenceGain -> current DPD output';
    info.inputReferencePlane = 'after interpolation and PA-input scaling';
    if strcmp(indirectFitMode, 'cumulative')
        info.estimator = 'cumulative Tikhonov-regularized batch LS';
    else
        info.estimator = 'blockwise Tikhonov-regularized batch LS';
    end
    info.indirectFitMode = indirectFitMode;
    info.indirectFitForgettingFactor = indirectFitForgettingFactor;
    info.regularizationCenterMode = regularizationCenterMode;
    info.initialization = initialization;
    info.identityInitialization = strcmp(initialization, 'identity');
    info.coefficientTransferDelayBlocks = 1;
    info.referenceGain = referenceGain;
    info.monitorSampleCount = monitorSampleCount;
    info.monitorFirstSample = monitorFirstSample;
    info.monitorLastSample = monitorLastSample;
    info.monitorIsIndependent = monitorIsIndependent;
    info.stressSampleCount = length(stressInput);
    info.stressIsIndependent = stressIsIndependent;
    info.stressInputPeakMagnitude = max(abs(stressInput));
    info.permittedDpdOutputPeakMagnitude = permittedOutputMagnitude;
    info.lineSearchFactors = lineSearchFactors;
    info.minimumMonitorImprovementDb = minimumMonitorImprovementDb;
    info.blockLength = blockLength;
    info.sampleCount = sampleCount;
    info.blockCount = blockCount;
    info.dpdMemory = dpdMemory;
    info.paMemory = paMemory;
    info.combinedMemory = combinedMemory;
    info.coefficientMask = coefficientMask;
    info.activeCoefficientIndices = activeIndices;
    info.numberOfActiveCoefficients = length(activeIndices);
    info.blockHistory = blockHistory;
    info.coefficientHistory = coefficientHistory;
    info.acceptedUpdateCount = sum(blockHistory.updateAccepted);
    info.rejectedUpdateCount = blockCount - info.acceptedUpdateCount;
    info.finalCoefficientNorm = norm(coefficients(activeIndices));
    if isempty(indirectFitState)
        info.fitStateUpdateCount = 0;
        info.fitTrainingSampleCount = 0;
        info.fitValidationSampleCount = 0;
        info.fitTestSampleCount = 0;
    else
        info.fitStateUpdateCount = indirectFitState.updateCount;
        info.fitTrainingSampleCount = indirectFitState.training.N;
        info.fitValidationSampleCount = indirectFitState.validation.N;
        info.fitTestSampleCount = indirectFitState.test.N;
    end
end

function validateModelInput(modelInput)
    if ~isnumeric(modelInput) || ~isvector(modelInput) || isempty(modelInput) || any(~isfinite(modelInput(:)))
        error('dpdIndirectLearn:InvalidModelInput', 'modelInput must be a nonempty finite numeric vector.');
    end
end

function validateModelConfiguration(modelCfg, description)
    requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'coefficients'};
    for fieldIndex = 1:length(requiredFields)
        fieldName = requiredFields{fieldIndex};
        if ~isfield(modelCfg, fieldName)
            error('dpdIndirectLearn:MissingConfigurationField', '%s configuration must contain %s.', description, fieldName);
        end
    end

    if strcmp(description, 'DPD') && ~isfield(modelCfg, 'diagonalCount')
        error('dpdIndirectLearn:MissingDiagonalCount', 'DPD configuration must contain diagonalCount.');
    end

    expectedColumns = 1 + (length(modelCfg.orders) - 1) * length(modelCfg.envelopeDelays);
    expectedSize = [length(modelCfg.signalDelays), expectedColumns];
    if ~isequal(size(modelCfg.coefficients), expectedSize)
        error('dpdIndirectLearn:InvalidCoefficientSize', ['%s coefficients have size %d-by-%d; expected ' '%d-by-%d.'], description, size(modelCfg.coefficients, 1), size(modelCfg.coefficients, 2), expectedSize(1), expectedSize(2));
    end
    if any(~isfinite(modelCfg.coefficients(:)))
        error('dpdIndirectLearn:NonfiniteModelCoefficients', '%s coefficients must be finite.', description);
    end
end

function validateReferenceGain(referenceGain)
    if ~isnumeric(referenceGain) || ~isscalar(referenceGain) || ~isfinite(referenceGain) || abs(referenceGain) <= realmin
        error('dpdIndirectLearn:InvalidReferenceGain', 'referenceGain must be a finite nonzero numeric scalar.');
    end
end

function blockLength = resolveBlockLength(dpdCfg)
    if isfield(dpdCfg, 'blockLength')
        blockLength = dpdCfg.blockLength;
    elseif isfield(dpdCfg, 'learningBlockLength')
        blockLength = dpdCfg.learningBlockLength;
    elseif isfield(dpdCfg, 'adaptationBlockLength')
        blockLength = dpdCfg.adaptationBlockLength;
    else
        blockLength = 30000;
    end

    validateattributes(blockLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    blockLength = double(blockLength);
end

function maximumDelay = maximumModelDelay(modelCfg)
    maximumDelay = max([ modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end

function coefficients = identityCoefficients(modelCfg)
    numberOfColumns = 1 + (length(modelCfg.orders) - 1) * length(modelCfg.envelopeDelays);
    coefficients = complex(zeros( length(modelCfg.signalDelays), numberOfColumns));
    coefficients(1, 1) = 1;
end

function [coefficients, initialization] = initialCoefficients(dpdCfg)
    initialization = configurationValue(dpdCfg, 'learningInitialization', 'identity');
    if ~ischar(initialization)
        error('dpdIndirectLearn:InvalidInitialization', 'learningInitialization must be a character vector.');
    end
    initialization = lower(strtrim(initialization));
    if strcmp(initialization, 'configured')
        coefficients = dpdCfg.coefficients;
    elseif strcmp(initialization, 'identity')
        coefficients = identityCoefficients(dpdCfg);
    else
        error('dpdIndirectLearn:UnknownInitialization', ['learningInitialization must be identity or ' 'configured.']);
    end
end

function validateLimit(value, fieldName, allowZero)
    valid = isnumeric(value) && isscalar(value) && isreal(value) && ~isnan(value);
    if allowZero
        valid = valid && value >= 0;
    else
        valid = valid && value > 0;
    end
    if ~valid
        error('dpdIndirectLearn:InvalidLimit', '%s must be a positive real scalar or Inf.', fieldName);
    end
end

function [accepted, reason, candidatePeak, monitorNmseDb] = validateCandidate(candidateCoefficients, candidateNorm, updateNorm, monitorInput, stressInput, dpdCfg, paCfg, referenceGain, combinedMemory, maximumCoefficientNorm, maximumCoefficientUpdateNorm, permittedOutputMagnitude)
    reasons = cell(0, 1);

    if any(~isfinite(candidateCoefficients(:))) || ~isfinite(candidateNorm)
        reasons{end + 1, 1} = 'nonfinite coefficients';
    elseif candidateNorm > maximumCoefficientNorm
        reasons{end + 1, 1} = 'coefficient norm limit exceeded';
    end

    if ~isfinite(updateNorm) || updateNorm > maximumCoefficientUpdateNorm
        reasons{end + 1, 1} = 'coefficient update norm limit exceeded';
    end

    candidatePeak = Inf;
    monitorNmseDb = NaN;
    if isempty(reasons)
        candidateCfg = dpdCfg;
        candidateCfg.coefficients = candidateCoefficients;
        candidateCfg.enabled = true;
        stressOutput = gmpCore(stressInput, candidateCfg);
        stressValidRows = combinedMemory + 1:length(stressInput);
        if any(~isfinite(stressOutput(stressValidRows)))
            reasons{end + 1, 1} = 'candidate stress output is nonfinite';
        else
            candidatePeak = max(abs(stressOutput(stressValidRows)));
            if candidatePeak > permittedOutputMagnitude
                reasons{end + 1, 1} = 'candidate stress output magnitude limit exceeded';
            else
                candidateOutput = gmpCore(monitorInput, candidateCfg);
                validRows = combinedMemory + 1:length(monitorInput);
                if any(~isfinite(candidateOutput(validRows)))
                    reasons{end + 1, 1} = 'candidate monitor output is nonfinite';
                else
                    monitorPaOutput = gmpCore(candidateOutput, paCfg);
                    monitorNmseDb = fixedGainNmseDb(monitorInput(validRows), monitorPaOutput(validRows), referenceGain);
                    if ~isfinite(monitorNmseDb)
                        reasons{end + 1, 1} = 'candidate fixed-monitor NMSE is nonfinite';
                    end
                end
            end
        end
    end

    accepted = isempty(reasons);
    if accepted
        reason = 'accepted';
    else
        reason = strjoin(reasons, '; ');
    end
end

function permittedPeak = permittedDpdPeak(dpdCfg, stressInput, combinedMemory)
    validRows = combinedMemory + 1:length(stressInput);
    inputPeak = max(abs(stressInput(validRows)));
    maximumPeakGainDb = configurationValue(dpdCfg, 'maximumDpdPeakGainDb', 12);
    maximumOutputMagnitude = configurationValue(dpdCfg, 'maximumOutputMagnitude', Inf);
    validateattributes(maximumPeakGainDb, {'numeric'}, {'scalar', 'real'});
    validateLimit(maximumOutputMagnitude, 'maximumOutputMagnitude', false);
    relativePeakLimit = inputPeak * 10^(maximumPeakGainDb / 20);
    permittedPeak = min(relativePeakLimit, maximumOutputMagnitude);
end

function nmseDb = monitorCascadeNmseDb(coefficients, monitorInput, dpdCfg, paCfg, referenceGain, combinedMemory)
    monitorCfg = dpdCfg;
    monitorCfg.coefficients = coefficients;
    monitorCfg.enabled = true;
    monitorDpdOutput = gmpCore(monitorInput, monitorCfg);
    if any(~isfinite(monitorDpdOutput(:)))
        nmseDb = Inf;
        return
    end
    monitorPaOutput = gmpCore(monitorDpdOutput, paCfg);
    validRows = combinedMemory + 1:length(monitorInput);
    nmseDb = fixedGainNmseDb(monitorInput(validRows), monitorPaOutput(validRows), referenceGain);
end

function validateLineSearchFactors(lineSearchFactors)
    valid = isnumeric(lineSearchFactors) && isvector(lineSearchFactors) && ~isempty(lineSearchFactors) && isreal(lineSearchFactors) && all(isfinite(lineSearchFactors(:))) && all(lineSearchFactors(:) > 0) && all(lineSearchFactors(:) <= 1);
    if ~valid
        error('dpdIndirectLearn:InvalidLineSearchFactors', 'indirectLineSearchFactors must be a finite numeric vector with values in the interval (0, 1].');
    end
end

function validateRegularizationConfiguration(dpdCfg)
    if isfield(dpdCfg, 'regularizationLambdaValues')
        validateattributes(dpdCfg.regularizationLambdaValues, {'numeric'}, {'vector', 'real', 'finite', 'positive', 'nonempty'});
    end
end

function indirectFitMode = resolveIndirectFitMode(dpdCfg)
    indirectFitMode = configurationValue(dpdCfg, 'indirectFitMode', 'block');
    if ~ischar(indirectFitMode)
        error('dpdIndirectLearn:InvalidFitMode', 'indirectFitMode must be a character vector.');
    end
    indirectFitMode = lower(strtrim(indirectFitMode));
    if ~strcmp(indirectFitMode, 'block') && ~strcmp(indirectFitMode, 'cumulative')
        error('dpdIndirectLearn:UnknownFitMode', 'indirectFitMode must be block or cumulative.');
    end
end

function regularizationCenterMode = resolveRegularizationCenterMode(dpdCfg)
    regularizationCenterMode = configurationValue(dpdCfg, 'indirectRegularizationCenter', 'zero');
    if ~ischar(regularizationCenterMode)
        error('dpdIndirectLearn:InvalidRegularizationCenter', 'indirectRegularizationCenter must be a character vector.');
    end
    regularizationCenterMode = lower(strtrim(regularizationCenterMode));
    if ~strcmp(regularizationCenterMode, 'zero') && ~strcmp(regularizationCenterMode, 'current')
        error('dpdIndirectLearn:UnknownRegularizationCenter', 'indirectRegularizationCenter must be zero or current.');
    end
end

function validateIndirectFitForgettingFactor(forgettingFactor)
    valid = isnumeric(forgettingFactor) && isscalar(forgettingFactor) && isreal(forgettingFactor) && isfinite(forgettingFactor) && forgettingFactor == 1;
    if ~valid
        error('dpdIndirectLearn:InvalidFitForgettingFactor', 'indirectFitForgettingFactor must equal 1.');
    end
end

function nmseDb = fixedGainNmseDb(reference, output, referenceGain)
    target = referenceGain * reference(:);
    output = output(:);
    targetEnergy = sum(abs(target).^2);
    errorEnergy = sum(abs(output - target).^2);
    if targetEnergy <= realmin || ~isfinite(targetEnergy) || ~isfinite(errorEnergy)
        nmseDb = NaN;
    else
        nmseDb = 10 * log10(max( errorEnergy / targetEnergy, realmin));
    end
end

function recoverable = isRecoverableFitFailure(fitError)
    identifier = lower(fitError.identifier);
    message = lower(fitError.message);
    eigFailure = contains(identifier, 'eig') || contains(message, 'input to eig') || contains(message, 'eigenvalue');
    invalidNumericalData = contains(message, 'nan') || contains(message, 'inf') || contains(message, 'finite') || contains(message, 'lapack');
    recoverable = eigFailure && invalidNumericalData;
end

function reason = numericalFitFailureReason(fitError)
    if isempty(fitError.identifier)
        detail = regexprep(fitError.message, '\s+', ' ');
    else
        detail = fitError.identifier;
    end
    reason = ['numerical batch fit failure: ' detail];
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
