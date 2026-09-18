function [coefficients, learningInfo] = dpdDirectLearn(modelInput, paCfg, dpdCfg, referenceGain, externalMonitorInput, externalStressInput, externalAclrMonitorInput)

    validateInputs(modelInput, paCfg, dpdCfg, referenceGain);
    modelInput = modelInput(:);

    blockLength = round(configurationValue(dpdCfg, 'blockLength', configurationValue(dpdCfg, 'adaptationBlockLength', 30000)));
    maximumEpochs = round(configurationValue( dpdCfg, 'directMaximumEpochs', 1));
    jacobianStride = round(configurationValue( dpdCfg, 'directJacobianStride', 4));
    lambdaValues = configurationValue( dpdCfg, 'directLambdaValues', logspace(-8, 2, 21));
    lineSearchFactors = configurationValue( dpdCfg, 'directLineSearchFactors', [1 0.5 0.25 0.125 0.0625]);
    minimumEpochImprovementDb = configurationValue( dpdCfg, 'directMinimumEpochImprovementDb', 0.01);
    monitorSampleCount = round(configurationValue( dpdCfg, 'directMonitorSampleCount', blockLength));
    minimumMonitorImprovementDb = configurationValue( dpdCfg, 'directMinimumMonitorImprovementDb', 0);
    selectionObjective = resolveDirectSelectionObjective(dpdCfg);
    selectionOptions.objective = selectionObjective;
    selectionOptions.maximumNmseDegradationDb = configurationValue(dpdCfg, 'directMaximumMonitorNmseDegradationDb', 0.25);
    selectionOptions.minimumAclr1ImprovementDb = configurationValue(dpdCfg, 'directMinimumMonitorAclr1ImprovementDb', 0.02);
    selectionOptions.maximumAclr2DegradationDb = configurationValue(dpdCfg, 'directMaximumMonitorAclr2DegradationDb', 1);
    selectionOptions.spectrumSegmentLength = round(configurationValue(dpdCfg, 'directAclrSpectrumSegmentLength', 8192));
    selectionOptions.candidateCount = round(configurationValue(dpdCfg, 'directAclrCandidateCount', 5));
    selectionOptions.channelBandwidth = configurationValue(dpdCfg, 'directAclrChannelBandwidth', NaN);
    selectionOptions.measurementBandwidth = configurationValue(dpdCfg, 'directAclrMeasurementBandwidth', NaN);
    selectionOptions.useSpectralObjective = configurationValue(dpdCfg, 'directTrainingUsesSpectralObjective', false);
    selectionOptions.spectralAclr1Weight = configurationValue(dpdCfg, 'directSpectralAclr1Weight', 0);
    selectionOptions.spectralAclr2Weight = configurationValue(dpdCfg, 'directSpectralAclr2Weight', 0);

    validateattributes(blockLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(maximumEpochs, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(jacobianStride, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(lambdaValues, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
    validateattributes(lineSearchFactors, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
    validateattributes(minimumEpochImprovementDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(monitorSampleCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(minimumMonitorImprovementDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(selectionOptions.maximumNmseDegradationDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(selectionOptions.minimumAclr1ImprovementDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(selectionOptions.maximumAclr2DegradationDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(selectionOptions.spectrumSegmentLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(selectionOptions.candidateCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    if ~(islogical(selectionOptions.useSpectralObjective) || isnumeric(selectionOptions.useSpectralObjective)) || ~isscalar(selectionOptions.useSpectralObjective) || ~isfinite(selectionOptions.useSpectralObjective)
        error('dpdDirectLearn:InvalidSpectralObjectiveFlag', 'directTrainingUsesSpectralObjective must be a finite logical scalar.');
    end
    selectionOptions.useSpectralObjective = logical(selectionOptions.useSpectralObjective) && strcmp(selectionObjective, 'aclr');
    validateattributes(selectionOptions.spectralAclr1Weight, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(selectionOptions.spectralAclr2Weight, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    if strcmp(selectionObjective, 'aclr')
        validateattributes(selectionOptions.channelBandwidth, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
        validateattributes(selectionOptions.measurementBandwidth, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    end

    lambdaValues = lambdaValues(:).';
    lineSearchFactors = lineSearchFactors(:).';
    coefficientMask = gmpDiagonalMask( dpdCfg, dpdCfg.diagonalCount);
    activeIndices = find(coefficientMask(:));
    coefficientCount = length(activeIndices);

    workingCfg = dpdCfg;
    workingCfg.enabled = true;
    [workingCfg.coefficients, initialization] = initialCoefficients(dpdCfg);
    workingCfg.coefficients(~coefficientMask) = 0;

    dpdMemory = maximumModelDelay(dpdCfg);
    paMemory = maximumModelDelay(paCfg);
    overlapLength = dpdMemory + paMemory;
    if blockLength <= overlapLength
        error('dpdDirectLearn:BlockTooShort', ['blockLength must exceed the combined DPD and PA memory ' 'of %d samples.'], overlapLength);
    end

    monitorIsIndependent = nargin >= 5 && ~isempty(externalMonitorInput);
    if monitorIsIndependent
        if ~isnumeric(externalMonitorInput) || ~isvector(externalMonitorInput) || isempty(externalMonitorInput) || any(~isfinite(externalMonitorInput(:)))
            error('dpdDirectLearn:InvalidMonitorInput', 'externalMonitorInput must be a nonempty finite numeric vector.');
        end
        monitorInput = externalMonitorInput(:);
        monitorSampleCount = length(monitorInput);
        monitorFirstSample = NaN;
        monitorLastSample = NaN;
    else
        monitorSampleCount = min(monitorSampleCount, length(modelInput));
        monitorFirstSample = length(modelInput) - monitorSampleCount + 1;
        monitorLastSample = length(modelInput);
        monitorInput = modelInput(monitorFirstSample:monitorLastSample);
    end
    if monitorSampleCount <= overlapLength
        error('dpdDirectLearn:MonitorSegmentTooShort', 'The fixed monitor segment must exceed the combined DPD and PA memory.');
    end
    stressIsIndependent = nargin >= 6 && ~isempty(externalStressInput);
    if stressIsIndependent
        if ~isnumeric(externalStressInput) || ~isvector(externalStressInput) || isempty(externalStressInput) || any(~isfinite(externalStressInput(:)))
            error('dpdDirectLearn:InvalidStressInput', 'externalStressInput must be a nonempty finite numeric vector.');
        end
        stressInput = externalStressInput(:);
    else
        stressInput = monitorInput;
    end
    if length(stressInput) <= overlapLength
        error('dpdDirectLearn:StressSegmentTooShort', 'The stress segment must exceed the combined DPD and PA memory.');
    end
    aclrMonitorIsIndependent = nargin >= 7 && ~isempty(externalAclrMonitorInput);
    if aclrMonitorIsIndependent
        if ~isnumeric(externalAclrMonitorInput) || ~isvector(externalAclrMonitorInput) || isempty(externalAclrMonitorInput) || any(~isfinite(externalAclrMonitorInput(:)))
            error('dpdDirectLearn:InvalidAclrMonitorInput', 'externalAclrMonitorInput must be a nonempty finite numeric vector.');
        end
        aclrMonitorInput = externalAclrMonitorInput(:);
        aclrMonitorSampleCount = round(configurationValue(dpdCfg, 'directAclrMonitorSampleCount', 262144));
        validateattributes(aclrMonitorSampleCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
        aclrMonitorSampleCount = min(aclrMonitorSampleCount, length(aclrMonitorInput));
        aclrMonitorFirstSample = floor((length(aclrMonitorInput) - aclrMonitorSampleCount) / 2) + 1;
        aclrMonitorInput = aclrMonitorInput(aclrMonitorFirstSample:aclrMonitorFirstSample + aclrMonitorSampleCount - 1);
    else
        aclrMonitorInput = monitorInput;
        aclrMonitorSampleCount = length(aclrMonitorInput);
    end
    if aclrMonitorSampleCount <= overlapLength
        error('dpdDirectLearn:AclrMonitorSegmentTooShort', 'The ACLR monitor segment must exceed the combined DPD and PA memory.');
    end
    permittedOutputMagnitude = permittedDpdPeak(dpdCfg, stressInput, overlapLength);

    blocksPerEpoch = ceil(length(modelInput) / blockLength);
    maximumHistoryLength = blocksPerEpoch * maximumEpochs;
    history = initializeHistory(maximumHistoryLength);
    historyIndex = 0;
    coefficientHistory = complex(zeros(numel(workingCfg.coefficients), maximumHistoryLength + 1));
    coefficientHistory(:, 1) = workingCfg.coefficients(:);

    initialTrainingNmseDb = cascadeNmseDb( modelInput, workingCfg, paCfg, referenceGain, overlapLength);
    [initialMonitorNmseDb, ~] = blockCascadeNmseDb( monitorInput, workingCfg, paCfg, referenceGain, overlapLength);
    initialMonitorPeak = dpdPeak(stressInput, workingCfg, overlapLength);
    currentMonitorNmseDb = initialMonitorNmseDb;
    currentMonitorPeak = initialMonitorPeak;
    if strcmp(selectionObjective, 'aclr')
        [currentMonitorAclr1Db, currentMonitorAclr2Db] = cascadeAclrDb(aclrMonitorInput, workingCfg, paCfg, overlapLength, selectionOptions);
    else
        currentMonitorAclr1Db = NaN;
        currentMonitorAclr2Db = NaN;
    end
    previousEpochNmseDb = initialTrainingNmseDb;
    completedEpochCount = 0;

    for epochIndex = 1:maximumEpochs
        epochAcceptedUpdateCount = 0;

        for firstSample = 1:blockLength:length(modelInput)
            lastSample = min( firstSample + blockLength - 1, length(modelInput));
            extensionFirstSample = max(1, firstSample - overlapLength);
            extendedInput = modelInput( extensionFirstSample:lastSample);

            [normalMatrix, gradient, usedRowCount, sampledNmseBeforeDb, blockNmseBeforeDb, blockPeakBefore] = directNormalEquations( extendedInput, workingCfg, paCfg, referenceGain, activeIndices, jacobianStride, selectionOptions);

            if ~isfinite(blockNmseBeforeDb) || ~isfinite(blockPeakBefore)
                [blockNmseBeforeDb, blockPeakBefore] = blockCascadeNmseDb( extendedInput, workingCfg, paCfg, referenceGain, overlapLength);
            end

            monitorNmseBeforeDb = currentMonitorNmseDb;
            if strcmp(selectionObjective, 'aclr')
                monitorAclr1BeforeDb = currentMonitorAclr1Db;
                monitorAclr2BeforeDb = currentMonitorAclr2Db;
            else
                monitorAclr1BeforeDb = NaN;
                monitorAclr2BeforeDb = NaN;
            end
            monitorPeakBefore = currentMonitorPeak;
            normalEquationsFinite = all(isfinite(normalMatrix(:))) && all(isfinite(gradient(:))) && isfinite(sampledNmseBeforeDb);
            monitorCandidateNmseDb = NaN;
            monitorNmseAfterDb = monitorNmseBeforeDb;
            monitorPeakAfter = monitorPeakBefore;
            selectedLambda = NaN;
            selectedLineSearchFactor = 0;
            accepted = false;
            candidateCfg = workingCfg;
            monitorAclr1CandidateDb = NaN;
            monitorAclr2CandidateDb = NaN;
            monitorAclr1AfterDb = monitorAclr1BeforeDb;
            monitorAclr2AfterDb = monitorAclr2BeforeDb;

            if normalEquationsFinite && isfinite(monitorNmseBeforeDb) && isfinite(monitorPeakBefore)
                [candidateCfg, accepted, selectedLambda, selectedLineSearchFactor, monitorCandidateNmseDb, monitorNmseAfterDb, monitorPeakAfter, monitorAclr1CandidateDb, monitorAclr2CandidateDb, monitorAclr1AfterDb, monitorAclr2AfterDb] = selectDirectUpdate( monitorInput, aclrMonitorInput, stressInput, workingCfg, paCfg, referenceGain, activeIndices, normalMatrix, gradient, usedRowCount, lambdaValues, lineSearchFactors, overlapLength, monitorNmseBeforeDb, monitorPeakBefore, minimumMonitorImprovementDb, permittedOutputMagnitude, monitorAclr1BeforeDb, monitorAclr2BeforeDb, selectionOptions);
            end

            if accepted
                workingCfg = candidateCfg;
                epochAcceptedUpdateCount = epochAcceptedUpdateCount + 1;
                [blockNmseAfterDb, blockPeakAfter] = blockCascadeNmseDb( extendedInput, workingCfg, paCfg, referenceGain, overlapLength);
            else
                blockNmseAfterDb = blockNmseBeforeDb;
                blockPeakAfter = blockPeakBefore;
            end
            currentMonitorNmseDb = monitorNmseAfterDb;
            currentMonitorPeak = monitorPeakAfter;
            currentMonitorAclr1Db = monitorAclr1AfterDb;
            currentMonitorAclr2Db = monitorAclr2AfterDb;

            historyIndex = historyIndex + 1;
            history.epochIndex(historyIndex) = epochIndex;
            history.blockIndex(historyIndex) = ceil(firstSample / blockLength);
            history.firstSample(historyIndex) = firstSample;
            history.lastSample(historyIndex) = lastSample;
            history.sampleCount(historyIndex) = lastSample - firstSample + 1;
            history.jacobianRowCount(historyIndex) = usedRowCount;
            history.sampledNmseBeforeDb(historyIndex) = sampledNmseBeforeDb;
            history.nmseBeforeDb(historyIndex) = blockNmseBeforeDb;
            history.nmseAfterDb(historyIndex) = blockNmseAfterDb;
            history.normalEquationsFinite(historyIndex) = normalEquationsFinite;
            history.monitorNmseBeforeDb(historyIndex) = monitorNmseBeforeDb;
            history.monitorNmseCandidateDb(historyIndex) = monitorCandidateNmseDb;
            history.monitorNmseAfterDb(historyIndex) = monitorNmseAfterDb;
            history.monitorAclr1BeforeDb(historyIndex) = monitorAclr1BeforeDb;
            history.monitorAclr1CandidateDb(historyIndex) = monitorAclr1CandidateDb;
            history.monitorAclr1AfterDb(historyIndex) = monitorAclr1AfterDb;
            history.monitorAclr2BeforeDb(historyIndex) = monitorAclr2BeforeDb;
            history.monitorAclr2CandidateDb(historyIndex) = monitorAclr2CandidateDb;
            history.monitorAclr2AfterDb(historyIndex) = monitorAclr2AfterDb;
            history.monitorDpdPeakBefore(historyIndex) = monitorPeakBefore;
            history.monitorDpdPeakAfter(historyIndex) = monitorPeakAfter;
            history.accepted(historyIndex) = accepted;
            history.lambda(historyIndex) = selectedLambda;
            history.lineSearchFactor(historyIndex) = selectedLineSearchFactor;
            history.coefficientNorm(historyIndex) = norm( workingCfg.coefficients(activeIndices));
            history.dpdPeakBefore(historyIndex) = blockPeakBefore;
            history.dpdPeakAfter(historyIndex) = blockPeakAfter;
            coefficientHistory(:, historyIndex + 1) = workingCfg.coefficients(:);
            fprintf('Direct block %d/%d: accepted=%d, monitor NMSE %.6f to %.6f dB, stress peak %.6f.\n', historyIndex, maximumHistoryLength, accepted, monitorNmseBeforeDb, monitorNmseAfterDb, monitorPeakAfter);
        end

        clear normalMatrix gradient extendedInput candidateCfg
        completedEpochCount = epochIndex;
        epochNmseDb = cascadeNmseDb( modelInput, workingCfg, paCfg, referenceGain, overlapLength);
        epochImprovementDb = previousEpochNmseDb - epochNmseDb;

        if epochAcceptedUpdateCount == 0 || epochImprovementDb < minimumEpochImprovementDb
            break;
        end
        previousEpochNmseDb = epochNmseDb;
    end

    coefficients = workingCfg.coefficients;
    requireFiniteCoefficients(coefficients);
    history = trimHistory(history, historyIndex);
    coefficientHistory = coefficientHistory(:, 1:historyIndex + 1);

    finalTrainingNmseDb = cascadeNmseDb( modelInput, workingCfg, paCfg, referenceGain, overlapLength);
    learningInfo.learningArchitecture = 'direct';
    learningInfo.initialization = initialization;
    learningInfo.identityInitialization = strcmp(initialization, 'identity');
    learningInfo.referenceSignalPoint = 'after interpolation and PA-input scaling, before DPD';
    learningInfo.objective = 'PA(DPD(x)) - referenceGain*x';
    learningInfo.referenceGain = referenceGain;
    learningInfo.blockLength = blockLength;
    learningInfo.trainingSampleCount = length(modelInput);
    learningInfo.numberOfBlocks = blocksPerEpoch;
    learningInfo.maximumEpochs = maximumEpochs;
    learningInfo.completedEpochs = completedEpochCount;
    learningInfo.totalProcessedBlocks = historyIndex;
    learningInfo.jacobianStride = jacobianStride;
    learningInfo.monitorSampleCount = monitorSampleCount;
    learningInfo.monitorFirstSample = monitorFirstSample;
    learningInfo.monitorLastSample = monitorLastSample;
    learningInfo.monitorIsIndependent = monitorIsIndependent;
    learningInfo.stressSampleCount = length(stressInput);
    learningInfo.stressIsIndependent = stressIsIndependent;
    learningInfo.stressInputPeakMagnitude = max(abs(stressInput));
    learningInfo.aclrMonitorSampleCount = aclrMonitorSampleCount;
    learningInfo.aclrMonitorIsIndependent = aclrMonitorIsIndependent;
    learningInfo.permittedDpdOutputPeakMagnitude = permittedOutputMagnitude;
    learningInfo.minimumMonitorImprovementDb = minimumMonitorImprovementDb;
    learningInfo.selectionObjective = selectionObjective;
    learningInfo.initialMonitorNmseDb = initialMonitorNmseDb;
    learningInfo.initialMonitorDpdPeak = initialMonitorPeak;
    learningInfo.activeCoefficientCount = coefficientCount;
    learningInfo.initialTrainingNmseDb = initialTrainingNmseDb;
    learningInfo.finalTrainingNmseDb = finalTrainingNmseDb;
    learningInfo.trainingImprovementDb = initialTrainingNmseDb - finalTrainingNmseDb;
    learningInfo.acceptedUpdateCount = nnz(history.accepted);
    learningInfo.blockHistory = history;
    learningInfo.coefficientHistory = coefficientHistory;
end

function [normalMatrix, gradient, usedRowCount, nmseDb, fullNmseDb, fullPeakMagnitude] = directNormalEquations(inputSignal, dpdCfg, paCfg, referenceGain, activeIndices, jacobianStride, selectionOptions)
    inputSignal = inputSignal(:);
    overlapLength = maximumModelDelay(dpdCfg) + maximumModelDelay(paCfg);
    sampledRows = (overlapLength + 1:jacobianStride:length(inputSignal)).';
    if isempty(sampledRows)
        error('dpdDirectLearn:InsufficientBlockSamples', 'A direct-learning block has no valid cascade samples.');
    end

    rowsPerChunk = round(configurationValue(dpdCfg, 'directJacobianRowsPerChunk', 512));
    validateattributes(rowsPerChunk, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    parameterCount = 2 * length(activeIndices);
    normalMatrix = zeros(parameterCount, parameterCount);
    gradient = zeros(parameterCount, 1);
    targetEnergy = 0;
    errorEnergy = 0;
    spectralGradientEnabled = selectionOptions.useSpectralObjective && (selectionOptions.spectralAclr1Weight > 0 || selectionOptions.spectralAclr2Weight > 0);
    fullNmseDb = NaN;
    fullPeakMagnitude = NaN;
    if spectralGradientEnabled
        dpdOutput = gmpCore(inputSignal, dpdCfg);
        paOutput = gmpCore(dpdOutput, paCfg);
        fullValidRows = overlapLength + 1:length(inputSignal);
        fullTarget = referenceGain * inputSignal(fullValidRows);
        fullResidual = paOutput(fullValidRows) - fullTarget;
        fullNmseDb = errorNmseDb(fullTarget, fullResidual);
        fullPeakMagnitude = max(abs(dpdOutput(fullValidRows)));
        [~, effectiveResidual] = dpdSpectralResidualObjective(fullResidual, paCfg.modelSampleRate, selectionOptions.channelBandwidth, selectionOptions.measurementBandwidth, selectionOptions.spectralAclr1Weight, selectionOptions.spectralAclr2Weight);
        clear dpdOutput paOutput fullTarget fullResidual
    else
        effectiveResidual = [];
    end

    for firstRowIndex = 1:rowsPerChunk:length(sampledRows)
        lastRowIndex = min(firstRowIndex + rowsPerChunk - 1, length(sampledRows));
        currentRows = sampledRows(firstRowIndex:lastRowIndex);
        extensionFirstSample = currentRows(1) - overlapLength;
        extendedInput = inputSignal(extensionFirstSample:currentRows(end));
        [outputSignal, derivativeB, derivativeConjugateB, localRows] = cascadeJacobian(extendedInput, dpdCfg, paCfg, activeIndices, jacobianStride);
        globalRows = localRows + extensionFirstSample - 1;
        if ~isequal(globalRows, currentRows)
            error('dpdDirectLearn:JacobianChunkAlignment', 'A direct-learning Jacobian chunk is misaligned.');
        end

        target = referenceGain * inputSignal(globalRows);
        residual = outputSignal - target;
        realParameterJacobian = derivativeB + derivativeConjugateB;
        imaginaryParameterJacobian = 1i * (derivativeB - derivativeConjugateB);
        complexJacobian = [realParameterJacobian imaginaryParameterJacobian];
        normalMatrix = normalMatrix + real(complexJacobian' * complexJacobian);
        if spectralGradientEnabled
            gradientResidual = effectiveResidual(globalRows - overlapLength);
        else
            gradientResidual = residual;
        end
        gradient = gradient + real(complexJacobian' * gradientResidual);
        targetEnergy = targetEnergy + sum(abs(target).^2);
        errorEnergy = errorEnergy + sum(abs(residual).^2);
    end

    normalMatrix = (normalMatrix + normalMatrix.') / 2;
    usedRowCount = length(sampledRows);
    if targetEnergy <= realmin || ~isfinite(targetEnergy) || ~isfinite(errorEnergy)
        nmseDb = NaN;
    else
        nmseDb = 10 * log10(max(errorEnergy / targetEnergy, realmin));
    end
end

function [candidateCfg, accepted, selectedLambda, selectedLineSearchFactor, bestCandidateNmseDb, acceptedNmseDb, acceptedPeak, bestCandidateAclr1Db, bestCandidateAclr2Db, acceptedAclr1Db, acceptedAclr2Db] = selectDirectUpdate(monitorInput, aclrMonitorInput, stressInput, currentCfg, paCfg, referenceGain, activeIndices, normalMatrix, gradient, usedRowCount, lambdaValues, lineSearchFactors, overlapLength, currentNmseDb, currentPeak, minimumMonitorImprovementDb, permittedOutputMagnitude, currentAclr1Db, currentAclr2Db, selectionOptions)
    activeCoefficientCount = length(activeIndices);
    diagonalFloor = configurationValue( currentCfg, 'directJacobianScaleFloor', 1e-12);
    normalizedNormalMatrix = normalMatrix / max(usedRowCount, 1);
    normalizedGradient = gradient / max(usedRowCount, 1);
    diagonalScale = sqrt(max( diag(normalizedNormalMatrix), diagonalFloor));
    scaledNormalMatrix = bsxfun(@rdivide, bsxfun(@rdivide, normalizedNormalMatrix, diagonalScale), diagonalScale.');
    scaledGradient = normalizedGradient ./ diagonalScale;
    candidateCfg = currentCfg;
    accepted = false;
    selectedLambda = NaN;
    selectedLineSearchFactor = 0;
    bestCandidateNmseDb = NaN;
    acceptedNmseDb = currentNmseDb;
    acceptedPeak = currentPeak;
    bestCandidateAclr1Db = NaN;
    bestCandidateAclr2Db = NaN;
    acceptedAclr1Db = currentAclr1Db;
    acceptedAclr2Db = currentAclr2Db;

    if any(~isfinite(normalizedNormalMatrix(:))) || any(~isfinite(normalizedGradient(:))) || any(~isfinite(diagonalScale(:))) || any(~isfinite(scaledNormalMatrix(:))) || any(~isfinite(scaledGradient(:)))
        return
    end

    [normalEigenvectors, normalEigenvalueMatrix] = eig(scaledNormalMatrix);
    normalEigenvalues = max(real(diag(normalEigenvalueMatrix)), 0);
    projectedGradient = normalEigenvectors.' * scaledGradient;

    if any(~isfinite(normalEigenvectors(:))) || any(~isfinite(normalEigenvalues(:))) || any(~isfinite(projectedGradient(:)))
        return
    end

    coefficientVector = currentCfg.coefficients(:);
    activeCoefficients = coefficientVector(activeIndices);
    maximumCoefficientNorm = configurationValue( currentCfg, 'maximumCoefficientNorm', Inf);
    bestCfg = currentCfg;
    bestLambda = NaN;
    bestLineSearchFactor = 0;
    bestPeak = NaN;
    maximumAclrTrialCount = length(lambdaValues) * length(lineSearchFactors);
    aclrTrialConfigs = cell(maximumAclrTrialCount, 1);
    aclrTrialNmseDb = NaN(maximumAclrTrialCount, 1);
    aclrTrialSpectralObjective = NaN(maximumAclrTrialCount, 1);
    aclrTrialLambda = NaN(maximumAclrTrialCount, 1);
    aclrTrialLineSearchFactor = NaN(maximumAclrTrialCount, 1);
    aclrTrialPeak = NaN(maximumAclrTrialCount, 1);
    aclrTrialCount = 0;

    for lambdaIndex = 1:length(lambdaValues)
        lambda = lambdaValues(lambdaIndex);
        scaledStep = -normalEigenvectors * (projectedGradient ./ (normalEigenvalues + lambda));
        parameterStep = scaledStep ./ diagonalScale;
        coefficientStep = parameterStep(1:activeCoefficientCount) + 1i * parameterStep(activeCoefficientCount+1:end);

        if any(~isfinite(coefficientStep))
            continue;
        end

        for scaleIndex = 1:length(lineSearchFactors)
            lineSearchFactor = lineSearchFactors(scaleIndex);
            trialActiveCoefficients = activeCoefficients + lineSearchFactor * coefficientStep;
            if any(~isfinite(trialActiveCoefficients)) || norm(trialActiveCoefficients) > maximumCoefficientNorm
                continue;
            end

            trialVector = complex(zeros(size(coefficientVector)));
            trialVector(activeIndices) = trialActiveCoefficients;
            trialCfg = currentCfg;
            trialCfg.coefficients = reshape( trialVector, size(currentCfg.coefficients));
            trialPeak = dpdPeak(stressInput, trialCfg, overlapLength);
            if ~candidateIsSafe(trialCfg, trialPeak, permittedOutputMagnitude)
                continue;
            end
            [trialNmseDb, ~, trialSpectralObjective] = blockCascadeNmseDb( monitorInput, trialCfg, paCfg, referenceGain, overlapLength, selectionOptions);

            if strcmp(selectionOptions.objective, 'aclr')
                if ~isfinite(trialNmseDb) || trialNmseDb > currentNmseDb + selectionOptions.maximumNmseDegradationDb
                    continue;
                end
                aclrTrialCount = aclrTrialCount + 1;
                aclrTrialConfigs{aclrTrialCount} = trialCfg;
                aclrTrialNmseDb(aclrTrialCount) = trialNmseDb;
                aclrTrialSpectralObjective(aclrTrialCount) = trialSpectralObjective;
                aclrTrialLambda(aclrTrialCount) = lambda;
                aclrTrialLineSearchFactor(aclrTrialCount) = lineSearchFactor;
                aclrTrialPeak(aclrTrialCount) = trialPeak;
            elseif isfinite(trialNmseDb) && (~isfinite(bestCandidateNmseDb) || trialNmseDb < bestCandidateNmseDb)
                    bestCandidateNmseDb = trialNmseDb;
                    bestCfg = trialCfg;
                    bestLambda = lambda;
                    bestLineSearchFactor = lineSearchFactor;
                    bestPeak = trialPeak;
            end
        end
    end

    if strcmp(selectionOptions.objective, 'aclr') && aclrTrialCount > 0
        [~, aclrTrialOrder] = sortrows([aclrTrialSpectralObjective(1:aclrTrialCount) aclrTrialNmseDb(1:aclrTrialCount)], [1 2]);
        aclrTrialOrder = aclrTrialOrder(1:min(selectionOptions.candidateCount, length(aclrTrialOrder)));
        for aclrTrialIndex = aclrTrialOrder(:).'
            trialCfg = aclrTrialConfigs{aclrTrialIndex};
            trialNmseDb = aclrTrialNmseDb(aclrTrialIndex);
            [trialAclr1Db, trialAclr2Db] = cascadeAclrDb(aclrMonitorInput, trialCfg, paCfg, overlapLength, selectionOptions);
            aclr2Accepted = isfinite(trialAclr2Db) && trialAclr2Db >= currentAclr2Db - selectionOptions.maximumAclr2DegradationDb;
            aclr1Better = isfinite(trialAclr1Db) && (~isfinite(bestCandidateAclr1Db) || trialAclr1Db > bestCandidateAclr1Db || (abs(trialAclr1Db - bestCandidateAclr1Db) <= 1e-12 && trialNmseDb < bestCandidateNmseDb));
            if aclr2Accepted && aclr1Better
                bestCandidateNmseDb = trialNmseDb;
                bestCandidateAclr1Db = trialAclr1Db;
                bestCandidateAclr2Db = trialAclr2Db;
                bestCfg = trialCfg;
                bestLambda = aclrTrialLambda(aclrTrialIndex);
                bestLineSearchFactor = aclrTrialLineSearchFactor(aclrTrialIndex);
                bestPeak = aclrTrialPeak(aclrTrialIndex);
            end
        end
    end

    nmseUpdateAccepted = strcmp(selectionOptions.objective, 'nmse') && isfinite(bestCandidateNmseDb) && currentNmseDb - bestCandidateNmseDb > minimumMonitorImprovementDb;
    aclrUpdateAccepted = strcmp(selectionOptions.objective, 'aclr') && isfinite(bestCandidateAclr1Db) && bestCandidateAclr1Db - currentAclr1Db > selectionOptions.minimumAclr1ImprovementDb;
    if nmseUpdateAccepted || aclrUpdateAccepted
        candidateCfg = bestCfg;
        accepted = true;
        selectedLambda = bestLambda;
        selectedLineSearchFactor = bestLineSearchFactor;
        acceptedNmseDb = bestCandidateNmseDb;
        acceptedPeak = bestPeak;
        acceptedAclr1Db = bestCandidateAclr1Db;
        acceptedAclr2Db = bestCandidateAclr2Db;
    end
end

function [worstAclr1Db, worstAclr2Db] = cascadeAclrDb(inputSignal, dpdCfg, paCfg, overlapLength, selectionOptions)
    dpdOutput = gmpCore(inputSignal, dpdCfg);
    paOutput = gmpCore(dpdOutput, paCfg);
    validRows = overlapLength + 1:length(inputSignal);
    [frequency, outputPsd] = welch(paOutput(validRows), paCfg.modelSampleRate, selectionOptions.spectrumSegmentLength);
    outputAclr = aclr(frequency, outputPsd, selectionOptions.channelBandwidth, selectionOptions.measurementBandwidth);
    worstAclr1Db = outputAclr.worstAclr1Db;
    worstAclr2Db = outputAclr.worstAclr2Db;
end

function [outputSignal, derivativeB, derivativeConjugateB, validRows] = cascadeJacobian(inputSignal, dpdCfg, paCfg, activeIndices, rowStride)
    inputSignal = inputSignal(:);
    dpdMemory = maximumModelDelay(dpdCfg);
    paMemory = maximumModelDelay(paCfg);
    firstValidSample = dpdMemory + paMemory + 1;
    validRows = (firstValidSample:rowStride:length(inputSignal)).';
    if isempty(validRows)
        error('dpdDirectLearn:InsufficientBlockSamples', 'A direct-learning block has no valid cascade samples.');
    end

    dpdOutput = gmpCore(inputSignal, dpdCfg);
    paOutput = gmpCore(dpdOutput, paCfg);
    outputSignal = paOutput(validRows);
    coefficientMask = false(size(dpdCfg.coefficients));
    coefficientMask(activeIndices) = true;
    coefficientCount = length(activeIndices);
    derivativeB = complex(zeros(length(validRows), coefficientCount));
    derivativeConjugateB = complex(zeros( length(validRows), coefficientCount));
    paDelays = unique([paCfg.signalDelays(:); paCfg.envelopeDelays(:)]).';
    directWeights = complex(zeros(length(validRows), length(paDelays)));
    conjugateWeights = complex(zeros(length(validRows), length(paDelays)));

    for signalIndex = 1:length(paCfg.signalDelays)
        signalDelay = paCfg.signalDelays(signalIndex);
        delayIndex = find(paDelays == signalDelay, 1);
        directWeights(:, delayIndex) = directWeights(:, delayIndex) + paCfg.coefficients(signalIndex, 1);
    end

    for orderIndex = 2:length(paCfg.orders)
        exponent = paCfg.orders(orderIndex) - 1;
        halfExponent = exponent / 2;
        for envelopeIndex = 1:length(paCfg.envelopeDelays)
            envelopeDelay = paCfg.envelopeDelays(envelopeIndex);
            envelopeDelayIndex = find(paDelays == envelopeDelay, 1);
            coefficientColumn = 2 + (orderIndex - 2) * length(paCfg.envelopeDelays) + envelopeIndex - 1;
            envelopeSamples = dpdOutput(validRows - envelopeDelay);
            envelopePower = abs(envelopeSamples).^exponent;
            envelopeDerivativePower = abs(envelopeSamples).^(exponent - 2);

            for signalIndex = 1:length(paCfg.signalDelays)
                signalDelay = paCfg.signalDelays(signalIndex);
                signalDelayIndex = find(paDelays == signalDelay, 1);
                coefficient = paCfg.coefficients( signalIndex, coefficientColumn);
                if coefficient == 0
                    continue
                end
                signalSamples = dpdOutput(validRows - signalDelay);

                directEnvelopeFactor = signalSamples .* halfExponent .* envelopeDerivativePower .* conj(envelopeSamples);
                conjugateEnvelopeFactor = signalSamples .* halfExponent .* envelopeDerivativePower .* envelopeSamples;

                directWeights(:, signalDelayIndex) = directWeights(:, signalDelayIndex) + coefficient .* envelopePower;
                directWeights(:, envelopeDelayIndex) = directWeights(:, envelopeDelayIndex) + coefficient .* directEnvelopeFactor;
                conjugateWeights(:, envelopeDelayIndex) = conjugateWeights(:, envelopeDelayIndex) + coefficient .* conjugateEnvelopeFactor;
            end
        end
    end

    for delayIndex = 1:length(paDelays)
        basisOutputRows = validRows - paDelays(delayIndex);
        dpdBasis = activeDpdBasisAtRows(inputSignal, dpdCfg, coefficientMask, basisOutputRows);
        directWeight = directWeights(:, delayIndex);
        conjugateWeight = conjugateWeights(:, delayIndex);
        if any(directWeight ~= 0)
            derivativeB = derivativeB + bsxfun(@times, dpdBasis, directWeight);
        end
        if any(conjugateWeight ~= 0)
            derivativeConjugateB = derivativeConjugateB + bsxfun(@times, conj(dpdBasis), conjugateWeight);
        end
    end
end

function basis = activeDpdBasisAtRows(inputSignal, modelCfg, coefficientMask, outputRows)
    outputRows = outputRows(:);
    orders = modelCfg.orders;
    signalDelays = modelCfg.signalDelays;
    envelopeDelays = modelCfg.envelopeDelays;
    maximumDelay = max([signalDelays envelopeDelays]);
    if any(outputRows < maximumDelay + 1) || any(outputRows > length(inputSignal))
        error('dpdDirectLearn:UnexpectedBasisAlignment', 'The DPD basis alignment is inconsistent with its memory.');
    end

    basis = complex(zeros(length(outputRows), nnz(coefficientMask)));
    activeColumnIndex = 1;
    coefficientColumn = 1;
    for signalIndex = 1:length(signalDelays)
        if coefficientMask(signalIndex, coefficientColumn)
            signalDelay = signalDelays(signalIndex);
            basis(:, activeColumnIndex) = inputSignal(outputRows - signalDelay);
            activeColumnIndex = activeColumnIndex + 1;
        end
    end

    for orderIndex = 2:length(orders)
        exponent = orders(orderIndex) - 1;
        for envelopeIndex = 1:length(envelopeDelays)
            coefficientColumn = 2 + (orderIndex - 2) * length(envelopeDelays) + envelopeIndex - 1;
            envelopeDelay = envelopeDelays(envelopeIndex);
            envelopeTerm = abs(inputSignal(outputRows - envelopeDelay)).^exponent;
            for signalIndex = 1:length(signalDelays)
                if coefficientMask(signalIndex, coefficientColumn)
                    signalDelay = signalDelays(signalIndex);
                    basis(:, activeColumnIndex) = inputSignal(outputRows - signalDelay) .* envelopeTerm;
                    activeColumnIndex = activeColumnIndex + 1;
                end
            end
        end
    end
end

function [nmseDb, peakMagnitude, spectralObjective] = blockCascadeNmseDb(inputSignal, dpdCfg, paCfg, referenceGain, overlapLength, selectionOptions)
    dpdOutput = gmpCore(inputSignal, dpdCfg);
    paOutput = gmpCore(dpdOutput, paCfg);
    validRows = overlapLength + 1:length(inputSignal);
    target = referenceGain * inputSignal(validRows);
    residual = paOutput(validRows) - target;
    nmseDb = errorNmseDb(target, residual);
    peakMagnitude = max(abs(dpdOutput(validRows)));
    if nargout >= 3
        if any(~isfinite(residual))
            spectralObjective = Inf;
        elseif nargin >= 6 && selectionOptions.useSpectralObjective
            spectralObjective = dpdSpectralResidualObjective(residual, paCfg.modelSampleRate, selectionOptions.channelBandwidth, selectionOptions.measurementBandwidth, selectionOptions.spectralAclr1Weight, selectionOptions.spectralAclr2Weight);
        else
            spectralObjective = sum(abs(residual).^2);
        end
    end
end

function nmseDb = cascadeNmseDb( inputSignal, dpdCfg, paCfg, referenceGain, overlapLength)
    inputSignal = inputSignal(:);
    chunkLength = round(configurationValue(dpdCfg, 'directMetricChunkLength', 30000));
    validateattributes(chunkLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    targetEnergy = 0;
    errorEnergy = 0;
    firstValidSample = overlapLength + 1;
    for firstSample = firstValidSample:chunkLength:length(inputSignal)
        lastSample = min(firstSample + chunkLength - 1, length(inputSignal));
        extensionFirstSample = firstSample - overlapLength;
        extendedInput = inputSignal(extensionFirstSample:lastSample);
        dpdOutput = gmpCore(extendedInput, dpdCfg);
        paOutput = gmpCore(dpdOutput, paCfg);
        localRows = overlapLength + 1:length(extendedInput);
        target = referenceGain * inputSignal(firstSample:lastSample);
        residual = paOutput(localRows) - target;
        targetEnergy = targetEnergy + sum(abs(target).^2);
        errorEnergy = errorEnergy + sum(abs(residual).^2);
    end
    if targetEnergy <= realmin || ~isfinite(targetEnergy) || ~isfinite(errorEnergy)
        nmseDb = NaN;
    else
        nmseDb = 10 * log10(max(errorEnergy / targetEnergy, realmin));
    end
end

function nmseDb = errorNmseDb(target, residual)
    targetEnergy = sum(abs(target).^2);
    errorEnergy = sum(abs(residual).^2);
    if targetEnergy <= realmin || ~isfinite(targetEnergy) || ~isfinite(errorEnergy)
        nmseDb = NaN;
    else
        nmseDb = 10 * log10(max( errorEnergy / targetEnergy, realmin));
    end
end

function peakMagnitude = dpdPeak(inputSignal, dpdCfg, overlapLength)
    dpdOutput = gmpCore(inputSignal, dpdCfg);
    validRows = overlapLength + 1:length(inputSignal);
    peakMagnitude = max(abs(dpdOutput(validRows)));
end

function accepted = candidateIsSafe(dpdCfg, outputPeak, permittedOutputMagnitude)
    accepted = isfinite(outputPeak) && outputPeak <= permittedOutputMagnitude && all(isfinite(dpdCfg.coefficients(:)));
end

function permittedPeak = permittedDpdPeak(dpdCfg, stressInput, overlapLength)
    validRows = overlapLength + 1:length(stressInput);
    inputPeak = max(abs(stressInput(validRows)));
    maximumPeakGainDb = configurationValue(dpdCfg, 'maximumDpdPeakGainDb', 12);
    maximumOutputMagnitude = configurationValue(dpdCfg, 'maximumOutputMagnitude', Inf);
    validateattributes(maximumPeakGainDb, {'numeric'}, {'scalar', 'real'});
    validateattributes(maximumOutputMagnitude, {'numeric'}, {'scalar', 'real', 'positive'});
    relativePeakLimit = inputPeak * 10^(maximumPeakGainDb / 20);
    permittedPeak = min(relativePeakLimit, maximumOutputMagnitude);
end

function [coefficients, initialization] = initialCoefficients(dpdCfg)
    initialization = configurationValue( dpdCfg, 'learningInitialization', 'identity');
    if ~ischar(initialization)
        error('dpdDirectLearn:InvalidInitialization', 'learningInitialization must be a character vector.');
    end

    initialization = lower(strtrim(initialization));
    if strcmp(initialization, 'configured')
        coefficients = dpdCfg.coefficients;
    elseif strcmp(initialization, 'identity')
        numberOfColumns = 1 + (length(dpdCfg.orders) - 1) * length(dpdCfg.envelopeDelays);
        coefficients = complex(zeros( length(dpdCfg.signalDelays), numberOfColumns));
        coefficients(1, 1) = 1;
    else
        error('dpdDirectLearn:UnknownInitialization', ['learningInitialization must be identity or ' 'configured.']);
    end
end

function history = initializeHistory(historyLength)
    history.epochIndex = zeros(historyLength, 1);
    history.blockIndex = zeros(historyLength, 1);
    history.firstSample = zeros(historyLength, 1);
    history.lastSample = zeros(historyLength, 1);
    history.sampleCount = zeros(historyLength, 1);
    history.jacobianRowCount = zeros(historyLength, 1);
    history.sampledNmseBeforeDb = NaN(historyLength, 1);
    history.nmseBeforeDb = NaN(historyLength, 1);
    history.nmseAfterDb = NaN(historyLength, 1);
    history.normalEquationsFinite = false(historyLength, 1);
    history.monitorNmseBeforeDb = NaN(historyLength, 1);
    history.monitorNmseCandidateDb = NaN(historyLength, 1);
    history.monitorNmseAfterDb = NaN(historyLength, 1);
    history.monitorAclr1BeforeDb = NaN(historyLength, 1);
    history.monitorAclr1CandidateDb = NaN(historyLength, 1);
    history.monitorAclr1AfterDb = NaN(historyLength, 1);
    history.monitorAclr2BeforeDb = NaN(historyLength, 1);
    history.monitorAclr2CandidateDb = NaN(historyLength, 1);
    history.monitorAclr2AfterDb = NaN(historyLength, 1);
    history.monitorDpdPeakBefore = NaN(historyLength, 1);
    history.monitorDpdPeakAfter = NaN(historyLength, 1);
    history.accepted = false(historyLength, 1);
    history.lambda = NaN(historyLength, 1);
    history.lineSearchFactor = zeros(historyLength, 1);
    history.coefficientNorm = NaN(historyLength, 1);
    history.dpdPeakBefore = NaN(historyLength, 1);
    history.dpdPeakAfter = NaN(historyLength, 1);
end

function history = trimHistory(history, historyLength)
    fieldNames = fieldnames(history);
    for fieldIndex = 1:length(fieldNames)
        fieldName = fieldNames{fieldIndex};
        history.(fieldName) = history.(fieldName)(1:historyLength);
    end
end

function validateInputs(modelInput, paCfg, dpdCfg, referenceGain)
    if ~isnumeric(modelInput) || ~isvector(modelInput) || isempty(modelInput) || any(~isfinite(modelInput(:)))
        error('dpdDirectLearn:InvalidModelInput', 'modelInput must be a finite nonempty numeric vector.');
    end
    if ~isstruct(paCfg) || ~isstruct(dpdCfg)
        error('dpdDirectLearn:InvalidConfiguration', 'paCfg and dpdCfg must be structures.');
    end
    if ~isnumeric(referenceGain) || ~isscalar(referenceGain) || ~isfinite(referenceGain) || abs(referenceGain) <= realmin
        error('dpdDirectLearn:InvalidReferenceGain', 'referenceGain must be a finite nonzero scalar.');
    end
end

function selectionObjective = resolveDirectSelectionObjective(dpdCfg)
    selectionObjective = configurationValue(dpdCfg, 'directSelectionObjective', 'nmse');
    if ~ischar(selectionObjective)
        error('dpdDirectLearn:InvalidSelectionObjective', 'directSelectionObjective must be a character vector.');
    end
    selectionObjective = lower(strtrim(selectionObjective));
    if ~strcmp(selectionObjective, 'nmse') && ~strcmp(selectionObjective, 'aclr')
        error('dpdDirectLearn:UnknownSelectionObjective', 'directSelectionObjective must be nmse or aclr.');
    end
end

function maximumDelay = maximumModelDelay(modelCfg)
    maximumDelay = max([ modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end

function requireFiniteCoefficients(coefficients)
    if isempty(coefficients) || any(~isfinite(coefficients(:)))
        error('dpdDirectLearn:NonfiniteCoefficients', 'Direct learning produced nonfinite coefficients.');
    end
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
