function [coefficients, fitInfo, fitState] = dpdBatchFit( paInput, normalizedPaOutput, dpdCfg, fitState, centerCoefficients)
    cumulativeRequested = nargin >= 4;
    if ~cumulativeRequested
        fitState = [];
    end
    if nargin < 5
        centerCoefficients = [];
    end
    paInput = paInput(:);
    normalizedPaOutput = normalizedPaOutput(:);

    sampleCount = min(length(paInput), length(normalizedPaOutput));
    paInput = paInput(1:sampleCount);
    normalizedPaOutput = normalizedPaOutput(1:sampleCount);

    coefficientMask = gmpDiagonalMask( dpdCfg, dpdCfg.diagonalCount);
    firstModelSample = maximumModelDelay(dpdCfg) + 1;
    target = paInput(firstModelSample:end);
    activeIndices = find(coefficientMask(:));
    centeredRegularization = ~isempty(centerCoefficients);
    if centeredRegularization
        if ~isnumeric(centerCoefficients) || ~isequal(size(centerCoefficients), size(coefficientMask)) || any(~isfinite(centerCoefficients(:)))
            error('dpdBatchFit:InvalidCenterCoefficients', 'centerCoefficients must be a finite numeric matrix matching the DPD coefficient matrix.');
        end
        centerActiveCoefficients = centerCoefficients(activeIndices);
    else
        centerActiveCoefficients = complex(zeros(length(activeIndices), 1));
    end

    rowCount = length(normalizedPaOutput) - firstModelSample + 1;
    if rowCount < 5
        error('dpdBatchFit:InsufficientSamples', ['At least five usable rows are required for the ' 'train/validation/test split.']);
    end
    trainingStop = max(1, floor(0.6 * rowCount));
    validationStop = max(trainingStop + 1, floor(0.8 * rowCount));
    validationStop = min(validationStop, rowCount - 1);

    trainingRows = 1:trainingStop;
    validationRows = trainingStop+1:validationStop;
    testRows = validationStop+1:rowCount;

    correlationChunkRows = round(configurationValue(dpdCfg, 'batchCorrelationChunkRows', 2048));
    validateattributes(correlationChunkRows, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    [trainingCorrelationRawCurrent, trainingCrossCorrelationRawCurrent, trainingTargetPowerCurrent] = rawStatistics(normalizedPaOutput, target, dpdCfg, coefficientMask, trainingRows, correlationChunkRows);
    [validationCorrelationRawCurrent, validationCrossCorrelationRawCurrent, validationTargetPowerCurrent] = rawStatistics(normalizedPaOutput, target, dpdCfg, coefficientMask, validationRows, correlationChunkRows);
    [testCorrelationRawCurrent, testCrossCorrelationRawCurrent, testTargetPowerCurrent] = rawStatistics(normalizedPaOutput, target, dpdCfg, coefficientMask, testRows, correlationChunkRows);

    if cumulativeRequested
        forgettingFactor = configurationValue(dpdCfg, 'indirectFitForgettingFactor', 1);
        validateForgettingFactor(forgettingFactor);
        [fitState, cumulativeStatistics] = accumulateFitState(fitState, dpdCfg, coefficientMask, forgettingFactor, trainingCorrelationRawCurrent, trainingCrossCorrelationRawCurrent, trainingTargetPowerCurrent, length(trainingRows), validationCorrelationRawCurrent, validationCrossCorrelationRawCurrent, validationTargetPowerCurrent, length(validationRows), testCorrelationRawCurrent, testCrossCorrelationRawCurrent, testTargetPowerCurrent, length(testRows));
        trainingCorrelationRaw = cumulativeStatistics.training.H;
        trainingCrossCorrelationRaw = cumulativeStatistics.training.g;
        validationCorrelationRaw = cumulativeStatistics.validation.H;
        validationCrossCorrelationRaw = cumulativeStatistics.validation.g;
        validationTargetPower = cumulativeStatistics.validation.E;
        testCorrelationRaw = cumulativeStatistics.test.H;
        testCrossCorrelationRaw = cumulativeStatistics.test.g;
        testTargetPower = cumulativeStatistics.test.E;
    else
        forgettingFactor = 1;
        trainingCorrelationRaw = trainingCorrelationRawCurrent;
        trainingCrossCorrelationRaw = trainingCrossCorrelationRawCurrent;
        validationCorrelationRaw = validationCorrelationRawCurrent;
        validationCrossCorrelationRaw = validationCrossCorrelationRawCurrent;
        validationTargetPower = validationTargetPowerCurrent;
        testCorrelationRaw = testCorrelationRawCurrent;
        testCrossCorrelationRaw = testCrossCorrelationRawCurrent;
        testTargetPower = testTargetPowerCurrent;
    end
    columnScale = correlationColumnScale(trainingCorrelationRaw);
    [trainingCorrelation, trainingCrossCorrelation] = normalizeStatistics(trainingCorrelationRaw, trainingCrossCorrelationRaw, columnScale);
    [validationCorrelation, validationCrossCorrelation] = normalizeStatistics(validationCorrelationRaw, validationCrossCorrelationRaw, columnScale);

    if isfield(dpdCfg, 'regularizationLambdaValues')
        lambdaValues = dpdCfg.regularizationLambdaValues(:).';
    else
        lambdaValues = logspace(-12, 2, 200);
    end

    [correlationEigenvectors, correlationEigenvalueMatrix] = eig(trainingCorrelation);
    correlationEigenvalues = max(real(diag(correlationEigenvalueMatrix)), 0);
    normalizedCenterCoefficients = centerActiveCoefficients .* columnScale.';
    projectedCrossCorrelation = correlationEigenvectors' * (trainingCrossCorrelation - trainingCorrelation * normalizedCenterCoefficients);
    validationNmseDb = zeros(size(lambdaValues));

    for lambdaIndex = 1:length(lambdaValues)
        normalizedCoefficients = normalizedCenterCoefficients + correlationEigenvectors * (projectedCrossCorrelation ./ (correlationEigenvalues + lambdaValues(lambdaIndex)));
        validationErrorPower = validationTargetPower - 2 * real(normalizedCoefficients' * validationCrossCorrelation) + real(normalizedCoefficients' * validationCorrelation * normalizedCoefficients);
        validationNmseDb(lambdaIndex) = 10 * log10(max(validationErrorPower / validationTargetPower, realmin));
    end

    [bestValidationNmseDb, bestLambdaIndex] = min(validationNmseDb);
    bestLambda = lambdaValues(bestLambdaIndex);

    fitRows = 1:validationStop;
    if cumulativeRequested && fitState.updateCount > 1
        fitRowCount = fitState.training.N + fitState.validation.N;
        fitCorrelationRaw = (fitState.training.H + fitState.validation.H) / fitRowCount;
        fitCrossCorrelationRaw = (fitState.training.g + fitState.validation.g) / fitRowCount;
    else
        fitRowCount = length(fitRows);
        fitCorrelationRaw = (length(trainingRows) * trainingCorrelationRaw + length(validationRows) * validationCorrelationRaw) / fitRowCount;
        fitCrossCorrelationRaw = (length(trainingRows) * trainingCrossCorrelationRaw + length(validationRows) * validationCrossCorrelationRaw) / fitRowCount;
    end
    fitColumnScale = correlationColumnScale(fitCorrelationRaw);
    [fitCorrelation, fitCrossCorrelation] = normalizeStatistics(fitCorrelationRaw, fitCrossCorrelationRaw, fitColumnScale);

    normalizedFitCenterCoefficients = centerActiveCoefficients .* fitColumnScale.';
    normalizedCoefficients = ( fitCorrelation + bestLambda * eye(size(fitCorrelation))) \ (fitCrossCorrelation + bestLambda * normalizedFitCenterCoefficients);
    activeCoefficients = normalizedCoefficients ./ fitColumnScale.';

    fullCoefficientVector = complex(zeros(numel(coefficientMask), 1));
    fullCoefficientVector(activeIndices) = activeCoefficients;
    coefficients = reshape(fullCoefficientVector, size(coefficientMask));

    testErrorPower = testTargetPower - 2 * real(activeCoefficients' * testCrossCorrelationRaw) + real(activeCoefficients' * testCorrelationRaw * activeCoefficients);
    testNmseDb = 10 * log10(max(testErrorPower / testTargetPower, realmin));

    fitInfo.firstModelSample = firstModelSample;
    fitInfo.sampleCount = rowCount;
    fitInfo.trainingRows = trainingRows;
    fitInfo.validationRows = validationRows;
    fitInfo.testRows = testRows;
    fitInfo.activeIndices = activeIndices;
    fitInfo.coefficientMask = coefficientMask;
    fitInfo.coefficientCount = length(activeIndices);
    fitInfo.lambdaValues = lambdaValues;
    fitInfo.validationNmseDb = validationNmseDb;
    fitInfo.bestLambda = bestLambda;
    fitInfo.bestValidationNmseDb = bestValidationNmseDb;
    fitInfo.testNmseDb = testNmseDb;
    fitInfo.columnScale = fitColumnScale;
    fitInfo.activeCoefficients = activeCoefficients;
    fitInfo.centeredRegularization = centeredRegularization;
    fitInfo.centerActiveCoefficients = centerActiveCoefficients;
    if cumulativeRequested
        fitInfo.fitMode = 'cumulative';
        fitInfo.forgettingFactor = forgettingFactor;
        fitInfo.fitStateUpdateCount = fitState.updateCount;
        fitInfo.accumulatedTrainingSampleCount = fitState.training.N;
        fitInfo.accumulatedValidationSampleCount = fitState.validation.N;
        fitInfo.accumulatedTestSampleCount = fitState.test.N;
    end
end

function [state, statistics] = accumulateFitState(state, modelCfg, coefficientMask, forgettingFactor, trainingH, trainingG, trainingE, trainingN, validationH, validationG, validationE, validationN, testH, testG, testE, testN)
    stateWasEmpty = isempty(state);
    if stateWasEmpty
        state = initializeFitState(modelCfg, coefficientMask);
    else
        validateFitState(state, modelCfg, coefficientMask);
    end

    state.training = accumulatePartition(state.training, trainingH, trainingG, trainingE, trainingN, forgettingFactor);
    state.validation = accumulatePartition(state.validation, validationH, validationG, validationE, validationN, forgettingFactor);
    state.test = accumulatePartition(state.test, testH, testG, testE, testN, forgettingFactor);
    state.updateCount = state.updateCount + 1;
    state.forgettingFactor = forgettingFactor;

    if stateWasEmpty
        statistics.training.H = trainingH;
        statistics.training.g = trainingG;
        statistics.training.E = trainingE;
        statistics.validation.H = validationH;
        statistics.validation.g = validationG;
        statistics.validation.E = validationE;
        statistics.test.H = testH;
        statistics.test.g = testG;
        statistics.test.E = testE;
    else
        statistics.training = partitionAverages(state.training);
        statistics.validation = partitionAverages(state.validation);
        statistics.test = partitionAverages(state.test);
    end
end

function state = initializeFitState(modelCfg, coefficientMask)
    coefficientCount = nnz(coefficientMask);
    state.version = 1;
    state.orders = modelCfg.orders(:).';
    state.signalDelays = modelCfg.signalDelays(:).';
    state.envelopeDelays = modelCfg.envelopeDelays(:).';
    state.diagonalCount = modelCfg.diagonalCount;
    state.coefficientMask = coefficientMask;
    state.activeIndices = find(coefficientMask(:));
    state.updateCount = 0;
    state.forgettingFactor = 1;
    state.training = emptyPartition(coefficientCount);
    state.validation = emptyPartition(coefficientCount);
    state.test = emptyPartition(coefficientCount);
end

function partition = emptyPartition(coefficientCount)
    partition.H = complex(zeros(coefficientCount, coefficientCount));
    partition.g = complex(zeros(coefficientCount, 1));
    partition.E = 0;
    partition.N = 0;
end

function partition = accumulatePartition(partition, correlation, crossCorrelation, targetPower, sampleCount, forgettingFactor)
    partition.H = forgettingFactor * partition.H + sampleCount * correlation;
    partition.H = (partition.H + partition.H') / 2;
    partition.g = forgettingFactor * partition.g + sampleCount * crossCorrelation;
    partition.E = forgettingFactor * partition.E + sampleCount * targetPower;
    partition.N = forgettingFactor * partition.N + sampleCount;
end

function statistics = partitionAverages(partition)
    statistics.H = partition.H / partition.N;
    statistics.H = (statistics.H + statistics.H') / 2;
    statistics.g = partition.g / partition.N;
    statistics.E = partition.E / partition.N;
end

function validateFitState(state, modelCfg, coefficientMask)
    requiredFields = {'version', 'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount', 'coefficientMask', 'activeIndices', 'updateCount', 'forgettingFactor', 'training', 'validation', 'test'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(state, requiredFields{fieldIndex})
            error('dpdBatchFit:InvalidFitState', 'The cumulative fit state must contain %s.', requiredFields{fieldIndex});
        end
    end
    configurationChanged = state.version ~= 1 || ~isequal(state.orders, modelCfg.orders(:).') || ~isequal(state.signalDelays, modelCfg.signalDelays(:).') || ~isequal(state.envelopeDelays, modelCfg.envelopeDelays(:).') || ~isequal(state.diagonalCount, modelCfg.diagonalCount) || ~isequal(state.coefficientMask, coefficientMask) || ~isequal(state.activeIndices, find(coefficientMask(:)));
    if configurationChanged
        error('dpdBatchFit:FitStateConfigurationChanged', 'The GMP configuration changed while a cumulative fit state was active.');
    end
    validateattributes(state.updateCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
    coefficientCount = nnz(coefficientMask);
    validatePartition(state.training, coefficientCount, 'training');
    validatePartition(state.validation, coefficientCount, 'validation');
    validatePartition(state.test, coefficientCount, 'test');
end

function validatePartition(partition, coefficientCount, description)
    requiredFields = {'H', 'g', 'E', 'N'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(partition, requiredFields{fieldIndex})
            error('dpdBatchFit:InvalidFitStatePartition', 'The %s partition must contain %s.', description, requiredFields{fieldIndex});
        end
    end
    validH = isnumeric(partition.H) && isequal(size(partition.H), [coefficientCount coefficientCount]) && all(isfinite(partition.H(:)));
    validG = isnumeric(partition.g) && isequal(size(partition.g), [coefficientCount 1]) && all(isfinite(partition.g(:)));
    validE = isnumeric(partition.E) && isscalar(partition.E) && isreal(partition.E) && isfinite(partition.E) && partition.E >= 0;
    validN = isnumeric(partition.N) && isscalar(partition.N) && isreal(partition.N) && isfinite(partition.N) && partition.N > 0;
    if ~validH || ~validG || ~validE || ~validN
        error('dpdBatchFit:InvalidFitStatePartition', 'The %s cumulative fit partition is invalid.', description);
    end
end

function validateForgettingFactor(forgettingFactor)
    valid = isnumeric(forgettingFactor) && isscalar(forgettingFactor) && isreal(forgettingFactor) && isfinite(forgettingFactor) && forgettingFactor == 1;
    if ~valid
        error('dpdBatchFit:InvalidForgettingFactor', 'indirectFitForgettingFactor must equal 1.');
    end
end

function [correlation, crossCorrelation, targetPower] = rawStatistics(inputSignal, target, modelCfg, coefficientMask, rows, chunkRowCount)
    coefficientCount = nnz(coefficientMask);
    rowCount = length(rows);
    correlation = complex(zeros(coefficientCount, coefficientCount));
    crossCorrelation = complex(zeros(coefficientCount, 1));
    targetEnergy = 0;
    for firstRowIndex = 1:chunkRowCount:rowCount
        lastRowIndex = min(firstRowIndex + chunkRowCount - 1, rowCount);
        currentRows = rows(firstRowIndex:lastRowIndex);
        basisChunk = modelBasisRows(inputSignal, modelCfg, coefficientMask, currentRows);
        correlation = correlation + basisChunk' * basisChunk;
        crossCorrelation = crossCorrelation + basisChunk' * target(currentRows);
        targetEnergy = targetEnergy + sum(abs(target(currentRows)).^2);
    end
    correlation = correlation / rowCount;
    correlation = (correlation + correlation') / 2;
    crossCorrelation = crossCorrelation / rowCount;
    targetPower = targetEnergy / rowCount;
end

function columnScale = correlationColumnScale(correlation)
    columnScale = max(sqrt(max(real(diag(correlation)), 0)).', sqrt(eps));
end

function [normalizedCorrelation, normalizedCrossCorrelation] = normalizeStatistics(correlation, crossCorrelation, columnScale)
    normalizedCorrelation = bsxfun(@rdivide, bsxfun(@rdivide, correlation, columnScale.'), columnScale);
    normalizedCorrelation = (normalizedCorrelation + normalizedCorrelation') / 2;
    normalizedCrossCorrelation = crossCorrelation ./ columnScale.';
end

function basis = modelBasisRows(inputSignal, modelCfg, coefficientMask, rows)
    if any(diff(rows) ~= 1)
        error('dpdBatchFit:NoncontiguousBasisRows', 'Basis rows must be contiguous.');
    end
    maximumDelay = maximumModelDelay(modelCfg);
    firstInputSample = rows(1);
    lastInputSample = rows(end) + maximumDelay;
    basis = gmpBasis(inputSignal(firstInputSample:lastInputSample), modelCfg, coefficientMask);
end

function maximumDelay = maximumModelDelay(modelCfg)
    maximumDelay = max([modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
