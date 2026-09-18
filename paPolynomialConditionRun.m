function results = paPolynomialConditionRun(fileName)
    if nargin < 1
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    data = load(fileName);

    if isfield(data, 'signal_in')
        inputSignal = data.signal_in(:);
        measuredOutput = data.signal_out(:);
    else
        inputSignal = data.txSignal(:);
        measuredOutput = data.rxSignal(:);
    end

    memoryOrder = 7;
    alignmentMemoryOrder = 30;
    polynomialOrder = (13:2:23).';
    candidateOrders = 1:2:polynomialOrder(end);

    sampleCount = min(length(inputSignal), length(measuredOutput));
    blockCount = 5;
    blockLength = floor(sampleCount / blockCount);
    sampleCount = blockCount * blockLength;

    inputBlocks = reshape(inputSignal(1:sampleCount), blockLength, blockCount);
    measuredOutputBlocks = reshape(measuredOutput(1:sampleCount), blockLength, blockCount);

    inputBlock = inputBlocks(:, 1);
    firstModelSample = alignmentMemoryOrder + 1;
    modelSampleIndices = (firstModelSample:blockLength).';
    maximumBasis = memoryPolynomialBasis(inputBlock, candidateOrders, memoryOrder, modelSampleIndices);
    modelOutputBlocks = measuredOutputBlocks(firstModelSample:end, :);

    modelSampleCount = size(maximumBasis, 1);
    trainingStop = floor(0.6 * modelSampleCount);
    validationStop = floor(0.8 * modelSampleCount);
    trainingRows = 1:trainingStop;
    validationRows = trainingStop+1:validationStop;
    fitRows = 1:validationStop;

    lambdaValues = logspace(-12, 2, 300);
    conditionNumber = zeros(size(polynomialOrder));
    bestLambda = zeros(size(polynomialOrder));

    for orderIndex = 1:length(polynomialOrder)
        basisOrderCount = (polynomialOrder(orderIndex) + 1) / 2;
        coefficientCount = basisOrderCount * (memoryOrder + 1);
        currentBasis = maximumBasis(:, 1:coefficientCount);

        trainingBasis = currentBasis(trainingRows, :);
        validationBasis = currentBasis(validationRows, :);
        trainingOutput = modelOutputBlocks(trainingRows, :);
        validationOutput = modelOutputBlocks(validationRows, :);

        columnScale = sqrt(mean(abs(trainingBasis).^2, 1));
        normalizedTrainingBasis = trainingBasis ./ columnScale;
        normalizedValidationBasis = validationBasis ./ columnScale;

        trainingCorrelation = normalizedTrainingBasis' * normalizedTrainingBasis / length(trainingRows);
        trainingCrossCorrelation = normalizedTrainingBasis' * mean(trainingOutput, 2) / length(trainingRows);
        validationCorrelation = normalizedValidationBasis' * normalizedValidationBasis / length(validationRows);
        validationCrossCorrelation = normalizedValidationBasis' * mean(validationOutput, 2) / length(validationRows);
        validationOutputPower = mean(abs(validationOutput(:)).^2);

        trainingCorrelation = (trainingCorrelation + trainingCorrelation') / 2;
        [eigenvectors, eigenvalueMatrix] = eig(trainingCorrelation);
        eigenvalues = max(real(diag(eigenvalueMatrix)), 0);
        coefficientCoordinates = (eigenvectors' * trainingCrossCorrelation) ./ (eigenvalues + lambdaValues);
        normalizedCoefficients = eigenvectors * coefficientCoordinates;
        validationErrorPower = validationOutputPower - 2 * real(sum(conj(normalizedCoefficients) .* validationCrossCorrelation, 1)) + real(sum(conj(normalizedCoefficients) .* (validationCorrelation * normalizedCoefficients), 1));
        validationNmse = max(validationErrorPower / validationOutputPower, realmin);

        [~, bestLambdaIndex] = min(validationNmse);
        bestLambda(orderIndex) = lambdaValues(bestLambdaIndex);

        fitBasis = currentBasis(fitRows, :);
        fitColumnScale = sqrt(mean(abs(fitBasis).^2, 1));
        normalizedFitBasis = fitBasis ./ fitColumnScale;
        fitCorrelation = normalizedFitBasis' * normalizedFitBasis / length(fitRows);
        fitCorrelation = (fitCorrelation + fitCorrelation') / 2;
        conditionNumber(orderIndex) = cond(fitCorrelation);
    end

    figure('Name', 'Memory Polynomial conditioning', 'NumberTitle', 'off', 'Color', 'w');

    subplot(2, 1, 1);
    semilogy(polynomialOrder, conditionNumber, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xticks(polynomialOrder);
    xlabel('Polynomial order');
    ylabel('matrix conditioning');

    subplot(2, 1, 2);
    semilogy(polynomialOrder, bestLambda, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xticks(polynomialOrder);
    xlabel('Polynomial order');
    ylabel('Best lambda');

    results.fileName = fileName;
    results.model = 'Memory Polynomial';
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryOrder + 1;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.polynomialOrder = polynomialOrder;
    results.conditionNumber = conditionNumber;
    results.bestLambda = bestLambda;
end

function basis = memoryPolynomialBasis(inputSignal, orders, memoryOrder, sampleIndices)
    delayCount = memoryOrder + 1;
    coefficientCount = length(orders) * delayCount;
    basis = complex(zeros(length(sampleIndices), coefficientCount));
    columnIndex = 1;

    for orderIndex = 1:length(orders)
        exponent = orders(orderIndex) - 1;

        for delay = 0:memoryOrder
            delayedSignal = inputSignal(sampleIndices - delay);
            basis(:, columnIndex) = delayedSignal .* abs(delayedSignal).^exponent;
            columnIndex = columnIndex + 1;
        end
    end
end
