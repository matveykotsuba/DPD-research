function results = paPolynomialSpectrumRun(fileName)
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

    memoryOrder = 4;
    alignmentMemoryOrder = 30;
    polynomialOrder = (13:2:23).';
    candidateOrders = 1:2:polynomialOrder(end);
    paCfg = paConfig(false);
    sampleRate = paCfg.modelSampleRate;
    segmentLength = 1024;

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
    selectionRows = validationStop+1:modelSampleCount;

    lambdaValues = logspace(-12, 2, 300);
    bestLambda = zeros(size(polynomialOrder));
    modelOutput = cell(length(polynomialOrder), 1);

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

        fitRows = 1:validationStop;
        fitBasis = currentBasis(fitRows, :);
        fitOutput = modelOutputBlocks(fitRows, :);
        fitColumnScale = sqrt(mean(abs(fitBasis).^2, 1));
        normalizedFitBasis = fitBasis ./ fitColumnScale;
        normalizedSelectionBasis = currentBasis(selectionRows, :) ./ fitColumnScale;
        fitCorrelation = normalizedFitBasis' * normalizedFitBasis / length(fitRows);
        fitCrossCorrelation = normalizedFitBasis' * mean(fitOutput, 2) / length(fitRows);
        identityMatrix = eye(size(fitCorrelation));
        fittedCoefficients = (fitCorrelation + bestLambda(orderIndex) * identityMatrix) \ fitCrossCorrelation;
        modelOutput{orderIndex} = normalizedSelectionBasis * fittedCoefficients;
    end

    [frequency, firstModelPsd] = welch(modelOutput{1}, sampleRate, segmentLength);
    modelPsd = zeros(length(firstModelPsd), length(polynomialOrder));
    modelPsd(:, 1) = firstModelPsd;

    for orderIndex = 2:length(polynomialOrder)
        [~, modelPsd(:, orderIndex)] = welch(modelOutput{orderIndex}, sampleRate, segmentLength);
    end

    measuredPsd = zeros(size(firstModelPsd));

    for blockIndex = 1:blockCount
        [~, measuredBlockPsd] = welch(modelOutputBlocks(selectionRows, blockIndex), sampleRate, segmentLength);
        measuredPsd = measuredPsd + measuredBlockPsd;
    end

    measuredPsd = measuredPsd / blockCount;
    psdReference = max(measuredPsd);
    modelPsdDb = 10 * log10(modelPsd / psdReference + eps);
    measuredPsdDb = 10 * log10(measuredPsd / psdReference + eps);

    figure('Name', 'Memory Polynomial spectra / polynomial order', 'NumberTitle', 'off', 'Color', 'w');
    colors = lines(length(polynomialOrder));
    legendText = cell(length(polynomialOrder) + 1, 1);

    for orderIndex = 1:length(polynomialOrder)
        plot(frequency / 1e6, modelPsdDb(:, orderIndex), 'Color', colors(orderIndex, :), 'LineWidth', 1.2);
        hold on;
        legendText{orderIndex} = sprintf('P = %d', polynomialOrder(orderIndex));
    end

    plot(frequency / 1e6, measuredPsdDb, 'k', 'LineWidth', 1.8);
    legendText{end} = 'Real PA output';

    grid on;
    xlabel('Frequency, MHz');
    ylabel(' PSD, dB');
    legend(legendText, 'Location', 'best');
    xlim([-100 100]);
    ylim([-80 5]);

    results.fileName = fileName;
    results.model = 'Memory Polynomial';
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryOrder + 1;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.sampleRate = sampleRate;
    results.segmentLength = segmentLength;
    results.blockCount = blockCount;
    results.polynomialOrder = polynomialOrder;
    results.frequency = frequency;
    results.modelPsd = modelPsd;
    results.modelPsdDb = modelPsdDb;
    results.measuredPsd = measuredPsd;
    results.measuredPsdDb = measuredPsdDb;
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
