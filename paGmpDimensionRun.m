function results = paGmpDimensionRun(fileName, maxPolynomialOrder, maxMemoryOrder, initialMemoryOrder, alignmentMemoryOrder, showPlots, calculatePrediction)
    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2 || isempty(maxPolynomialOrder)
        maxPolynomialOrder = 25;
    end

    if nargin < 3 || isempty(maxMemoryOrder)
        maxMemoryOrder = 12;
    end

    if nargin < 4 || isempty(initialMemoryOrder)
        initialMemoryOrder = min(4, maxMemoryOrder);
    end

    if nargin < 5 || isempty(alignmentMemoryOrder)
        alignmentMemoryOrder = max(30, maxMemoryOrder);
    end

    if nargin < 6 || isempty(showPlots)
        showPlots = true;
    end

    if nargin < 7 || isempty(calculatePrediction)
        calculatePrediction = true;
    end

    polynomialThresholdDb = 0.20;
    memoryThresholdDb = 0.05;
    diagonalThresholdDb = 0.05;
    alignmentMemoryOrder = max(alignmentMemoryOrder, maxMemoryOrder);
    dimensionResults = paDimensionRun(fileName, maxPolynomialOrder, maxMemoryOrder, polynomialThresholdDb, memoryThresholdDb, false);

    if ~dimensionResults.coordinateConvergenceReached
        error('paGmpDimensionRun:NoConvergence', 'The Memory Polynomial P and M search did not converge.');
    end

    selectedPolynomialOrder = dimensionResults.selectedPolynomialOrder;
    selectedOrders = dimensionResults.selectedOrders;
    selectedMemoryOrder = dimensionResults.selectedMemoryOrder;
    measuredData = paMeasuredCaptures(fileName);
    lambdaValues = logspace(-12, 2, 300);
    diagonalResults = paGmpDiagonalProfile(measuredData, selectedOrders, selectedMemoryOrder, alignmentMemoryOrder, lambdaValues, diagonalThresholdDb);
    selectedDiagonalCount = diagonalResults.selectedDiagonalCount;
    [selectedBasis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'gmp', selectedOrders, selectedMemoryOrder, selectedDiagonalCount, alignmentMemoryOrder);
    outputBlocks = measuredData.outputBlocks(sampleIndices, :);
    selectedScore = paCrossValidate(selectedBasis, outputBlocks, lambdaValues, calculatePrediction);

    if showPlots
        plotPolynomialStudy(dimensionResults.polynomial, selectedPolynomialOrder);
        plotMemoryStudy(dimensionResults.memory, selectedMemoryOrder, selectedPolynomialOrder);
        plotDiagonalStudy(diagonalResults, selectedDiagonalCount, selectedPolynomialOrder, selectedMemoryOrder);
    end

    results.fileName = fileName;
    results.model = 'Generalized Memory Polynomial';
    results.searchMethod = 'P and M from Memory Polynomial basis, then d from fixed P and M GMP profile';
    results.dimensionBasis = 'Memory Polynomial';
    results.diagonalBasis = 'Generalized Memory Polynomial';
    results.polynomialThresholdDb = polynomialThresholdDb;
    results.memoryThresholdDb = memoryThresholdDb;
    results.diagonalThresholdDb = diagonalThresholdDb;
    results.maximumPolynomialOrder = maxPolynomialOrder;
    results.maximumMemoryOrder = maxMemoryOrder;
    results.requestedInitialMemoryOrder = initialMemoryOrder;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.dimension = dimensionResults;
    results.polynomial = dimensionResults.polynomial;
    results.memory = dimensionResults.memory;
    results.diagonal = diagonalResults;
    results.iterationHistory = dimensionResults.iterationHistory;
    results.coordinateIterationCount = dimensionResults.coordinateIterationCount;
    results.coordinateConvergenceReached = dimensionResults.coordinateConvergenceReached;
    results.polynomialOrder = selectedPolynomialOrder;
    results.selectedPolynomialOrder = selectedPolynomialOrder;
    results.orders = selectedOrders;
    results.selectedOrders = selectedOrders;
    results.memoryOrder = selectedMemoryOrder;
    results.selectedMemoryOrder = selectedMemoryOrder;
    results.memoryDepth = selectedMemoryOrder + 1;
    results.signalDelays = 0:selectedMemoryOrder;
    results.envelopeDelays = 0:selectedMemoryOrder;
    results.diagonalCount = selectedDiagonalCount;
    results.sideDiagonalCount = selectedDiagonalCount;
    results.diagonalSelection = 'smallest d within 0.05 dB of the minimum repetition 5 validation NMSE';
    results.sampleIndices = sampleIndices;
    results.coefficientCount = size(selectedBasis, 2);
    results.bestLambda = selectedScore.bestLambda;
    results.validationNmseDb = selectedScore.validationNmseDb;
    results.selectedNmseDb = selectedScore.validationNmseDb;
    results.bestNmseDb = selectedScore.validationNmseDb;
    results.selectedValidationNmseDb = selectedScore.validationNmseDb;
    results.selectedTestNmseDb = selectedScore.testNmseDb;
    results.bestTestNmseDb = selectedScore.testNmseDb;
    results.selectedScore = selectedScore;
    results.bestScore = selectedScore;
    results.trainingRepetitions = 1:4;
    results.validationRepetition = 5;
    results.splitMethod = '80/20 holdout';

    fprintf('Selected GMP polynomial order = %d.\n', selectedPolynomialOrder);
    fprintf('Selected GMP memory order = %d.\n', selectedMemoryOrder);
    fprintf('Selected GMP diagonal count = %d.\n', selectedDiagonalCount);
    fprintf('Repetition 5 validation NMSE = %.6f dB.\n', selectedScore.validationNmseDb);
end

function plotPolynomialStudy(polynomialResults, selectedPolynomialOrder)
    figure('Name', sprintf('GMP dimensions / polynomial order, M = %d', polynomialResults.memoryOrder), 'NumberTitle', 'off', 'Color', 'w');
    plot(polynomialResults.polynomialOrder, polynomialResults.selectionNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    plot(selectedPolynomialOrder, polynomialResults.selectedNmseDb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    grid on;
    xticks(polynomialResults.polynomialOrder);
    xlabel('Polynomial order');
    ylabel('NMSE, dB');
end

function plotMemoryStudy(memoryResults, selectedMemoryOrder, selectedPolynomialOrder)
    figure('Name', sprintf('GMP dimensions / memory order, P = %d', selectedPolynomialOrder), 'NumberTitle', 'off', 'Color', 'w');
    plot(memoryResults.memoryOrder, memoryResults.selectionNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    plot(selectedMemoryOrder, memoryResults.selectedNmseDb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    grid on;
    xticks(memoryResults.memoryOrder);
    xlabel('Memory');
    ylabel('NMSE, dB');
end

function plotDiagonalStudy(diagonalResults, selectedDiagonalCount, selectedPolynomialOrder, selectedMemoryOrder)
    figure('Name', sprintf('GMP dimensions / diagonal count, P = %d, M = %d', selectedPolynomialOrder, selectedMemoryOrder), 'NumberTitle', 'off', 'Color', 'w');
    plot(diagonalResults.diagonalCount, diagonalResults.validationNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    plot(selectedDiagonalCount, diagonalResults.selectedNmseDb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    grid on;
    xticks(diagonalResults.diagonalCount);
    xlabel('Diagonal count');
    ylabel('NMSE, dB');
end
