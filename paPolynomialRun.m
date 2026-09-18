function results = paPolynomialRun(fileName, maxPolynomialOrder, improvementThresholdDb, memoryOrder, alignmentMemoryOrder, showPlot, candidateOrders, includeLinearModel, calculateTest)
    if nargin < 1
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2
        maxPolynomialOrder = 25;
    end

    if nargin < 3
        improvementThresholdDb = 0.20;
    end

    if nargin < 4
        memoryOrder = 0;
    end

    if nargin < 5
        alignmentMemoryOrder = memoryOrder;
    end

    if nargin < 6
        showPlot = true;
    end

    if nargin < 7
        candidateOrders = 1:2:maxPolynomialOrder;
    end

    if nargin < 8
        includeLinearModel = false;
    end

    if nargin < 9
        calculateTest = true;
    end

    candidateOrders = candidateOrders(:).';
    measuredData = paMeasuredCaptures(fileName);
    lambdaValues = logspace(-12, 2, 300);

    if includeLinearModel
        modelStartIndex = 1;
    else
        modelStartIndex = 2;
    end

    polynomialOrder = candidateOrders(modelStartIndex:end).';
    coefficientCount = zeros(size(polynomialOrder));
    bestLambda = zeros(size(polynomialOrder));
    crossValidationNmseDb = zeros(size(polynomialOrder));

    if memoryOrder == 0
        modelType = 'memoryless';
        modelName = 'Memoryless Polynomial';
    else
        modelType = 'memory';
        modelName = 'Memory Polynomial';
    end

    for orderIndex = 1:length(polynomialOrder)
        basisOrderCount = orderIndex + modelStartIndex - 1;
        currentOrders = candidateOrders(1:basisOrderCount);
        [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, modelType, currentOrders, memoryOrder, memoryOrder, alignmentMemoryOrder);
        score = paCrossValidate(basis, measuredData.outputBlocks(sampleIndices, :), lambdaValues, false);
        coefficientCount(orderIndex) = size(basis, 2);
        bestLambda(orderIndex) = score.bestLambda;
        crossValidationNmseDb(orderIndex) = score.validationNmseDb;
    end

    [bestNmseDb, bestIndex] = min(crossValidationNmseDb);
    nmseImprovementDb = [NaN; crossValidationNmseDb(1:end-1) - crossValidationNmseDb(2:end)];
    stoppingIndex = find(nmseImprovementDb < improvementThresholdDb, 1, 'first');
    convergenceReached = ~isempty(stoppingIndex);

    if convergenceReached
        selectedIndex = stoppingIndex - 1;
    else
        selectedIndex = length(polynomialOrder);
    end

    selectedOrderCount = selectedIndex + modelStartIndex - 1;
    selectedOrders = candidateOrders(1:selectedOrderCount);
    [selectedBasis, selectedSampleIndices] = paBuildBasis(measuredData.inputBlock, modelType, selectedOrders, memoryOrder, memoryOrder, alignmentMemoryOrder);
    selectedScore = paCrossValidate(selectedBasis, measuredData.outputBlocks(selectedSampleIndices, :), lambdaValues, calculateTest);
    summary = table(polynomialOrder, coefficientCount, bestLambda, crossValidationNmseDb, nmseImprovementDb, 'VariableNames', {'PolynomialOrder', 'CoefficientCount', 'BestLambda', 'ValidationNMSE_dB', 'NMSEImprovement_dB'});
    disp(summary);

    if showPlot
        figure('Name', sprintf('NMSE / polynomial order, M = %d', memoryOrder), 'NumberTitle', 'off', 'Color', 'w');
        plot(polynomialOrder, crossValidationNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(polynomialOrder(selectedIndex), crossValidationNmseDb(selectedIndex), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        grid on;
        xticks(polynomialOrder);
        xlabel('Polynomial order');
        ylabel('NMSE, dB');
    end

    fprintf('Minimum cross-validation NMSE = %.6f dB at polynomial order %d.\n', bestNmseDb, polynomialOrder(bestIndex));
    fprintf('Selected polynomial order = %d, cross-validation NMSE = %.6f dB.\n', polynomialOrder(selectedIndex), crossValidationNmseDb(selectedIndex));
    results.fileName = fileName;
    results.model = modelName;
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryOrder + 1;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.maximumPolynomialOrder = maxPolynomialOrder;
    results.candidateOrders = candidateOrders;
    results.polynomialOrder = polynomialOrder;
    results.coefficientCount = coefficientCount;
    results.lambdaValues = lambdaValues;
    results.bestLambda = bestLambda;
    results.crossValidationNmseDb = crossValidationNmseDb;
    results.validationNmseDb = crossValidationNmseDb;
    results.selectionNmseDb = crossValidationNmseDb;
    results.nmseImprovementDb = nmseImprovementDb;
    results.summary = summary;
    results.bestPolynomialOrder = polynomialOrder(bestIndex);
    results.bestNmseDb = bestNmseDb;
    results.selectedPolynomialOrder = polynomialOrder(selectedIndex);
    results.selectedOrders = selectedOrders;
    results.selectedNmseDb = crossValidationNmseDb(selectedIndex);
    results.selectedTestNmseDb = selectedScore.testNmseDb;
    results.selectedScore = selectedScore;
    results.bestScore = selectedScore;
    results.improvementThresholdDb = improvementThresholdDb;
    results.convergenceReached = convergenceReached;

    if convergenceReached
        results.stoppingPolynomialOrder = polynomialOrder(stoppingIndex);
    else
        results.stoppingPolynomialOrder = [];
    end
end
