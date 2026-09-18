function results = paGmpPolynomialRun(fileName, maxPolynomialOrder, improvementThresholdDb, memoryOrder, alignmentMemoryOrder, showPlot, calculateTest)
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
        memoryOrder = 4;
    end

    if nargin < 5
        alignmentMemoryOrder = 30;
    end

    if nargin < 6
        showPlot = true;
    end

    if nargin < 7
        calculateTest = true;
    end

    candidateOrders = 1:2:maxPolynomialOrder;
    polynomialOrder = candidateOrders(2:end).';
    measuredData = paMeasuredCaptures(fileName);
    lambdaValues = logspace(-12, 2, 300);
    coefficientCount = zeros(size(polynomialOrder));
    bestLambda = zeros(size(polynomialOrder));
    crossValidationNmseDb = zeros(size(polynomialOrder));
    profiledDiagonalCount = zeros(size(polynomialOrder));
    diagonalProfiles = cell(size(polynomialOrder));

    for orderIndex = 1:length(polynomialOrder)
        currentOrders = candidateOrders(1:orderIndex + 1);
        [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'memory', currentOrders, memoryOrder, 0, alignmentMemoryOrder);
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

    selectedOrders = candidateOrders(1:selectedIndex + 1);
    selectedDiagonalCount = 0;
    selectedDiagonalProfile = [];
    [selectedBasis, selectedSampleIndices] = paBuildBasis(measuredData.inputBlock, 'memory', selectedOrders, memoryOrder, selectedDiagonalCount, alignmentMemoryOrder);
    selectedScore = paCrossValidate(selectedBasis, measuredData.outputBlocks(selectedSampleIndices, :), lambdaValues, calculateTest);
    summary = table(polynomialOrder, profiledDiagonalCount, coefficientCount, bestLambda, crossValidationNmseDb, nmseImprovementDb, 'VariableNames', {'PolynomialOrder', 'DiagonalCount', 'CoefficientCount', 'BestLambda', 'CrossValidationNMSE_dB', 'NMSEImprovement_dB'});
    disp(summary);

    if showPlot
        figure('Name', sprintf('GMP NMSE / polynomial order, M = %d', memoryOrder), 'NumberTitle', 'off', 'Color', 'w');
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
    results.model = 'GMP preliminary polynomial-order search on MP basis';
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryOrder + 1;
    results.signalDelays = 0:memoryOrder;
    results.envelopeDelays = 0:memoryOrder;
    results.diagonalCount = selectedDiagonalCount;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.maximumPolynomialOrder = maxPolynomialOrder;
    results.polynomialOrder = polynomialOrder;
    results.profiledDiagonalCount = profiledDiagonalCount;
    results.diagonalProfiles = diagonalProfiles;
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
    results.selectedDiagonalCount = selectedDiagonalCount;
    results.selectedDiagonalProfile = selectedDiagonalProfile;
    results.selectedNmseDb = crossValidationNmseDb(selectedIndex);
    results.selectedTestNmseDb = selectedScore.testNmseDb;
    results.selectedScore = selectedScore;
    results.bestScore = selectedScore;
    results.improvementThresholdDb = improvementThresholdDb;
    results.selectionRule = 'previous polynomial order before the first incremental NMSE improvement below the threshold';
    results.convergenceReached = convergenceReached;

    if convergenceReached
        results.stoppingPolynomialOrder = polynomialOrder(stoppingIndex);
    else
        results.stoppingPolynomialOrder = [];
    end
end
