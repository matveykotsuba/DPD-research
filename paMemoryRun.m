function results = paMemoryRun(fileName, maxMemoryOrder, orders, improvementThresholdDb, showPlot, calculateTest)
    if nargin < 1
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2
        maxMemoryOrder = 30;
    end

    if nargin < 3
        paCfg = paConfig(false);
        orders = paCfg.orders;
    end

    if nargin < 4
        improvementThresholdDb = 0.05;
    end

    if nargin < 5
        showPlot = true;
    end

    if nargin < 6
        calculateTest = true;
    end

    orders = orders(:).';
    measuredData = paMeasuredCaptures(fileName);
    lambdaValues = logspace(-12, 2, 300);
    memoryOrder = (0:maxMemoryOrder).';
    memoryDepth = memoryOrder + 1;
    coefficientCount = zeros(size(memoryOrder));
    bestLambda = zeros(size(memoryOrder));
    crossValidationNmseDb = zeros(size(memoryOrder));

    for memoryIndex = 1:length(memoryOrder)
        currentMemoryOrder = memoryOrder(memoryIndex);
        [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'memory', orders, currentMemoryOrder, currentMemoryOrder, maxMemoryOrder);
        score = paCrossValidate(basis, measuredData.outputBlocks(sampleIndices, :), lambdaValues, false);
        coefficientCount(memoryIndex) = size(basis, 2);
        bestLambda(memoryIndex) = score.bestLambda;
        crossValidationNmseDb(memoryIndex) = score.validationNmseDb;
    end

    [bestNmseDb, bestIndex] = min(crossValidationNmseDb);
    nmseImprovementDb = [NaN; crossValidationNmseDb(1:end-1) - crossValidationNmseDb(2:end)];
    stoppingIndex = find(nmseImprovementDb < improvementThresholdDb, 1, 'first');
    convergenceReached = ~isempty(stoppingIndex);

    if convergenceReached
        selectedIndex = stoppingIndex - 1;
    else
        selectedIndex = length(memoryOrder);
    end

    selectedMemoryOrder = memoryOrder(selectedIndex);
    [selectedBasis, selectedSampleIndices] = paBuildBasis(measuredData.inputBlock, 'memory', orders, selectedMemoryOrder, selectedMemoryOrder, maxMemoryOrder);
    selectedScore = paCrossValidate(selectedBasis, measuredData.outputBlocks(selectedSampleIndices, :), lambdaValues, calculateTest);
    summary = table(memoryOrder,  coefficientCount, bestLambda, crossValidationNmseDb, nmseImprovementDb, 'VariableNames', {'MemoryOrder',  'CoefficientCount', 'BestLambda', 'ValidationNMSE_dB', 'NMSEImprovement_dB'});
    disp(summary);

    if showPlot
        figure('Name', sprintf('NMSE / memory order, P = %d', orders(end)), 'NumberTitle', 'off', 'Color', 'w');
        plot(memoryOrder, crossValidationNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(selectedMemoryOrder, crossValidationNmseDb(selectedIndex), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        grid on;
        xticks(0:2:maxMemoryOrder);
        xlabel('Memory');
        ylabel('NMSE, dB');
    end

    fprintf('Minimum cross-validation NMSE = %.6f dB at M = %d.\n', bestNmseDb, memoryOrder(bestIndex));
    fprintf('Selected memory order = %d, cross-validation NMSE = %.6f dB.\n', selectedMemoryOrder, crossValidationNmseDb(selectedIndex));
    results.fileName = fileName;
    results.model = 'Memory Polynomial';
    results.orders = orders;
    results.maximumMemoryOrder = maxMemoryOrder;
    results.blockCount = measuredData.blockCount;
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryDepth;
    results.coefficientCount = coefficientCount;
    results.lambdaValues = lambdaValues;
    results.bestLambda = bestLambda;
    results.crossValidationNmseDb = crossValidationNmseDb;
    results.validationNmseDb = crossValidationNmseDb;
    results.selectionNmseDb = crossValidationNmseDb;
    results.nmseImprovementDb = nmseImprovementDb;
    results.summary = summary;
    results.bestMemoryOrder = memoryOrder(bestIndex);
    results.bestNmseDb = bestNmseDb;
    results.selectedMemoryOrder = selectedMemoryOrder;
    results.selectedMemoryDepth = selectedMemoryOrder + 1;
    results.selectedMemoryDelays = 0:selectedMemoryOrder;
    results.selectedNmseDb = crossValidationNmseDb(selectedIndex);
    results.selectedTestNmseDb = selectedScore.testNmseDb;
    results.selectedScore = selectedScore;
    results.bestScore = selectedScore;
    results.improvementThresholdDb = improvementThresholdDb;
    results.convergenceReached = convergenceReached;

    if convergenceReached
        results.stoppingMemoryOrder = memoryOrder(stoppingIndex);
    else
        results.stoppingMemoryOrder = [];
    end
end
