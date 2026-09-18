function results = paDimensionRun(fileName, maxPolynomialOrder, maxMemoryOrder, nonlinearityThresholdDb, memoryThresholdDb, showPlots)
    if nargin < 1
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2
        maxPolynomialOrder = 25;
    end

    if nargin < 3
        maxMemoryOrder = 30;
    end

    if nargin < 4
        nonlinearityThresholdDb = 0.20;
    end

    if nargin < 5
        memoryThresholdDb = 0.05;
    end

    if nargin < 6
        showPlots = true;
    end

    initialPolynomialResults = paPolynomialRun(fileName, maxPolynomialOrder, nonlinearityThresholdDb, 0, maxMemoryOrder, false, 1:2:maxPolynomialOrder, false, false);

    if ~initialPolynomialResults.convergenceReached
        error('paDimensionRun:PolynomialBoundary', 'The polynomial-order stopping threshold was not reached. Increase maxPolynomialOrder.');
    end

    initialMemoryResults = paMemoryRun(fileName, maxMemoryOrder, initialPolynomialResults.selectedOrders, memoryThresholdDb, false, false);

    if ~initialMemoryResults.convergenceReached
        error('paDimensionRun:MemoryBoundary', 'The memory-order stopping threshold was not reached. Increase maxMemoryOrder.');
    end
    currentPolynomialOrder = initialPolynomialResults.selectedPolynomialOrder;
    currentMemoryOrder = initialMemoryResults.selectedMemoryOrder;
    polynomialResults = initialPolynomialResults;
    memoryResults = initialMemoryResults;
    maximumIterationCount = 10;
    polynomialOrderHistory = zeros(maximumIterationCount + 1, 1);
    memoryOrderHistory = zeros(maximumIterationCount + 1, 1);
    polynomialOrderHistory(1) = currentPolynomialOrder;
    memoryOrderHistory(1) = currentMemoryOrder;
    iterationCount = 0;
    convergenceReached = false;

    for iterationIndex = 1:maximumIterationCount
        polynomialResults = paPolynomialRun(fileName, maxPolynomialOrder, nonlinearityThresholdDb, currentMemoryOrder, maxMemoryOrder, false, 1:2:maxPolynomialOrder, false, false);

        if ~polynomialResults.convergenceReached
            error('paDimensionRun:PolynomialBoundary', 'The polynomial-order stopping threshold was not reached. Increase maxPolynomialOrder.');
        end

        memoryResults = paMemoryRun(fileName, maxMemoryOrder, polynomialResults.selectedOrders, memoryThresholdDb, false, false);

        if ~memoryResults.convergenceReached
            error('paDimensionRun:MemoryBoundary', 'The memory-order stopping threshold was not reached. Increase maxMemoryOrder.');
        end
        newPolynomialOrder = polynomialResults.selectedPolynomialOrder;
        newMemoryOrder = memoryResults.selectedMemoryOrder;
        iterationCount = iterationIndex;
        polynomialOrderHistory(iterationIndex + 1) = newPolynomialOrder;
        memoryOrderHistory(iterationIndex + 1) = newMemoryOrder;

        if newPolynomialOrder == currentPolynomialOrder && newMemoryOrder == currentMemoryOrder
            convergenceReached = true;
            break
        end

        currentPolynomialOrder = newPolynomialOrder;
        currentMemoryOrder = newMemoryOrder;
    end

    if ~convergenceReached
        error('paDimensionRun:NoConvergence', 'The MP P and M search did not converge in %d iterations.', maximumIterationCount);
    end

    iteration = (0:iterationCount).';
    polynomialOrderHistory = polynomialOrderHistory(1:iterationCount + 1);
    memoryOrderHistory = memoryOrderHistory(1:iterationCount + 1);
    iterationHistory = table(iteration, polynomialOrderHistory, memoryOrderHistory, 'VariableNames', {'Iteration', 'PolynomialOrder', 'MemoryOrder'});
    disp(iterationHistory);

    if showPlots
        figure('Name', 'NMSE / polynomial order, Memoryless Polynomial', 'NumberTitle', 'off', 'Color', 'w');
        plot(initialPolynomialResults.polynomialOrder, initialPolynomialResults.selectionNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(initialPolynomialResults.selectedPolynomialOrder, initialPolynomialResults.selectedNmseDb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        grid on;
        xticks(initialPolynomialResults.polynomialOrder);
        xlabel('Polynomial order');
        ylabel('NMSE, dB');

        figure('Name', sprintf('NMSE / polynomial order, Memory Polynomial, M = %d', polynomialResults.memoryOrder), 'NumberTitle', 'off', 'Color', 'w');
        plot(polynomialResults.polynomialOrder, polynomialResults.selectionNmseDb, 'bs-', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(polynomialResults.selectedPolynomialOrder, polynomialResults.selectedNmseDb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        grid on;
        xticks(polynomialResults.polynomialOrder);
        xlabel('Polynomial order');
        ylabel('NMSE, dB');

        figure('Name', sprintf('NMSE / memory order, Memory Polynomial, P = %d', polynomialResults.selectedPolynomialOrder), 'NumberTitle', 'off', 'Color', 'w');
        plot(memoryResults.memoryOrder, memoryResults.selectionNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(memoryResults.selectedMemoryOrder, memoryResults.selectedNmseDb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        grid on;
        xticks(0:2:maxMemoryOrder);
        xlabel('Memory');
        ylabel('NMSE, dB');
    end

    results.fileName = fileName;
    results.initialPolynomial = initialPolynomialResults;
    results.initialMemory = initialMemoryResults;
    results.polynomial = polynomialResults;
    results.memory = memoryResults;
    results.iterationHistory = iterationHistory;
    results.coordinateIterationCount = iterationCount;
    results.coordinateConvergenceReached = convergenceReached;
    results.memorylessPolynomialOrder = initialPolynomialResults.selectedPolynomialOrder;
    results.memorylessOrders = initialPolynomialResults.selectedOrders;
    results.memoryPolynomialOrder = polynomialResults.selectedPolynomialOrder;
    results.memoryPolynomialOrders = polynomialResults.selectedOrders;
    results.memoryOrder = memoryResults.selectedMemoryOrder;
    results.memoryDepth = memoryResults.selectedMemoryDepth;
    results.selectedPolynomialOrder = results.memoryPolynomialOrder;
    results.selectedOrders = results.memoryPolynomialOrders;
    results.selectedMemoryOrder = results.memoryOrder;
    results.selectedMemoryDepth = results.memoryDepth;

    fprintf('Selected Memoryless model: polynomial order %d.\n', results.memorylessPolynomialOrder);
    fprintf('Selected MP model: polynomial order %d, memory order %d, memory depth %d samples.\n', results.memoryPolynomialOrder, results.memoryOrder, results.memoryDepth);
end
