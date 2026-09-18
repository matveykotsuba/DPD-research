function results = paParityRun(fileName, maxPolynomialOrder, improvementThresholdDb, memoryOrder, alignmentMemoryOrder)

    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    selection = [];

    if nargin < 2 || isempty(maxPolynomialOrder) || nargin < 4 || isempty(memoryOrder)
        selection = loadSelection(fileName);
    end

    if nargin < 2 || isempty(maxPolynomialOrder)
        maxPolynomialOrder = 20;
    end

    if nargin < 3 || isempty(improvementThresholdDb)
        improvementThresholdDb = 0.20;
    end

    if nargin < 4 || isempty(memoryOrder)
        memoryOrder = selection.memory.memoryOrder;
    end

    if nargin < 5 || isempty(alignmentMemoryOrder)
        alignmentMemoryOrder = max(30, memoryOrder);
    end

    alignmentMemoryOrder = max(alignmentMemoryOrder, memoryOrder);

    oddCandidateOrders = 1:2:maxPolynomialOrder;
    allCandidateOrders = 1:maxPolynomialOrder;

    oddResults = paPolynomialRun(fileName, maxPolynomialOrder, improvementThresholdDb, memoryOrder, alignmentMemoryOrder, false, oddCandidateOrders, true, false);
    allResults = paPolynomialRun(fileName, maxPolynomialOrder, improvementThresholdDb, memoryOrder, alignmentMemoryOrder, false, allCandidateOrders, true, false);

    commonMaximumOrder = oddResults.polynomialOrder;
    [isPresent, allIndices] = ismember( commonMaximumOrder, allResults.polynomialOrder);
    oddIndices = find(isPresent);
    commonMaximumOrder = commonMaximumOrder(isPresent);
    allIndices = allIndices(isPresent);

    oddCoefficientCount = oddResults.coefficientCount(oddIndices);
    allCoefficientCount = allResults.coefficientCount(allIndices);
    oddNmseDb = oddResults.selectionNmseDb(oddIndices);
    allNmseDb = allResults.selectionNmseDb(allIndices);
    allOrdersGainDb = oddNmseDb - allNmseDb;

    comparison = table(commonMaximumOrder, oddCoefficientCount, allCoefficientCount, oddNmseDb, allNmseDb, allOrdersGainDb, 'VariableNames', {'MaximumOrder', 'OddCoefficientCount', 'AllCoefficientCount', 'OddNMSE_dB', 'AllNMSE_dB', 'AllOrdersGain_dB'});

    disp(comparison);

    figure('Name', sprintf('Odd-only and all polynomial orders, M = %d', memoryOrder), 'NumberTitle', 'off', 'Color', 'w');
    plot(oddResults.polynomialOrder, oddResults.selectionNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
    plot(allResults.polynomialOrder, allResults.selectionNmseDb, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xticks(1:maxPolynomialOrder);
    xlabel(' polynomial order');
    ylabel('NMSE, dB');
    legend( 'Odd orders only', 'All integer orders', 'Location', 'best');

    results.fileName = fileName;
    results.memoryOrder = memoryOrder;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.oddCandidateOrders = oddCandidateOrders;
    results.allCandidateOrders = allCandidateOrders;
    results.oddOnly = oddResults;
    results.allOrders = allResults;
    results.comparison = comparison;
end

function selection = loadSelection(fileName)
    selectionFile = fullfile(fileparts(mfilename('fullpath')), 'paBestModels.mat');
    cacheIsCurrent = false;

    if exist(selectionFile, 'file') == 2
        savedSelection = load(selectionFile, 'selection');
        selection = savedSelection.selection;
        sourceInfo = dir(fileName);
        cacheIsCurrent = isfield(selection, 'methodVersion') && selection.methodVersion == 11;
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'memory');
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'fileName') && strcmpi(selection.fileName, fileName);
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'sourceBytes') && selection.sourceBytes == sourceInfo.bytes;
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'sourceDatenum') && selection.sourceDatenum == sourceInfo.datenum;
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'maxPolynomialOrder') && selection.maxPolynomialOrder == 25;
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'maxMemoryOrder') && selection.maxMemoryOrder == 30;
        cacheIsCurrent = cacheIsCurrent && isfield(selection, 'maxGmpMemoryOrder') && selection.maxGmpMemoryOrder == 8;
    end

    if ~cacheIsCurrent
        selection = paModelSelectionRun(fileName, 25, 30, false, 8);
    end
end
