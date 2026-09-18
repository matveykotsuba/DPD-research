function results = paParitySpectrumRun(fileName, polynomialOrder, memoryOrder, alignmentMemoryOrder)
    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    selection = [];

    if nargin < 2 || isempty(polynomialOrder) || nargin < 3 || isempty(memoryOrder)
        selection = loadSelection(fileName);
    end

    if nargin < 2 || isempty(polynomialOrder)
        polynomialOrder = selection.memory.polynomialOrder;
    end

    if nargin < 3 || isempty(memoryOrder)
        memoryOrder = selection.memory.memoryOrder;
    end

    if nargin < 4 || isempty(alignmentMemoryOrder)
        if isempty(selection)
            alignmentMemoryOrder = max(30, memoryOrder);
        else
            alignmentMemoryOrder = selection.commonAlignmentMemoryOrder;
        end
    end

    alignmentMemoryOrder = max(alignmentMemoryOrder, memoryOrder);
    oddOrders = 1:2:polynomialOrder;
    allOrders = 1:polynomialOrder;
    paCfg = paConfig(false);
    sampleRate = paCfg.modelSampleRate;
    segmentLength = 1024;
    measuredData = paMeasuredCaptures(fileName);
    lambdaValues = logspace(-12, 2, 300);

    [oddBasis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'memory', oddOrders, memoryOrder, memoryOrder, alignmentMemoryOrder);
    [allBasis, allSampleIndices] = paBuildBasis(measuredData.inputBlock, 'memory', allOrders, memoryOrder, memoryOrder, alignmentMemoryOrder);
    oddOutputBlocks = measuredData.outputBlocks(sampleIndices, :);
    allOutputBlocks = measuredData.outputBlocks(allSampleIndices, :);
    oddScore = paCrossValidate(oddBasis, oddOutputBlocks, lambdaValues, true);
    allScore = paCrossValidate(allBasis, allOutputBlocks, lambdaValues, true);
    oddOutput = oddScore.testPrediction;
    allOutput = allScore.testPrediction;
    measuredOutput = oddScore.testReference;

    [frequency, oddPsd] = welch(oddOutput, sampleRate, segmentLength);
    [~, allPsd] = welch(allOutput, sampleRate, segmentLength);
    [~, measuredPsd] = welch(measuredOutput, sampleRate, segmentLength);

    psdReference = max(measuredPsd);
    oddPsdDb = 10 * log10(oddPsd / psdReference + eps);
    allPsdDb = 10 * log10(allPsd / psdReference + eps);
    measuredPsdDb = 10 * log10(measuredPsd / psdReference + eps);

    spectrumFigure = figure('Name', 'Odd-only and all polynomial orders ', 'NumberTitle', 'off', 'Color', 'w');
    plot(frequency / 1e6, oddPsdDb, 'b', 'LineWidth', 1.5);
    hold on;
    plot(frequency / 1e6, allPsdDb, 'r', 'LineWidth', 1.5);
    plot(frequency / 1e6, measuredPsdDb, 'k', 'LineWidth', 1.8);
    grid on;
    xlabel('Frequency, MHz');
    ylabel(' PSD, dB');
    title(sprintf(' P = %d, M = %d', polynomialOrder, memoryOrder));
    legend('Odd orders ', 'All  orders', 'Real PA output', 'Location', 'best');
    xlim([-100 100]);
    ylim([-80 5]);

    selectionInput = measuredData.inputBlocks(sampleIndices, 5);
    amAmPointCount = min(24000, length(selectionInput));
    amAmPointIndices = round(linspace(1, length(selectionInput), amAmPointCount)).';
    amAmInputAmplitude = abs(selectionInput(amAmPointIndices));
    oddAmAmOutputAmplitude = abs(oddOutput(amAmPointIndices));
    allAmAmOutputAmplitude = abs(allOutput(amAmPointIndices));

    amAmFigure = figure('Name', 'Odd-only and all-orders AM-AM', 'NumberTitle', 'off', 'Color', 'w');
    hold on;
    oddAmAmHandle = plot(amAmInputAmplitude, oddAmAmOutputAmplitude, 'b.', 'LineStyle', 'none', 'MarkerSize', 4);
    allAmAmHandle = plot(amAmInputAmplitude, allAmAmOutputAmplitude, 'r.', 'LineStyle', 'none', 'MarkerSize', 4);
    maximumAmplitude = max([amAmInputAmplitude; oddAmAmOutputAmplitude; allAmAmOutputAmplitude]);
    plot([0 maximumAmplitude], [0 maximumAmplitude], 'k--', 'LineWidth', 1.0);
    grid on;
    axis equal;
    xlim([0 maximumAmplitude]);
    ylim([0 maximumAmplitude]);
    xlabel('input amplitude');
    ylabel('output amplitude');
    legend([oddAmAmHandle allAmAmHandle], {'Odd orders', 'All orders'}, 'Location', 'best');

    measuredPaInputAmplitude = amAmInputAmplitude;
    measuredPaOutputAmplitude = abs(measuredOutput(amAmPointIndices));
    measuredPaMaximum = max([measuredPaInputAmplitude; measuredPaOutputAmplitude]);

    measuredPaAmAmFigure = figure('Name', 'Real PA AM-AM', 'NumberTitle', 'off', 'Color', 'w');
    plot(measuredPaInputAmplitude, measuredPaOutputAmplitude, 'b.', 'LineStyle', 'none', 'MarkerSize', 4);
    hold on;
    plot([0 measuredPaMaximum], [0 measuredPaMaximum], 'k--', 'LineWidth', 1.0);
    grid on;
    axis equal;
    xlim([0 measuredPaMaximum]);
    ylim([0 measuredPaMaximum]);
    xlabel('input amplitude');
    ylabel('output amplitude');
    title('Real PA AM/AM');

    results.fileName = fileName;
    results.model = 'Memory Polynomial';
    results.polynomialOrder = polynomialOrder;
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryOrder + 1;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.sampleIndices = sampleIndices;
    results.oddSampleIndices = sampleIndices;
    results.allSampleIndices = allSampleIndices;
    results.oddOrders = oddOrders;
    results.allOrders = allOrders;
    results.oddCoefficientCount = size(oddBasis, 2);
    results.allCoefficientCount = size(allBasis, 2);
    results.oddLambda = oddScore.bestLambda;
    results.allLambda = allScore.bestLambda;
    results.oddValidationNmseDb = oddScore.validationNmseDb;
    results.allValidationNmseDb = allScore.validationNmseDb;
    results.oddTestNmseDb = oddScore.testNmseDb;
    results.allTestNmseDb = allScore.testNmseDb;
    results.oddScore = oddScore;
    results.allScore = allScore;
    results.sampleRate = sampleRate;
    results.segmentLength = segmentLength;
    results.frequency = frequency;
    results.oddPsd = oddPsd;
    results.allPsd = allPsd;
    results.measuredPsd = measuredPsd;
    results.oddPsdDb = oddPsdDb;
    results.allPsdDb = allPsdDb;
    results.measuredPsdDb = measuredPsdDb;
    results.spectrumFigure = spectrumFigure;
    results.amAmInputAmplitude = amAmInputAmplitude;
    results.oddAmAmOutputAmplitude = oddAmAmOutputAmplitude;
    results.allAmAmOutputAmplitude = allAmAmOutputAmplitude;
    results.amAmFigure = amAmFigure;
    results.measuredPaInputAmplitude = measuredPaInputAmplitude;
    results.measuredPaOutputAmplitude = measuredPaOutputAmplitude;
    results.measuredPaAmAmFigure = measuredPaAmAmFigure;
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
