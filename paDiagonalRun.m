function results = paDiagonalRun(fileName, orders, memoryOrder, alignmentMemoryOrder, showPlots)
    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    paCfg = paConfig(false);
    finalSelection = loadFinalSelection(fileName);

    if nargin < 2 || isempty(orders)
        if isempty(finalSelection)
            orders = paCfg.orders;
        else
            orders = finalSelection.gmp.orders;
        end
    end

    if nargin < 3 || isempty(memoryOrder)
        if isempty(finalSelection)
            memoryOrder = max([paCfg.signalDelays paCfg.envelopeDelays]);
        else
            memoryOrder = finalSelection.gmp.memoryOrder;
        end
    end

    if nargin < 4 || isempty(alignmentMemoryOrder)
        if isempty(finalSelection)
            alignmentMemoryOrder = 30;
        else
            alignmentMemoryOrder = finalSelection.commonAlignmentMemoryOrder;
        end
    end

    if nargin < 5
        showPlots = true;
    end

    orders = orders(:).';
    alignmentMemoryOrder = max(alignmentMemoryOrder, memoryOrder);
    measuredData = paMeasuredCaptures(fileName);

    if isempty(finalSelection)
        lambdaValues = logspace(-12, 2, 300);
        sampleRate = paCfg.modelSampleRate;
        parameterSource = 'paConfig';
    else
        lambdaValues = finalSelection.lambdaValues;
        sampleRate = finalSelection.modelSampleRate;
        parameterSource = 'paBestModels.mat';
    end

    testedDiagonalCount = (0:memoryOrder).';
    modelCount = length(testedDiagonalCount);
    coefficientCountAll = zeros(modelCount, 1);
    bestLambdaAll = zeros(modelCount, 1);
    crossValidationNmseDbAll = zeros(modelCount, 1);
    fullCoefficientCount = memoryOrder + 1 + (length(orders) - 1) * (memoryOrder + 1)^2;
    selectedColumnsAll = false(fullCoefficientCount, modelCount);
    coefficientVectorsAll = complex(zeros(fullCoefficientCount, modelCount));
    sampleIndices = (alignmentMemoryOrder + 1:measuredData.blockLength).';
    modelOutputsAll = complex(zeros(length(sampleIndices), modelCount));
    scoreAll = cell(modelCount, 1);

    for diagonalIndex = 1:modelCount
        currentDiagonalCount = testedDiagonalCount(diagonalIndex);
        [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'gmp', orders, memoryOrder, currentDiagonalCount, alignmentMemoryOrder);
        outputBlocks = measuredData.outputBlocks(sampleIndices, :);
        score = paCrossValidate(basis, outputBlocks, lambdaValues, true);

        selectedColumns = gmpColumnMask(orders, memoryOrder, currentDiagonalCount);
        coefficientCountAll(diagonalIndex) = size(basis, 2);
        bestLambdaAll(diagonalIndex) = score.bestLambda;
        crossValidationNmseDbAll(diagonalIndex) = score.validationNmseDb;
        selectedColumnsAll(:, diagonalIndex) = selectedColumns(:);
        coefficientVectorsAll(selectedColumns(:), diagonalIndex) = score.finalCoefficients(:);
        modelOutputsAll(:, diagonalIndex) = score.finalModelOutput;
        scoreAll{diagonalIndex} = score;
    end

    diagonalToleranceDb = 0.05;
    [bestNmseDb, bestIndex] = min(crossValidationNmseDbAll);
    bestDiagonalCount = testedDiagonalCount(bestIndex);
    selectedIndex = find(crossValidationNmseDbAll <= bestNmseDb + diagonalToleranceDb, 1, 'first');
    selectedDiagonalCount = testedDiagonalCount(selectedIndex);
    selectedNmseDb = crossValidationNmseDbAll(selectedIndex);
    selectedScore = scoreAll{selectedIndex};
    nmseImprovementDbAll = [NaN; crossValidationNmseDbAll(1:end-1) - crossValidationNmseDbAll(2:end)];
    fullSummary = table(testedDiagonalCount, coefficientCountAll, bestLambdaAll, crossValidationNmseDbAll, nmseImprovementDbAll, 'VariableNames', {'SideDiagonalsPerSide', 'CoefficientCount', 'BestLambda', 'NMSE_dB', 'NMSEImprovement_dB'});
    disp(fullSummary);

    sideDiagonalCount = testedDiagonalCount(2:end);
    coefficientCount = coefficientCountAll(2:end);
    bestLambda = bestLambdaAll(2:end);
    validationNmseDb = crossValidationNmseDbAll(2:end);
    nmseImprovementDb = nmseImprovementDbAll(2:end);
    summary = fullSummary(2:end, :);
    diagonalSelectedColumns = selectedColumnsAll(:, 2:end);
    diagonalCoefficientVectors = coefficientVectorsAll(:, 2:end);
    diagonalModelOutputs = modelOutputsAll(:, 2:end);
    alignedInput = measuredData.inputBlock(sampleIndices);
    measuredOutput = mean(measuredData.outputBlocks(sampleIndices, 1:4), 2);
    comparisonSegmentLength = 2048;
    [frequency, measuredPsd, welchInfo] = welch(measuredOutput, sampleRate, comparisonSegmentLength);
    diagonalPsd = zeros(length(measuredPsd), length(sideDiagonalCount));

    for modelIndex = 1:length(sideDiagonalCount)
        [~, diagonalPsd(:, modelIndex)] = welch(diagonalModelOutputs(:, modelIndex), sampleRate, comparisonSegmentLength);
    end

    psdReference = max(measuredPsd);
    measuredPsdDb = 10 * log10(measuredPsd / psdReference + eps);
    diagonalPsdDb = 10 * log10(diagonalPsd / psdReference + eps);
    amAmPointCount = min(24000, length(alignedInput));
    amAmPointIndices = round(linspace(1, length(alignedInput), amAmPointCount)).';
    amAmInputAmplitude = abs(alignedInput(amAmPointIndices));
    diagonalAmAm = abs(diagonalModelOutputs(amAmPointIndices, :));
    nmseFigure = [];
    spectrumFigure = [];
    amAmFigure = [];

    if showPlots
        nmseFigure = plotDiagonalNmse(testedDiagonalCount, crossValidationNmseDbAll, selectedDiagonalCount, orders(end), memoryOrder);

        if ~isempty(sideDiagonalCount)
            spectrumFigure = plotDiagonalSpectrumComparison(frequency, measuredPsdDb, diagonalPsdDb, sideDiagonalCount, validationNmseDb);
            amAmFigure = plotDiagonalAmAmComparison(amAmInputAmplitude, diagonalAmAm, sideDiagonalCount);
        end
    end

    fprintf('Memory Polynomial baseline NMSE: %.9f dB, %d coefficients.\n', crossValidationNmseDbAll(1), coefficientCountAll(1));
    fprintf('Global minimum: d = %d, NMSE = %.9f dB.\n', bestDiagonalCount, bestNmseDb);
    fprintf('Selected diagonal count: d = %d, NMSE = %.9f dB, %d coefficients.\n', selectedDiagonalCount, selectedNmseDb, coefficientCountAll(selectedIndex));

    finalSelectionMatch = false;
    modelCompareValidationDifferenceDb = NaN;

    if ~isempty(finalSelection) && isequal(orders, finalSelection.gmp.orders) && memoryOrder == finalSelection.gmp.memoryOrder
        finalSelectionMatch = selectedDiagonalCount == finalSelection.gmp.diagonalCount;
        finalIndex = find(testedDiagonalCount == finalSelection.gmp.diagonalCount, 1);

        if ~isempty(finalIndex)
            modelCompareValidationDifferenceDb = crossValidationNmseDbAll(finalIndex) - finalSelection.gmp.validationNmseDb;
            fprintf('Model Compare match at d = %d: NMSE difference = %.3e dB.\n', finalSelection.gmp.diagonalCount, modelCompareValidationDifferenceDb);
        end
    end

    results.fileName = fileName;
    results.model = 'Generalized Memory Polynomial';
    results.sampleRate = sampleRate;
    results.parameterSource = parameterSource;
    results.firstModelSample = sampleIndices(1);
    results.sampleCount = length(sampleIndices);
    results.orders = orders;
    results.memoryOrder = memoryOrder;
    results.memoryDepth = memoryOrder + 1;
    results.signalDelays = 0:memoryOrder;
    results.envelopeDelays = 0:memoryOrder;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.diagonalRule = '|m - q| <= d';
    results.nmseMetric = '80/20 holdout: repetitions 1-4 train the model; repetition 5 selects lambda and reports NMSE';
    results.diagonalSelectionToleranceDb = diagonalToleranceDb;
    results.diagonalSelectionRule = 'smallest d within 0.05 dB of the minimum validation NMSE';
    results.psdFitScope = 'repetitions 1-4 fit with repetition-5-selected lambda';
    results.modelOutputFitScope = results.psdFitScope;
    results.lambdaValues = lambdaValues;
    results.testedDiagonalCount = testedDiagonalCount;
    results.allCoefficientCount = coefficientCountAll;
    results.allBestLambda = bestLambdaAll;
    results.allCrossValidationNmseDb = crossValidationNmseDbAll;
    results.allNmseImprovementDb = nmseImprovementDbAll;
    results.fullSummary = fullSummary;
    results.bestDiagonalCount = bestDiagonalCount;
    results.selectedDiagonalCount = selectedDiagonalCount;
    results.bestNmseDb = bestNmseDb;
    results.selectedNmseDb = selectedNmseDb;
    results.bestTestNmseDb = bestNmseDb;
    results.selectedTestNmseDb = selectedNmseDb;
    results.selectedScore = selectedScore;
    results.bestScore = scoreAll{bestIndex};
    results.finalSelectionMatch = finalSelectionMatch;
    results.modelCompareValidationDifferenceDb = modelCompareValidationDifferenceDb;
    results.memoryPolynomialCoefficientCount = coefficientCountAll(1);
    results.memoryPolynomialNmseDb = crossValidationNmseDbAll(1);
    results.memoryPolynomialValidationNmseDb = crossValidationNmseDbAll(1);
    results.sideDiagonalCount = sideDiagonalCount;
    results.coefficientCount = coefficientCount;
    results.bestLambda = bestLambda;
    results.validationNmseDb = validationNmseDb;
    results.crossValidationNmseDb = validationNmseDb;
    results.testNmseDb = validationNmseDb;
    results.nmseImprovementDb = nmseImprovementDb;
    results.summary = summary;
    results.alignedInput = alignedInput;
    results.measuredOutput = measuredOutput;
    results.diagonalSelectedColumns = diagonalSelectedColumns;
    results.diagonalCoefficientVectors = diagonalCoefficientVectors;
    results.diagonalModelOutputs = diagonalModelOutputs;
    results.frequency = frequency;
    results.measuredPsd = measuredPsd;
    results.diagonalPsd = diagonalPsd;
    results.psdReference = psdReference;
    results.measuredPsdDb = measuredPsdDb;
    results.diagonalPsdDb = diagonalPsdDb;
    results.comparisonSegmentLength = comparisonSegmentLength;
    results.welchInfo = welchInfo;
    results.amAmInputAmplitude = amAmInputAmplitude;
    results.diagonalAmAm = diagonalAmAm;
    results.amAmSideDiagonalCount = sideDiagonalCount;
    results.nmseFigure = nmseFigure;
    results.spectrumFigure = spectrumFigure;
    results.amAmFigure = amAmFigure;
end

function selectedColumns = gmpColumnMask(orders, memoryOrder, diagonalCount)
    delayCount = memoryOrder + 1;
    coefficientCount = delayCount + (length(orders) - 1) * delayCount^2;
    selectedColumns = false(coefficientCount, 1);
    selectedColumns(1:delayCount) = true;
    columnIndex = delayCount + 1;

    for orderIndex = 2:length(orders)
        for envelopeDelay = 0:memoryOrder
            for signalDelay = 0:memoryOrder
                selectedColumns(columnIndex) = abs(signalDelay - envelopeDelay) <= diagonalCount;
                columnIndex = columnIndex + 1;
            end
        end
    end
end

function nmseFigure = plotDiagonalNmse(diagonalCount, validationNmseDb, selectedDiagonalCount, polynomialOrder, memoryOrder)
    nmseFigure = figure('Name', 'NMSE / GMP side diagonals', 'NumberTitle', 'off', 'Color', 'w');
    plot(diagonalCount, validationNmseDb, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 7);
    hold on;

    selectedIndex = find(diagonalCount == selectedDiagonalCount, 1);
    plot(selectedDiagonalCount, validationNmseDb(selectedIndex), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);

    grid on;
    xticks(diagonalCount);
    xlabel('diagonals per side');
    ylabel('NMSE, dB');
    title(sprintf('diagonals , P = %d, M = %d', polynomialOrder, memoryOrder));
end

function selection = loadFinalSelection(fileName)
    selection = [];
    selectionFile = fullfile(fileparts(mfilename('fullpath')), 'paBestModels.mat');

    if exist(selectionFile, 'file') ~= 2
        return
    end

    savedSelection = load(selectionFile, 'selection');
    candidate = savedSelection.selection;
    sourceInfo = dir(fileName);
    isCurrent = isfield(candidate, 'methodVersion') && candidate.methodVersion == 11;
    isCurrent = isCurrent && strcmpi(canonicalPath(candidate.fileName), canonicalPath(fileName));
    isCurrent = isCurrent && isfield(candidate, 'sourceBytes') && candidate.sourceBytes == sourceInfo.bytes;
    isCurrent = isCurrent && isfield(candidate, 'sourceDatenum') && candidate.sourceDatenum == sourceInfo.datenum;

    if isCurrent
        selection = candidate;
    end
end

function pathName = canonicalPath(fileName)
    [success, attributes] = fileattrib(fileName);

    if success
        pathName = attributes.Name;
    else
        pathName = fileName;
    end
end

function spectrumFigure = plotDiagonalSpectrumComparison(frequency, measuredPsdDb, diagonalPsdDb, sideDiagonalCount, validationNmseDb)
    spectrumFigure = figure('Name', 'Measured PA / GMP diagonal spectra', 'NumberTitle', 'off', 'Color', 'w');
    hold on;
    modelCount = length(sideDiagonalCount);
    colors = lines(modelCount);
    lineStyles = {'--', '-.', ':', '-'};
    legendEntries = cell(modelCount + 1, 1);
    legendEntries{1} = 'Real PA output';
    lineHandles = gobjects(modelCount + 1, 1);

    for modelIndex = 1:modelCount
        styleIndex = mod(modelIndex - 1, length(lineStyles)) + 1;
        lineHandles(modelIndex + 1) = plot(frequency / 1e6, diagonalPsdDb(:, modelIndex), 'Color', colors(modelIndex, :), 'LineStyle', lineStyles{styleIndex}, 'LineWidth', 1.3);
        legendEntries{modelIndex + 1} = sprintf('GMP, d = %d (CV NMSE %.2f dB)', sideDiagonalCount(modelIndex), validationNmseDb(modelIndex));
    end

    lineHandles(1) = plot(frequency / 1e6, measuredPsdDb, 'k', 'LineWidth', 1.8);
    grid on;
    xlabel('Frequency, MHz');
    ylabel(' PSD, dB');
    legend(lineHandles, legendEntries, 'Location', 'best');
    xlim([-100 100]);
    ylim([-80 5]);
end

function amAmFigure = plotDiagonalAmAmComparison(inputAmplitude, diagonalAmAm, sideDiagonalCount)
    amAmFigure = figure('Name', 'GMP diagonal AM-AM', 'NumberTitle', 'off', 'Color', 'w');
    hold on;
    modelCount = length(sideDiagonalCount);
    colors = lines(modelCount);
    lineHandles = gobjects(modelCount, 1);
    legendEntries = cell(modelCount, 1);

    for modelIndex = 1:modelCount
        lineHandles(modelIndex) = plot(inputAmplitude, diagonalAmAm(:, modelIndex), 'Color', colors(modelIndex, :), 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 4);
        legendEntries{modelIndex} = sprintf('D%d', sideDiagonalCount(modelIndex));
    end

    maximumAmplitude = max([inputAmplitude; diagonalAmAm(:)]);
    plot([0 maximumAmplitude], [0 maximumAmplitude], 'k--', 'LineWidth', 1.0);
    grid on;
    axis equal;
    xlim([0 maximumAmplitude]);
    ylim([0 maximumAmplitude]);
    xlabel('input amplitude');
    ylabel('output amplitude');
    legend(lineHandles, legendEntries, 'Location', 'best');
end
