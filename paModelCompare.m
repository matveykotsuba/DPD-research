function results = paModelCompare(fileName, selection)
    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2 || isempty(selection)
        selectionFile = fullfile(fileparts(mfilename('fullpath')), 'paBestModels.mat');
        expectedPolynomialOrder = 25;
        expectedMemoryOrder = 30;
        sourceInfo = dir(fileName);

        if exist(selectionFile, 'file') == 2
            savedSelection = load(selectionFile, 'selection');
            selection = savedSelection.selection;

            cacheIsCurrent = isfield(selection, 'methodVersion') && selection.methodVersion == 11;
            cacheIsCurrent = cacheIsCurrent && strcmpi(canonicalPath(selection.fileName), canonicalPath(fileName));
            cacheIsCurrent = cacheIsCurrent && isfield(selection, 'sourceBytes') && selection.sourceBytes == sourceInfo.bytes;
            cacheIsCurrent = cacheIsCurrent && isfield(selection, 'sourceDatenum') && selection.sourceDatenum == sourceInfo.datenum;
            cacheIsCurrent = cacheIsCurrent && isfield(selection, 'maxPolynomialOrder') && selection.maxPolynomialOrder == expectedPolynomialOrder;
            cacheIsCurrent = cacheIsCurrent && isfield(selection, 'maxMemoryOrder') && selection.maxMemoryOrder == expectedMemoryOrder;

            if ~cacheIsCurrent
                selection = paModelSelectionRun(fileName, expectedPolynomialOrder, expectedMemoryOrder, false);
            end
        else
            selection = paModelSelectionRun(fileName, expectedPolynomialOrder, expectedMemoryOrder, false);
        end
    end

    measuredData = paMeasuredCaptures(fileName);
    alignmentMemoryOrder = selection.commonAlignmentMemoryOrder;
    memorylessResult = evaluateModel(selection.memoryless, measuredData, alignmentMemoryOrder);
    memoryResult = evaluateModel(selection.memory, measuredData, alignmentMemoryOrder);
    gmpResult = evaluateModel(selection.gmp, measuredData, alignmentMemoryOrder);
    sampleIndices = memorylessResult.sampleIndices;
    validationReference = measuredData.outputBlocks(sampleIndices, 5);
    sampleRate = selection.modelSampleRate;
    comparisonSegmentLength = 2048;

    [frequency, measuredPsd] = welch(validationReference, sampleRate, comparisonSegmentLength);
    [memorylessFrequency, memorylessPsd] = welch(memorylessResult.score.validationPrediction, sampleRate, comparisonSegmentLength);
    [memoryFrequency, memoryPsd] = welch(memoryResult.score.validationPrediction, sampleRate, comparisonSegmentLength);
    [gmpFrequency, gmpPsd] = welch(gmpResult.score.validationPrediction, sampleRate, comparisonSegmentLength);

    psdReference = max(measuredPsd);
    measuredPsdDb = 10 * log10(measuredPsd / psdReference + eps);
    memorylessPsdDb = 10 * log10(memorylessPsd / psdReference + eps);
    memoryPsdDb = 10 * log10(memoryPsd / psdReference + eps);
    gmpPsdDb = 10 * log10(gmpPsd / psdReference + eps);

    model = {selection.memoryless.name; selection.memory.name; selection.gmp.name};
    polynomialOrder = [selection.memoryless.polynomialOrder; selection.memory.polynomialOrder; selection.gmp.polynomialOrder];
    memoryOrder = [selection.memoryless.memoryOrder; selection.memory.memoryOrder; selection.gmp.memoryOrder];
    diagonalCount = [selection.memoryless.diagonalCount; selection.memory.diagonalCount; selection.gmp.diagonalCount];
    coefficientCount = [memorylessResult.coefficientCount; memoryResult.coefficientCount; gmpResult.coefficientCount];
    bestLambda = [memorylessResult.score.bestLambda; memoryResult.score.bestLambda; gmpResult.score.bestLambda];
    validationNmseDb = [memorylessResult.score.validationNmseDb; memoryResult.score.validationNmseDb; gmpResult.score.validationNmseDb];
    summary = table(model, polynomialOrder, memoryOrder, diagonalCount, coefficientCount, bestLambda, validationNmseDb, 'VariableNames', {'Model', 'P', 'M', 'd', 'K', 'Lambda', 'NMSE_dB'});
    disp(summary);

    spectrumFigure = figure('Name', 'Memoryless Memory Polynomial GMP', 'NumberTitle', 'off', 'Color', 'w');
    plot(frequency / 1e6, measuredPsdDb, 'k', 'LineWidth', 1.6);
    hold on;
    plot(memorylessFrequency / 1e6, memorylessPsdDb, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
    plot(memoryFrequency / 1e6, memoryPsdDb, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
    plot(gmpFrequency / 1e6, gmpPsdDb, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);
    grid on;
    xlabel('Frequency, MHz');
    ylabel(' PSD, dB');
    titleText = sprintf('Validation NMSE: Memoryless = %.2f dB, MP = %.2f dB, GMP = %.2f dB', validationNmseDb(1), validationNmseDb(2), validationNmseDb(3));
    title(titleText);
    legend('122M dpd 775 output', 'Memoryless model', 'Memory Polynomial model', 'GMP model', 'Location', 'best');
    xlim([-100 100]);
    ylim([-80 5]);

    results.fileName = fileName;
    results.selection = selection;
    results.sampleRate = sampleRate;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.sampleIndices = sampleIndices;
    results.sampleCount = length(validationReference);
    results.summary = summary;
    results.validationReference = validationReference;
    results.testReference = validationReference;
    results.measuredOutput = validationReference;
    results.memoryless = memorylessResult;
    results.memory = memoryResult;
    results.gmp = gmpResult;
    results.memorylessSpec = selection.memoryless;
    results.memorySpec = selection.memory;
    results.gmpSpec = selection.gmp;
    results.memorylessOutput = memorylessResult.score.validationPrediction;
    results.memoryOutput = memoryResult.score.validationPrediction;
    results.gmpOutput = gmpResult.score.validationPrediction;
    results.frequency = frequency;
    results.measuredPsd = measuredPsd;
    results.memorylessPsd = memorylessPsd;
    results.memoryPsd = memoryPsd;
    results.gmpPsd = gmpPsd;
    results.psdReference = psdReference;
    results.measuredPsdDb = measuredPsdDb;
    results.memorylessPsdDb = memorylessPsdDb;
    results.memoryPsdDb = memoryPsdDb;
    results.gmpPsdDb = gmpPsdDb;
    results.comparisonSegmentLength = comparisonSegmentLength;
    results.spectrumFigure = spectrumFigure;
end

function modelResult = evaluateModel(spec, measuredData, alignmentMemoryOrder)
    [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, spec.modelType, spec.orders, spec.memoryOrder, spec.diagonalCount, alignmentMemoryOrder);
    outputBlocks = measuredData.outputBlocks(sampleIndices, :);
    score = paCrossValidate(basis, outputBlocks, spec.bestLambda, true);
    modelResult.spec = spec;
    modelResult.sampleIndices = sampleIndices;
    modelResult.coefficientCount = size(basis, 2);
    modelResult.score = score;
end

function pathName = canonicalPath(fileName)
    [success, attributes] = fileattrib(fileName);

    if success
        pathName = attributes.Name;
    else
        pathName = fileName;
    end
end
