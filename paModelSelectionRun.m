function selection = paModelSelectionRun(fileName, maxPolynomialOrder, maxMemoryOrder, showPlots, maxGmpMemoryOrder)
    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2 || isempty(maxPolynomialOrder)
        maxPolynomialOrder = 25;
    end

    if nargin < 3 || isempty(maxMemoryOrder)
        maxMemoryOrder = 30;
    end

    if nargin < 4 || isempty(showPlots)
        showPlots = true;
    end

    if nargin < 5
        maxGmpMemoryOrder = [];
    end

    commonAlignmentMemoryOrder = max(maxMemoryOrder, 30);
    memorySearch = paDimensionRun(fileName, maxPolynomialOrder, maxMemoryOrder, 0.20, 0.05, showPlots);
    memorylessSearch = memorySearch.initialPolynomial;

    if ~memorySearch.coordinateConvergenceReached
        error('paModelSelectionRun:MpNoConvergence', 'The MP parameter search did not converge.');
    end

    memorylessSpec = makeSpec('Memoryless', 'memoryless', memorylessSearch.selectedPolynomialOrder, 0, 0);
    memorySpec = makeSpec('Memory Polynomial', 'memory', memorySearch.selectedPolynomialOrder, memorySearch.selectedMemoryOrder, 0);
    clear memorySearch memorylessSearch
    measuredData = paMeasuredCaptures(fileName);
    lambdaValues = logspace(-12, 2, 300);
    gmpDiagonalSearch = paGmpDiagonalProfile(measuredData, memorySpec.orders, memorySpec.memoryOrder, commonAlignmentMemoryOrder, lambdaValues, 0.05);
    gmpSpec = makeSpec('GMP', 'gmp', memorySpec.polynomialOrder, memorySpec.memoryOrder, gmpDiagonalSearch.selectedDiagonalCount);

    commonAlignmentMemoryOrder = max([commonAlignmentMemoryOrder memorySpec.memoryOrder gmpSpec.memoryOrder]);

    memorylessSpec = evaluateSpec(memorylessSpec, measuredData, lambdaValues, commonAlignmentMemoryOrder);
    memorySpec = evaluateSpec(memorySpec, measuredData, lambdaValues, commonAlignmentMemoryOrder);
    gmpSpec = evaluateSpec(gmpSpec, measuredData, lambdaValues, commonAlignmentMemoryOrder);

    model = {memorylessSpec.name; memorySpec.name; gmpSpec.name};
    polynomialOrder = [memorylessSpec.polynomialOrder; memorySpec.polynomialOrder; gmpSpec.polynomialOrder];
    memoryOrder = [memorylessSpec.memoryOrder; memorySpec.memoryOrder; gmpSpec.memoryOrder];
    diagonalCount = [memorylessSpec.diagonalCount; memorySpec.diagonalCount; gmpSpec.diagonalCount];
    coefficientCount = [memorylessSpec.coefficientCount; memorySpec.coefficientCount; gmpSpec.coefficientCount];
    bestLambda = [memorylessSpec.bestLambda; memorySpec.bestLambda; gmpSpec.bestLambda];
    validationNmseDb = [memorylessSpec.validationNmseDb; memorySpec.validationNmseDb; gmpSpec.validationNmseDb];
    summary = table(model, polynomialOrder, memoryOrder, diagonalCount, coefficientCount, bestLambda, validationNmseDb, 'VariableNames', {'Model', 'P', 'M', 'd', 'K', 'Lambda', 'NMSE_dB'});

    paCfg = paConfig(false);
    selection.fileName = fileName;
    sourceInfo = dir(fileName);
    selection.methodVersion = 11;
    selection.sourceBytes = sourceInfo.bytes;
    selection.sourceDatenum = sourceInfo.datenum;
    selection.maxPolynomialOrder = maxPolynomialOrder;
    selection.maxMemoryOrder = maxMemoryOrder;
    selection.maxGmpMemoryOrder = maxMemoryOrder;
    selection.legacyMaxGmpMemoryOrderArgument = maxGmpMemoryOrder;
    selection.polynomialThresholdDb = 0.20;
    selection.memoryThresholdDb = 0.05;
    selection.diagonalSelectionRule = 'after GMP inherits the selected MP P and M, select the smallest d within 0.05 dB of the minimum validation NMSE';
    selection.selectionRule = 'Memoryless and MP select P and M; GMP inherits the selected MP P and M and selects only d; all NMSE calculations use one 80/20 holdout';
    selection.modelSampleRate = paCfg.modelSampleRate;
    selection.commonAlignmentMemoryOrder = commonAlignmentMemoryOrder;
    selection.lambdaValues = lambdaValues;
    selection.memoryless = memorylessSpec;
    selection.memory = memorySpec;
    selection.gmp = gmpSpec;
    selection.gmpDiagonal = gmpDiagonalSearch;
    selection.models = [memorylessSpec memorySpec gmpSpec];
    selection.summary = summary;
    selection.trainingRepetitions = 1:4;
    selection.validationRepetition = 5;
    selection.lambdaSelectionRepetition = 5;
    selection.independentTestRepetition = [];
    selection.foldCount = 0;
    selection.splitMethod = '80/20 holdout';
    selection.created = datestr(now, 30);
    selectionFile = fullfile(fileparts(mfilename('fullpath')), 'paBestModels.mat');
    selection.selectionFile = selectionFile;

    save(selectionFile, 'selection', '-v7');
    disp(summary);

end

function spec = makeSpec(name, modelType, polynomialOrder, memoryOrder, diagonalCount)
    spec.name = name;
    spec.modelType = modelType;
    spec.polynomialOrder = polynomialOrder;
    spec.orders = 1:2:polynomialOrder;
    spec.memoryOrder = memoryOrder;
    spec.memoryDepth = memoryOrder + 1;
    spec.signalDelays = 0:memoryOrder;
    spec.envelopeDelays = 0:memoryOrder;
    spec.diagonalCount = diagonalCount;
end

function spec = evaluateSpec(spec, measuredData, lambdaValues, alignmentMemoryOrder)
    [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, spec.modelType, spec.orders, spec.memoryOrder, spec.diagonalCount, alignmentMemoryOrder);
    outputBlocks = measuredData.outputBlocks(sampleIndices, :);
    score = paCrossValidate(basis, outputBlocks, lambdaValues, true);
    spec.alignmentMemoryOrder = alignmentMemoryOrder;
    spec.sampleIndices = sampleIndices;
    spec.coefficientCount = size(basis, 2);
    spec.bestLambda = score.bestLambda;
    spec.validationNmseDb = score.validationNmseDb;
    spec.nmseDb = score.validationNmseDb;
    spec.testNmseDb = score.testNmseDb;
end
