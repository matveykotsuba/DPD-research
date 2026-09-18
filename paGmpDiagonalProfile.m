function results = paGmpDiagonalProfile(measuredData, orders, memoryOrder, alignmentMemoryOrder, lambdaValues, selectionToleranceDb)
    if nargin < 6
        selectionToleranceDb = 0.05;
    end

    orders = orders(:).';
    diagonalCount = (0:memoryOrder).';
    coefficientCount = zeros(size(diagonalCount));
    bestLambda = zeros(size(diagonalCount));
    crossValidationNmseDb = zeros(size(diagonalCount));

    for diagonalIndex = 1:length(diagonalCount)
        currentDiagonalCount = diagonalCount(diagonalIndex);
        [basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'gmp', orders, memoryOrder, currentDiagonalCount, alignmentMemoryOrder);
        score = paCrossValidate(basis, measuredData.outputBlocks(sampleIndices, :), lambdaValues, false);
        coefficientCount(diagonalIndex) = size(basis, 2);
        bestLambda(diagonalIndex) = score.bestLambda;
        crossValidationNmseDb(diagonalIndex) = score.validationNmseDb;
    end

    [bestNmseDb, bestIndex] = min(crossValidationNmseDb);
    selectedIndex = find(crossValidationNmseDb <= bestNmseDb + selectionToleranceDb, 1, 'first');
    selectedNmseDb = crossValidationNmseDb(selectedIndex);
    selectedDiagonalCount = diagonalCount(selectedIndex);
    summary = table(diagonalCount, coefficientCount, bestLambda, crossValidationNmseDb, 'VariableNames', {'DiagonalCount', 'CoefficientCount', 'BestLambda', 'CrossValidationNMSE_dB'});

    results.orders = orders;
    results.memoryOrder = memoryOrder;
    results.alignmentMemoryOrder = alignmentMemoryOrder;
    results.diagonalCount = diagonalCount;
    results.coefficientCount = coefficientCount;
    results.bestLambda = bestLambda;
    results.crossValidationNmseDb = crossValidationNmseDb;
    results.validationNmseDb = crossValidationNmseDb;
    results.summary = summary;
    results.bestIndex = bestIndex;
    results.bestDiagonalCount = diagonalCount(bestIndex);
    results.bestNmseDb = bestNmseDb;
    results.selectedIndex = selectedIndex;
    results.selectedDiagonalCount = selectedDiagonalCount;
    results.selectedCoefficientCount = coefficientCount(selectedIndex);
    results.selectedBestLambda = bestLambda(selectedIndex);
    results.selectedNmseDb = selectedNmseDb;
    results.selectionToleranceDb = selectionToleranceDb;
    results.selectionRule = 'smallest d within tolerance of the minimum validation NMSE';
end
