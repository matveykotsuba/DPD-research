function result = paCrossValidate(basis, outputBlocks, lambdaValues, calculatePrediction)
    if nargin < 4
        calculatePrediction = true;
    end

    lambdaValues = lambdaValues(:).';
    rowCount = size(basis, 1);
    trainingOutput = mean(outputBlocks(:, 1:4), 2);
    validationOutput = outputBlocks(:, 5);
    columnScale = sqrt(mean(abs(basis).^2, 1));
    columnScale(columnScale == 0) = 1;
    normalizedBasis = bsxfun(@rdivide, basis, columnScale);
    correlation = normalizedBasis' * normalizedBasis / rowCount;
    correlation = (correlation + correlation') / 2;
    trainingCrossCorrelation = normalizedBasis' * trainingOutput / rowCount;
    validationCrossCorrelation = normalizedBasis' * validationOutput / rowCount;
    validationPower = mean(abs(validationOutput).^2);

    [eigenvectors, eigenvalueMatrix] = eig(correlation);
    eigenvalues = max(real(diag(eigenvalueMatrix)), 0);
    denominator = bsxfun(@plus, eigenvalues, lambdaValues);
    coefficientCoordinates = bsxfun(@rdivide, eigenvectors' * trainingCrossCorrelation, denominator);
    normalizedCoefficients = eigenvectors * coefficientCoordinates;
    validationErrorPower = validationPower - 2 * real(sum(conj(normalizedCoefficients) .* validationCrossCorrelation, 1));
    validationErrorPower = validationErrorPower + real(sum(conj(normalizedCoefficients) .* (correlation * normalizedCoefficients), 1));
    validationNmseValuesDb = 10 * log10(max(validationErrorPower / validationPower, realmin));
    [validationNmseDb, bestLambdaIndex] = min(validationNmseValuesDb);
    bestLambda = lambdaValues(bestLambdaIndex);

    finalCoefficients = [];
    validationPrediction = [];

    if calculatePrediction
        identityMatrix = eye(size(correlation));
        finalNormalizedCoefficients = (correlation + bestLambda * identityMatrix) \ trainingCrossCorrelation;
        finalCoefficients = finalNormalizedCoefficients ./ columnScale.';
        validationPrediction = basis * finalCoefficients;
        validationNmseDb = calculateNmse(validationOutput, validationPrediction);
    end

    result.lambdaValues = lambdaValues;
    result.validationNmseValuesDb = validationNmseValuesDb;
    result.bestLambda = bestLambda;
    result.validationNmseDb = validationNmseDb;
    result.validationPrediction = validationPrediction;
    result.validationReference = validationOutput;
    result.finalCoefficients = finalCoefficients;
    result.finalModelOutput = validationPrediction;
    result.trainingOutput = trainingOutput;
    result.trainingRepetitions = 1:4;
    result.validationRepetition = 5;
    result.splitMethod = '80/20 holdout';
    result.testNmseDb = validationNmseDb;
    result.testPrediction = validationPrediction;
    result.testReference = validationOutput;
    result.foldNmseDb = [];
    result.foldEdges = [];
end

function nmseDb = calculateNmse(reference, estimate)
    nmseDb = 10 * log10(sum(abs(reference - estimate).^2) / sum(abs(reference).^2));
end
