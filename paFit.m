function [coefficientsBackslash, coefficientsTikhonov, fitInfo] = paFit(fileName)

measuredData = paMeasuredCaptures(fileName);
[modelSpec, modelSampleRate, alignmentMemoryOrder, selectionSource] = fitConfiguration(fileName);
lambdaValues = logspace(-30, 2, 300);

[basis, sampleIndices] = paBuildBasis(measuredData.inputBlock, 'gmp', modelSpec.orders, modelSpec.memoryOrder, modelSpec.diagonalCount, alignmentMemoryOrder);
outputBlocks = measuredData.outputBlocks(sampleIndices, :);
score = paCrossValidate(basis, outputBlocks, lambdaValues, true);

trainingOutput = mean(outputBlocks(:, 1:4), 2);
validationReference = outputBlocks(:, 5);
columnScaleFull = sqrt(mean(abs(basis).^2, 1));
columnScaleFull(columnScaleFull == 0) = 1;
normalizedBasis = bsxfun(@rdivide, basis, columnScaleFull);
coefficientsBackslashNormalized = normalizedBasis \ trainingOutput;
coefficientsBackslashActive = coefficientsBackslashNormalized ./ columnScaleFull.';
coefficientsTikhonovActive = score.finalCoefficients;
backslashPrediction = basis * coefficientsBackslashActive;

nmseBackslashValidationDb = calculateNmse(validationReference, backslashPrediction);
nmseTikhonovValidationDb = score.validationNmseDb;
nmseBackslashFullDb = calculateNmse(trainingOutput, backslashPrediction);
nmseTikhonovFullDb = calculateNmse(trainingOutput, score.finalModelOutput);

maskCfg.orders = modelSpec.orders;
maskCfg.signalDelays = modelSpec.signalDelays;
maskCfg.envelopeDelays = modelSpec.envelopeDelays;
[coefficientMask, activeIndices] = gmpDiagonalMask(maskCfg, modelSpec.diagonalCount);

coefficientsBackslash = complex(zeros(size(coefficientMask)));
coefficientsTikhonov = complex(zeros(size(coefficientMask)));
coefficientsBackslash(activeIndices) = coefficientsBackslashActive;
coefficientsTikhonov(activeIndices) = coefficientsTikhonovActive;

referenceInputRms = sqrt(mean(abs(measuredData.inputBlock).^2));

fitInfo.orders = modelSpec.orders;
fitInfo.signalDelays = modelSpec.signalDelays;
fitInfo.envelopeDelays = modelSpec.envelopeDelays;
fitInfo.diagonalCount = modelSpec.diagonalCount;
fitInfo.modelSampleRate = modelSampleRate;
fitInfo.reference_input_rms = referenceInputRms;
fitInfo.maximumInputMagnitude = max(abs(measuredData.inputBlock));
fitInfo.coefficientCount = numel(coefficientsTikhonov);
fitInfo.activeCoefficientCount = length(activeIndices);
fitInfo.coefficientMask = coefficientMask;
fitInfo.activeCoefficientIndices = activeIndices;
fitInfo.lambdaValues = lambdaValues;
fitInfo.nmseValuesDb = score.validationNmseValuesDb;
fitInfo.columnScaleTrain = columnScaleFull;
fitInfo.columnScaleFull = columnScaleFull;
fitInfo.bestLambda = score.bestLambda;
fitInfo.nmseBackslashValidationDb = nmseBackslashValidationDb;
fitInfo.nmseTikhonovValidationDb = nmseTikhonovValidationDb;
fitInfo.nmseBackslashTestDb = nmseBackslashValidationDb;
fitInfo.nmseTikhonovTestDb = nmseTikhonovValidationDb;
fitInfo.nmseBackslashFullDb = nmseBackslashFullDb;
fitInfo.nmseTikhonovFullDb = nmseTikhonovFullDb;
fitInfo.coefficientsBackslash = coefficientsBackslash;
fitInfo.coefficientsTikhonov = coefficientsTikhonov;
fitInfo.alignmentMemoryOrder = alignmentMemoryOrder;
fitInfo.trainingRepetitions = 1:4;
fitInfo.validationRepetition = 5;
fitInfo.lambdaSelectionRepetition = 5;
fitInfo.independentTestRepetition = [];
fitInfo.splitMethod = '80/20 holdout';
fitInfo.selectionSource = selectionSource;
fitInfo.selectionMethodVersion = 11;
fitInfo.nmseMetric = 'repetitions 1-4 training, repetition 5 validation';

modelFile = fullfile(fileparts(mfilename('fullpath')), 'paModel.mat');
fitInfo.modelFile = modelFile;
save(modelFile, 'coefficientsTikhonov', 'fitInfo', '-v7');

figure('Name', 'NMSE / lambda', 'NumberTitle', 'off', 'Color', 'w');
tikhonovLine = semilogx(lambdaValues, score.validationNmseValuesDb, 'LineWidth', 1.5);
hold on;
leastSquaresLine = semilogx(lambdaValues, nmseBackslashValidationDb * ones(size(lambdaValues)), 'k--', 'LineWidth', 1.2);
bestPoint = plot(score.bestLambda, nmseTikhonovValidationDb, 'ko', 'MarkerFaceColor', 'k');
bestIndex = find(lambdaValues == score.bestLambda, 1);
additionalIndices = max(1, bestIndex - 3):min(length(lambdaValues), bestIndex + 3);
additionalIndices(additionalIndices == bestIndex) = [];
plot(lambdaValues(additionalIndices), score.validationNmseValuesDb(additionalIndices), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
grid on;
xlabel('lambda');
ylabel('Validation NMSE, dB');
title(['NMSE / lambda: ' fileName], 'Interpreter', 'none');
legend([tikhonovLine leastSquaresLine bestPoint], {'Tikhonov', 'Least squares', 'Best lambda'}, 'Location', 'best');
leftIndex = max(1, bestIndex - 3);
rightIndex = min(length(lambdaValues), bestIndex + 3);
xlim([lambdaValues(leftIndex) lambdaValues(rightIndex)]);

fprintf('LS validation NMSE = %.9f dB\n', nmseBackslashValidationDb);
fprintf('Tikhonov validation NMSE = %.9f dB\n', nmseTikhonovValidationDb);
fprintf('Best lambda = %.9e\n', score.bestLambda);
fprintf('Active GMP coefficient count = %d\n', length(activeIndices));
disp('Tikhonov coefficients:');
disp(coefficientsTikhonov);

end

function [modelSpec, modelSampleRate, alignmentMemoryOrder, selectionSource] = fitConfiguration(fileName)
paCfg = paConfig(false);
modelSpec.orders = paCfg.orders;
modelSpec.signalDelays = paCfg.signalDelays;
modelSpec.envelopeDelays = paCfg.envelopeDelays;
modelSpec.memoryOrder = max([paCfg.signalDelays paCfg.envelopeDelays]);
modelSpec.diagonalCount = paCfg.diagonalCount;
modelSampleRate = paCfg.modelSampleRate;
alignmentMemoryOrder = max(30, modelSpec.memoryOrder);
selectionSource = 'paConfig';

selectionFile = fullfile(fileparts(mfilename('fullpath')), 'paBestModels.mat');
if exist(selectionFile, 'file') ~= 2
    return
end

savedSelection = load(selectionFile, 'selection');
selection = savedSelection.selection;
sourceInfo = dir(fileName);
selectionIsCurrent = isfield(selection, 'methodVersion') && selection.methodVersion == 11;
selectionIsCurrent = selectionIsCurrent && isfield(selection, 'fileName') && strcmpi(canonicalPath(selection.fileName), canonicalPath(fileName));
selectionIsCurrent = selectionIsCurrent && isfield(selection, 'sourceBytes') && selection.sourceBytes == sourceInfo.bytes;
selectionIsCurrent = selectionIsCurrent && isfield(selection, 'sourceDatenum') && selection.sourceDatenum == sourceInfo.datenum;
selectionIsCurrent = selectionIsCurrent && isfield(selection, 'gmp') && isfield(selection.gmp, 'orders');

if ~selectionIsCurrent
    return
end

modelSpec.orders = selection.gmp.orders;
modelSpec.signalDelays = selection.gmp.signalDelays;
modelSpec.envelopeDelays = selection.gmp.envelopeDelays;
modelSpec.memoryOrder = selection.gmp.memoryOrder;
modelSpec.diagonalCount = selection.gmp.diagonalCount;
modelSampleRate = selection.modelSampleRate;
alignmentMemoryOrder = selection.commonAlignmentMemoryOrder;
selectionSource = 'paBestModels methodVersion 11';
end

function pathName = canonicalPath(fileName)
[success, attributes] = fileattrib(fileName);

if success
    pathName = attributes.Name;
else
    pathName = fileName;
end
end

function nmseDb = calculateNmse(reference, estimate)
nmseDb = 10 * log10(sum(abs(reference - estimate).^2) / sum(abs(reference).^2));
end
