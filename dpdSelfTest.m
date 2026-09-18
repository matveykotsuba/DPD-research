function results = dpdSelfTest()

rng(23, 'twister');

dpdCfg = dpdConfig([], false);
expectedCounts = [20 44 62 74 80];
actualCounts = zeros(size(expectedCounts));

for diagonalCount = 0:4
    mask = gmpDiagonalMask(dpdCfg, diagonalCount);
    actualCounts(diagonalCount + 1) = nnz(mask);
end

assert(isequal(actualCounts, expectedCounts), 'Unexpected GMP diagonal-mask dimensions.');

sampleCount = 4096;
inputSignal = complex(randn(sampleCount, 1), randn(sampleCount, 1));
inputSignal = inputSignal / sqrt(mean(abs(inputSignal).^2));

offCfg = dpdCfg;
offCfg.enabled = false;
bypassOutput = dpdCore(inputSignal, offCfg);
assert(isequal(bypassOutput, inputSignal), 'Disabled DPD is not an exact bypass.');

modelCfg = dpdCfg;
modelCfg.enabled = true;
modelCfg.diagonalCount = dpdCfg.diagonalCount;
modelCfg.coefficients = complex(zeros(size(modelCfg.coefficients)));
modelCfg.coefficients(1, 1) = 0.92 + 0.03i;
modelCfg.coefficients(2, 1) = 0.04 - 0.01i;
activeMask = gmpDiagonalMask(modelCfg, modelCfg.diagonalCount);
nonlinearIndices = find(activeMask);
nonlinearIndices = nonlinearIndices(nonlinearIndices > length(modelCfg.signalDelays));
modelCfg.coefficients(nonlinearIndices(1:8)) = 0.002 * complex(randn(8, 1), randn(8, 1));

fullOutput = gmpCore(inputSignal, modelCfg);
blockOutput = complex(zeros(size(inputSignal)));
blockState = [];
firstSample = 1;
blockLengths = [257 991 64 1536 1248];

for blockIndex = 1:length(blockLengths)
    lastSample = firstSample + blockLengths(blockIndex) - 1;
    [blockOutput(firstSample:lastSample), blockState] = gmpCoreBlock(inputSignal(firstSample:lastSample), modelCfg, blockState);
    firstSample = lastSample + 1;
end

blockRelativeError = norm(blockOutput - fullOutput) / max(norm(fullOutput), realmin);
assert(blockRelativeError < 1e-12, 'Block GMP does not match full-record GMP.');

alignmentCfg = dpdCfg;
alignmentCfg.feedbackAlignmentMaximumDelay = 16;
alignmentCfg.feedbackRemoveDc = false;
knownDelay = 7;
knownGain = 0.8 * exp(1i * 0.35);
delayedOutput = [complex(zeros(knownDelay, 1)); knownGain * inputSignal(1:end-knownDelay)];
[alignedInput, normalizedOutput, alignmentInfo] = dpdFeedbackAlign(inputSignal, delayedOutput, alignmentCfg);
alignmentRelativeError = norm(normalizedOutput - alignedInput) / max(norm(alignedInput), realmin);
assert(alignmentInfo.integerDelay == knownDelay && alignmentRelativeError < 1e-12, 'Feedback alignment did not recover the known delay and gain.');

assert(dpdCfg.blockLength == 30000, 'The common DPD learning block must contain 30000 samples.');
assert(dpdCfg.trainingSampleCount == 300000, 'The DPD training record must contain 300000 samples.');
trainingBlockCount = ceil( dpdCfg.trainingSampleCount / dpdCfg.blockLength);
assert(trainingBlockCount == 10, 'The DPD training record must contain ten learning blocks.');
assert(strcmp(dpdCfg.learningArchitecture, 'indirect') || strcmp(dpdCfg.learningArchitecture, 'direct'), 'The default DPD learning architecture must be direct or indirect.');

results.diagonalCoefficientCounts = actualCounts;
results.blockRelativeError = blockRelativeError;
results.alignmentDelay = alignmentInfo.integerDelay;
results.alignmentRelativeError = alignmentRelativeError;
results.learningBlockLength = dpdCfg.blockLength;
results.trainingSampleCount = dpdCfg.trainingSampleCount;
results.trainingBlockCount = trainingBlockCount;
results.defaultLearningArchitecture = dpdCfg.learningArchitecture;
results.passed = true;

fprintf('DPD self-test passed.\n');
fprintf('  diagonal counts: [%s].\n', sprintf('%d ', actualCounts));
fprintf('  block relative error: %.3g.\n', blockRelativeError);
fprintf('  alignment relative error: %.3g.\n', alignmentRelativeError);
fprintf('  learning blocks: %d x %d samples.\n', trainingBlockCount, dpdCfg.blockLength);
end
