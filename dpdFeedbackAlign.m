function [alignedPaInput, normalizedPaOutput, alignInfo] = dpdFeedbackAlign(paInput, paOutput, dpdCfg)

if nargin < 3
    dpdCfg = struct();
end

if ~isnumeric(paInput) || ~isvector(paInput) || ~isnumeric(paOutput) || ~isvector(paOutput)
    error('dpdFeedbackAlign:InvalidSignal', 'paInput and paOutput must be numeric vectors.');
end

paInput = paInput(:);
paOutput = paOutput(:);

if isempty(paInput) || isempty(paOutput)
    error('dpdFeedbackAlign:EmptySignal', 'The feedback signals must not be empty.');
end

maximumDelay = configurationValue( dpdCfg, 'feedbackAlignmentMaximumDelay', 64);
minimumSamples = configurationValue( dpdCfg, 'minimumFeedbackSamples', 128);
removeDc = configurationValue(dpdCfg, 'feedbackRemoveDc', true);

validateattributes(maximumDelay, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative', 'integer'});
validateattributes(minimumSamples, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'});

maximumPossibleOverlap = min(length(paInput), length(paOutput));
if maximumPossibleOverlap < 2
    error('dpdFeedbackAlign:InsufficientSamples', 'At least two samples are required for alignment.');
end

minimumOverlap = min(minimumSamples, maximumPossibleOverlap);
maximumDelay = min(maximumDelay, maximumPossibleOverlap - minimumOverlap);

correlationInput = paInput - mean(paInput);
correlationOutput = paOutput - mean(paOutput);

candidateDelays = (-maximumDelay:maximumDelay).';
correlationScores = zeros(length(candidateDelays), 1);

for delayIndex = 1:length(candidateDelays)
    integerDelay = candidateDelays(delayIndex);
    [inputSegment, outputSegment] = alignedSegments( correlationInput, correlationOutput, integerDelay);

    inputEnergy = sum(abs(inputSegment).^2);
    outputEnergy = sum(abs(outputSegment).^2);
    energyProduct = inputEnergy * outputEnergy;

    if energyProduct > 0 && isfinite(energyProduct)
        correlationScores(delayIndex) = abs( inputSegment' * outputSegment) / sqrt(energyProduct);
    else
        correlationScores(delayIndex) = 0;
    end
end

[bestCorrelationScore, bestDelayIndex] = max(correlationScores);
integerDelay = candidateDelays(bestDelayIndex);
[alignedPaInput, alignedPaOutput] = alignedSegments(paInput, paOutput, integerDelay);

inputDc = mean(alignedPaInput);
outputDc = mean(alignedPaOutput);

if removeDc
    alignedPaInput = alignedPaInput - inputDc;
    alignedPaOutput = alignedPaOutput - outputDc;
end

inputEnergy = sum(abs(alignedPaInput).^2);
if ~isfinite(inputEnergy) || inputEnergy <= eps
    error('dpdFeedbackAlign:ZeroInputEnergy', 'The aligned PA input has insufficient energy.');
end

linearGain = (alignedPaInput' * alignedPaOutput) / inputEnergy;
if ~isfinite(linearGain) || abs(linearGain) <= sqrt(eps)
    error('dpdFeedbackAlign:InvalidLinearGain', 'The feedback linear gain is zero or nonfinite.');
end

normalizedPaOutput = alignedPaOutput / linearGain;
linearResidual = normalizedPaOutput - alignedPaInput;

alignInfo.integerDelay = integerDelay;
alignInfo.correlationScore = bestCorrelationScore;
alignInfo.candidateDelays = candidateDelays;
alignInfo.correlationScores = correlationScores;
alignInfo.linearGain = linearGain;
alignInfo.removeDc = logical(removeDc);
alignInfo.inputDc = inputDc;
alignInfo.outputDc = outputDc;
alignInfo.numberOfSamples = length(alignedPaInput);
alignInfo.linearNmseDb = 10 * log10( sum(abs(linearResidual).^2) / sum(abs(alignedPaInput).^2));

end

function [inputSegment, outputSegment] = alignedSegments(inputSignal, outputSignal, integerDelay)
if integerDelay >= 0
    overlapLength = min(length(inputSignal), length(outputSignal) - integerDelay);
    inputSegment = inputSignal(1:overlapLength);
    outputSegment = outputSignal( 1+integerDelay:integerDelay+overlapLength);
else
    inputAdvance = -integerDelay;
    overlapLength = min(length(inputSignal) - inputAdvance, length(outputSignal));
    inputSegment = inputSignal( 1+inputAdvance:inputAdvance+overlapLength);
    outputSegment = outputSignal(1:overlapLength);
end
end

function value = configurationValue(configuration, fieldName, defaultValue)
if isfield(configuration, fieldName)
    value = configuration.(fieldName);
else
    value = defaultValue;
end
end
