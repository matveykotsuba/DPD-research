function [outputSignal, dpdInfo] = dpdCore(inputSignal, dpdCfg)

if nargin < 2 || ~isstruct(dpdCfg)
    error('dpdCore:MissingConfiguration', 'A DPD configuration structure is required.');
end

if ~isfield(dpdCfg, 'enabled')
    error('dpdCore:MissingEnabledFlag', 'dpdCfg.enabled is required.');
end

if ~isnumeric(inputSignal) || ~isvector(inputSignal)
    error('dpdCore:InvalidInput', 'inputSignal must be a numeric vector.');
end

inputSignal = inputSignal(:);

if dpdCfg.enabled
    if isfield(dpdCfg, 'modelStale') && dpdCfg.modelStale
        error('dpdCore:StaleModel', '%s', dpdCfg.modelStaleReason);
    end
    if isfield(dpdCfg, 'modelLoaded') && dpdCfg.modelLoaded && isfield(dpdCfg, 'loadedModelArchitecture') && ~strcmp(dpdCfg.loadedModelArchitecture, dpdCfg.learningArchitecture)
        error('dpdCore:LearningArchitectureMismatch', 'The loaded DPD coefficients do not match learningArchitecture. Recreate the configuration with config(''direct'') or config(''indirect'').');
    end
    validateDpdModel(dpdCfg);
    activeDpdCfg = dpdCfg;
    if isfield(dpdCfg, 'diagonalCount')
        activeMask = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);
        activeDpdCfg.coefficients(~activeMask) = 0;
    end
    outputSignal = gmpCore(inputSignal, activeDpdCfg);
else
    outputSignal = inputSignal;
end

dpdInfo.enabled = logical(dpdCfg.enabled);
dpdInfo.bypassed = ~logical(dpdCfg.enabled);
dpdInfo.inputAveragePower = averagePower(inputSignal);
dpdInfo.outputAveragePower = averagePower(outputSignal);
dpdInfo.inputRms = sqrt(dpdInfo.inputAveragePower);
dpdInfo.outputRms = sqrt(dpdInfo.outputAveragePower);
dpdInfo.outputPeakMagnitude = peakMagnitude(outputSignal);
dpdInfo.outputPaprDb = paprDb(outputSignal);

if isfield(dpdCfg, 'diagonalCount') && isfield(dpdCfg, 'coefficients')
    currentMask = gmpDiagonalMask(dpdCfg, dpdCfg.diagonalCount);
    activeCoefficients = dpdCfg.coefficients(currentMask);
elseif isfield(dpdCfg, 'coefficientMask') && isequal(size(dpdCfg.coefficientMask), size(dpdCfg.coefficients))
    activeCoefficients = dpdCfg.coefficients(logical(dpdCfg.coefficientMask));
elseif isfield(dpdCfg, 'coefficients')
    activeCoefficients = dpdCfg.coefficients(:);
else
    activeCoefficients = complex(zeros(0, 1));
end
dpdInfo.coefficientNorm = norm(activeCoefficients);
dpdInfo.numberOfActiveCoefficients = length(activeCoefficients);
if isfield(dpdCfg, 'diagonalCount')
    dpdInfo.diagonalCount = dpdCfg.diagonalCount;
else
    dpdInfo.diagonalCount = NaN;
end

if isfield(dpdCfg, 'maximumOutputMagnitude') && isfinite(dpdCfg.maximumOutputMagnitude) && ~isempty(outputSignal)
    dpdInfo.outOfRangeFraction = mean( abs(outputSignal) > dpdCfg.maximumOutputMagnitude);
else
    dpdInfo.outOfRangeFraction = 0;
end

end

function validateDpdModel(dpdCfg)
requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'coefficients'};
for fieldIndex = 1:length(requiredFields)
    if ~isfield(dpdCfg, requiredFields{fieldIndex})
        error('dpdCore:MissingConfigurationField', 'dpdCfg must contain %s.', requiredFields{fieldIndex});
    end
end

expectedColumns = 1 + (length(dpdCfg.orders) - 1) * length(dpdCfg.envelopeDelays);
expectedSize = [length(dpdCfg.signalDelays), expectedColumns];

if ~isequal(size(dpdCfg.coefficients), expectedSize)
    error('dpdCore:InvalidCoefficientSize', 'dpdCfg.coefficients must have size %d-by-%d.', expectedSize(1), expectedSize(2));
end

if any(~isfinite(dpdCfg.coefficients(:)))
    error('dpdCore:NonfiniteCoefficients', 'DPD coefficients must be finite.');
end
end

function value = averagePower(signal)
if isempty(signal)
    value = NaN;
else
    value = mean(abs(signal).^2);
end
end

function value = peakMagnitude(signal)
if isempty(signal)
    value = NaN;
else
    value = max(abs(signal));
end
end

function value = paprDb(signal)
powerValue = averagePower(signal);
if isempty(signal) || ~isfinite(powerValue) || powerValue <= 0
    value = NaN;
else
    value = 10 * log10(max(abs(signal).^2) / powerValue);
end
end
