function signature = cfrSignature(cfrCfg)
    if nargin < 1 || isempty(cfrCfg)
        cfrCfg = cfrConfig(false);
    end
    if ~isstruct(cfrCfg) || ~isscalar(cfrCfg)
        error('cfrSignature:InvalidConfiguration', 'cfrCfg must be a scalar structure.');
    end

    enabled = configurationLogical(cfrCfg, 'enabled', false);
    signature.contractVersion = 1;
    signature.enabled = enabled;

    if ~enabled
        return
    end

    requiredFields = {'algorithm', 'algorithmVersion', 'pulseShape', 'targetPaprDb', 'maximumPulsesPerSymbol', 'cancellationFactor', 'lineSearchFactors', 'relativeTolerance', 'preserveAveragePower', 'fftSize', 'activeIndices'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(cfrCfg, requiredFields{fieldIndex})
            error('cfrSignature:MissingConfigurationField', 'Enabled cfrCfg must contain %s.', requiredFields{fieldIndex});
        end
    end

    if ~ischar(cfrCfg.algorithm) || isempty(strtrim(cfrCfg.algorithm))
        error('cfrSignature:InvalidAlgorithm', 'cfrCfg.algorithm must be a nonempty character vector.');
    end
    if ~ischar(cfrCfg.pulseShape) || isempty(strtrim(cfrCfg.pulseShape))
        error('cfrSignature:InvalidPulseShape', 'cfrCfg.pulseShape must be a nonempty character vector.');
    end
    validatePositiveInteger(cfrCfg.algorithmVersion, 'algorithmVersion');
    validatePositiveInteger(cfrCfg.maximumPulsesPerSymbol, 'maximumPulsesPerSymbol');
    validatePositiveInteger(cfrCfg.fftSize, 'fftSize');
    validateFiniteScalar(cfrCfg.targetPaprDb, 'targetPaprDb');
    validateFiniteScalar(cfrCfg.cancellationFactor, 'cancellationFactor');
    validateFiniteScalar(cfrCfg.relativeTolerance, 'relativeTolerance');
    if cfrCfg.targetPaprDb <= 0
        error('cfrSignature:InvalidTargetPapr', 'cfrCfg.targetPaprDb must be positive.');
    end
    if cfrCfg.cancellationFactor <= 0 || cfrCfg.cancellationFactor > 1
        error('cfrSignature:InvalidCancellationFactor', 'cfrCfg.cancellationFactor must be in (0, 1].');
    end
    if cfrCfg.relativeTolerance < 0
        error('cfrSignature:InvalidTolerance', 'cfrCfg.relativeTolerance must be nonnegative.');
    end
    if ~isnumeric(cfrCfg.lineSearchFactors) || isempty(cfrCfg.lineSearchFactors) || any(~isfinite(cfrCfg.lineSearchFactors(:))) || any(cfrCfg.lineSearchFactors(:) <= 0) || any(cfrCfg.lineSearchFactors(:) > 1)
        error('cfrSignature:InvalidLineSearchFactors', 'cfrCfg.lineSearchFactors must contain finite values in (0, 1].');
    end
    if ~isnumeric(cfrCfg.activeIndices) || isempty(cfrCfg.activeIndices) || any(~isfinite(cfrCfg.activeIndices(:))) || any(cfrCfg.activeIndices(:) ~= round(cfrCfg.activeIndices(:))) || any(cfrCfg.activeIndices(:) < 1) || any(cfrCfg.activeIndices(:) > cfrCfg.fftSize)
        error('cfrSignature:InvalidActiveIndices', 'cfrCfg.activeIndices must contain valid FFT-bin indices.');
    end

    signature.algorithm = lower(strtrim(cfrCfg.algorithm));
    signature.algorithmVersion = cfrCfg.algorithmVersion;
    signature.pulseShape = lower(strtrim(cfrCfg.pulseShape));
    signature.targetPaprDb = cfrCfg.targetPaprDb;
    signature.maximumPulsesPerSymbol = cfrCfg.maximumPulsesPerSymbol;
    signature.cancellationFactor = cfrCfg.cancellationFactor;
    signature.lineSearchFactors = cfrCfg.lineSearchFactors(:).';
    signature.relativeTolerance = cfrCfg.relativeTolerance;
    signature.preserveAveragePower = configurationLogical(cfrCfg, 'preserveAveragePower', true);
    signature.fftSize = cfrCfg.fftSize;
    signature.activeIndices = unique(cfrCfg.activeIndices(:));
    if cfrCfg.algorithmVersion >= 2
        frameFields = {'processingPoint', 'boundaryMode', 'cpLengths', 'wolaLength'};
        for k = 1:length(frameFields)
            if ~isfield(cfrCfg, frameFields{k})
                error('cfrSignature:MissingFrameConfiguration', 'Enabled frame CFR requires %s.', frameFields{k});
            end
        end
        if ~isnumeric(cfrCfg.cpLengths) || isempty(cfrCfg.cpLengths) || ~isreal(cfrCfg.cpLengths) || ...
                any(~isfinite(cfrCfg.cpLengths(:))) || any(cfrCfg.cpLengths(:) < 0) || ...
                any(cfrCfg.cpLengths(:) ~= round(cfrCfg.cpLengths(:))) || any(cfrCfg.cpLengths(:) > cfrCfg.fftSize)
            error('cfrSignature:InvalidCpLengths', 'Invalid CP lengths.');
        end
        validateFiniteScalar(cfrCfg.wolaLength, 'wolaLength');
        if cfrCfg.wolaLength < 0 || cfrCfg.wolaLength ~= round(cfrCfg.wolaLength) || cfrCfg.wolaLength > min(cfrCfg.cpLengths)
            error('cfrSignature:InvalidWolaLength', 'WOLA length must be an integer no greater than the shortest CP.');
        end
        signature.contractVersion = 2;
        signature.processingPoint = cfrCfg.processingPoint;
        signature.boundaryMode = cfrCfg.boundaryMode;
        signature.cpLengths = cfrCfg.cpLengths(:);
        signature.wolaLength = cfrCfg.wolaLength;
    end
end

function value = configurationLogical(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        candidate = configuration.(fieldName);
    else
        candidate = defaultValue;
    end
    if ~(islogical(candidate) || isnumeric(candidate)) || ~isscalar(candidate) || ~isreal(candidate) || ~isfinite(candidate) || (candidate ~= 0 && candidate ~= 1)
        error('cfrSignature:InvalidLogicalField', 'cfrCfg.%s must be a logical scalar.', fieldName);
    end
    value = logical(candidate);
end

function validatePositiveInteger(value, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0 || value ~= round(value)
        error('cfrSignature:InvalidPositiveInteger', 'cfrCfg.%s must be a positive integer.', fieldName);
    end
end

function validateFiniteScalar(value, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('cfrSignature:InvalidFiniteScalar', 'cfrCfg.%s must be a finite real scalar.', fieldName);
    end
end
