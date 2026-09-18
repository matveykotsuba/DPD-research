function cfrCfg = cfrConfig(enabledOverride)
    configuredEnabled = true;

    if nargin >= 1 && ~isempty(enabledOverride)
        if ~(islogical(enabledOverride) || isnumeric(enabledOverride)) || ~isscalar(enabledOverride) || ~isreal(enabledOverride) || ~isfinite(enabledOverride) || (enabledOverride ~= 0 && enabledOverride ~= 1)
            error('cfrConfig:InvalidEnabledOverride', 'enabledOverride must be a logical scalar.');
        end
        configuredEnabled = logical(enabledOverride);
    end

    cfrCfg.enabled = configuredEnabled;
    cfrCfg.algorithm = 'sequential-peak-cancellation';
    cfrCfg.algorithmVersion = 2;
    cfrCfg.processingPoint = 'after-cp-wola';
    cfrCfg.pulseShape = 'frame-occupied-band-rectangular-mask';
    cfrCfg.boundaryMode = 'periodic-frame';
    cfrCfg.targetPaprDb = 7.5;
    cfrCfg.maximumPulsesPerSymbol = 64;
    cfrCfg.cancellationFactor = 1.0;
    cfrCfg.lineSearchFactors = [1 0.5 0.25 0.125];
    cfrCfg.relativeTolerance = 1e-12;
    cfrCfg.preserveAveragePower = true;
    cfrCfg.analysisPlotEnabled = true;
    cfrCfg.analysisWelchSegmentLength = 8192;
end
