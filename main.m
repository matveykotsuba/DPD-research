function project = main(cfg, seed)
    if nargin < 1 || isempty(cfg)
        cfg = config();
    end

    if nargin < 2 || isempty(seed)
        seed = 7;
    end

    if ~isfield(cfg, 'cfr') || isempty(cfg.cfr)
        cfg.cfr = cfrConfig(false);
    end
    cfg.cfr.fftSize = cfg.fftSize;
    cfg.cfr.activeIndices = cfg.activeIndices;
    cfg.cfr.cpLengths = cfg.cpLengths;
    cfg.cfr.wolaLength = cfg.wolaLength;

    if isfield(cfg, 'dpd') && isstruct(cfg.dpd)
        [cfg.dpd, ~] = dpdRefreshCfrCompatibility(cfg.dpd, cfg.cfr);
        cfg.dpd.adaptationMetricSeed = seed;
        cfg.dpd.adaptationMetricOfdmSymbols = cfg.numOfdmSymbols;
    end

    if isfield(cfg, 'dpdLearningArchitecture') && ~strcmpi(strtrim(cfg.dpdLearningArchitecture), cfg.dpd.learningArchitecture)
        error('main:DpdArchitectureConfigurationMismatch', 'Use config(''direct'') or config(''indirect'') so the architecture and coefficient file are selected together.');
    end

    dpdTrainingInfo = struct();
    dpdEnabled = isfield(cfg, 'dpd') && isfield(cfg.dpd, 'enabled') && cfg.dpd.enabled;
    dpdAdaptationEnabled = dpdEnabled && isfield(cfg.dpd, 'adaptationEnabled') && cfg.dpd.adaptationEnabled;
    dpdModelStale = dpdEnabled && isfield(cfg.dpd, 'modelStale') && cfg.dpd.modelStale;
    dpdModelLoaded = dpdEnabled && isfield(cfg.dpd, 'modelLoaded') && cfg.dpd.modelLoaded;
    dpdInMemoryModelVerified = dpdEnabled && isfield(cfg.dpd, 'trained') && cfg.dpd.trained && isfield(cfg.dpd, 'modelVerified') && cfg.dpd.modelVerified;
    dpdModelUsable = ~dpdModelStale && (dpdModelLoaded || dpdInMemoryModelVerified);
    dpdTrainingRequired = dpdEnabled && (dpdAdaptationEnabled || ~dpdModelUsable);
    if dpdTrainingRequired
        if ~dpdAdaptationEnabled
            fprintf('No compatible %s DPD model is available; running automatic training.\n', upper(cfg.dpd.learningArchitecture));
        end
        [cfg.dpd, dpdTrainingInfo] = dpdTrainRun(cfg);
        cfg.dpd.enabled = true;
    end

    dpdAdaptationMetrics = struct();
    if dpdEnabled && isfield(dpdTrainingInfo, 'blockMetrics')
        dpdAdaptationMetrics = dpdTrainingInfo.blockMetrics;
    elseif dpdEnabled && isfield(cfg.dpd, 'savedTrainingInfo') && isfield(cfg.dpd.savedTrainingInfo, 'blockMetrics')
        dpdAdaptationMetrics = cfg.dpd.savedTrainingInfo.blockMetrics;
    end
    dpdArchitectureComparisonMetrics = struct();
    if dpdEnabled && isfield(dpdTrainingInfo, 'architectureComparisonMetrics')
        [dpdArchitectureComparisonMetrics, comparisonMetricsAvailable] = dpdMakeArchitectureComparisonMetrics(dpdTrainingInfo);
        if ~comparisonMetricsAvailable && isfield(cfg.dpd, 'savedTrainingInfo')
            [dpdArchitectureComparisonMetrics, comparisonMetricsAvailable] = dpdMakeArchitectureComparisonMetrics(cfg.dpd.savedTrainingInfo);
        end
        if ~comparisonMetricsAvailable
            dpdArchitectureComparisonMetrics = struct();
        end
    elseif dpdEnabled && isfield(cfg.dpd, 'savedTrainingInfo')
        [dpdArchitectureComparisonMetrics, comparisonMetricsAvailable] = dpdMakeArchitectureComparisonMetrics(cfg.dpd.savedTrainingInfo);
        if ~comparisonMetricsAvailable
            dpdArchitectureComparisonMetrics = struct();
        end
    end
    [txSignal, txInfo] = tx(cfg, seed);

    cfrPeakCancellationFigure = [];
    cfrAnalysisPlotEnabled = isfield(cfg, 'cfr') && isfield(cfg.cfr, 'analysisPlotEnabled') && logical(cfg.cfr.analysisPlotEnabled);
    if isfield(txInfo, 'cfr') && isfield(txInfo.cfr, 'enabled') && txInfo.cfr.enabled && cfrAnalysisPlotEnabled
        cfrPeakCancellationFigure = cfrPeakCancellationPlot(txInfo.cfr, cfg.cfr);
    end

    prePa = struct();
    prePa.signal = txSignal;
    prePa.frequency = txInfo.frequency;
    prePa.psd = txInfo.psd;
    prePa.psdInfo = txInfo.psdInfo;
    prePa.aclr = txInfo.aclr;
    prePa.averagePower = txInfo.averagePower;
    prePa.peakPower = txInfo.peakPower;
    prePa.paprDb = txInfo.paprDb;

    [paSignal, paInfo] = pa(txSignal, txInfo.cfg);

    postAclrWelchSegmentLength = 8192;
    if isfield(cfg, 'dpd') && isfield(cfg.dpd, 'adaptationMetricSpectrumSegmentLength')
        postAclrWelchSegmentLength = cfg.dpd.adaptationMetricSpectrumSegmentLength;
    end
    [postFrequency, postPsd, postPsdInfo] = welch(paSignal, txInfo.cfg.sampleRate, postAclrWelchSegmentLength);
    postAclr = aclr(postFrequency, postPsd, txInfo.cfg.channelBandwidth, txInfo.cfg.aclrMeasurementBandwidth);

    postPa = struct();
    postPa.signal = paSignal;
    postPa.frequency = postFrequency;
    postPa.psd = postPsd;
    postPa.psdInfo = postPsdInfo;
    postPa.aclr = postAclr;
    postPa.averagePower = mean(abs(paSignal).^2);
    postPa.peakPower = max(abs(paSignal).^2);
    postPa.paprDb = 10 * log10(postPa.peakPower / postPa.averagePower);

    if isfield(cfg, 'dpd')
        [dpdArchitectureComparisonFigure, dpdArchitectureComparison] = dpdArchitectureComparisonPlot(cfg.dpd, dpdArchitectureComparisonMetrics, postAclr);
        if isfield(dpdArchitectureComparison, 'available') && dpdArchitectureComparison.available
            activeArchitecture = lower(strtrim(cfg.dpd.learningArchitecture));
            if isfield(dpdArchitectureComparison, activeArchitecture)
                dpdAdaptationMetrics = dpdArchitectureComparison.(activeArchitecture);
            end
        end
        if ~blockMetricsMatchMainAclr(dpdAdaptationMetrics, cfg.dpd, postPsdInfo)
            if isstruct(dpdAdaptationMetrics) && ~isempty(fieldnames(dpdAdaptationMetrics))
                fprintf('DPD block plot was not created because its reporting history does not match the current full-frame Welch ACLR measurement.\n');
            end
            dpdAdaptationMetrics = struct();
        end
        dpdAdaptationFigure = dpdAdaptationPlot(dpdAdaptationMetrics, cfg.dpd, postAclr);
    else
        disabledDpdCfg = struct();
        disabledDpdCfg.enabled = false;
        [dpdArchitectureComparisonFigure, dpdArchitectureComparison] = dpdArchitectureComparisonPlot(disabledDpdCfg, struct(), postAclr);
        dpdAdaptationFigure = dpdAdaptationPlot(dpdAdaptationMetrics, disabledDpdCfg, postAclr);
    end

    txInfo.prePa = prePa;
    txInfo.postPa = postPa;
    txInfo.paInputSignal = txSignal;
    txInfo.paSignal = paSignal;
    txInfo.paInfo = paInfo;
    txInfo.txSignal = paSignal;
    txInfo.frequency = postFrequency;
    txInfo.psd = postPsd;
    txInfo.psdInfo = postPsdInfo;
    txInfo.aclr = postAclr;
    txInfo.averagePower = postPa.averagePower;
    txInfo.peakPower = postPa.peakPower;
    txInfo.paprDb = postPa.paprDb;

    [rxSignal, channelInfo] = channel(paSignal, txInfo.cfg);
    results = rx(rxSignal, txInfo, channelInfo);

    if ~isempty(dpdArchitectureComparisonFigure) && ishghandle(dpdArchitectureComparisonFigure)
        figure(dpdArchitectureComparisonFigure);
        drawnow;
    elseif dpdEnabled && isfield(dpdArchitectureComparison, 'available') && ~dpdArchitectureComparison.available
        fprintf('DPD Direct/Indirect comparison plot was not created: %s.\n', dpdArchitectureComparison.reason);
        if isfield(dpdArchitectureComparison, 'status')
            fprintf('  Indirect [%s]: %s.\n', dpdArchitectureComparison.status.indirect.source, dpdArchitectureComparison.status.indirect.reason);
            fprintf('  Direct [%s]: %s.\n', dpdArchitectureComparison.status.direct.source, dpdArchitectureComparison.status.direct.reason);
        end
    end

    project.txSignal = txSignal;
    project.txInfo = txInfo;
    project.cfrInfo = txInfo.cfr;
    project.cfrPeakCancellationFigure = cfrPeakCancellationFigure;
    project.prePa = prePa;
    project.paSignal = paSignal;
    project.paInfo = paInfo;
    project.dpdInfo = paInfo.dpd;
    project.dpdTrainingInfo = dpdTrainingInfo;
    project.dpdAdaptationMetrics = dpdAdaptationMetrics;
    project.dpdAdaptationFigure = dpdAdaptationFigure;
    project.dpdArchitectureComparisonMetrics = dpdArchitectureComparisonMetrics;
    project.dpdArchitectureComparison = dpdArchitectureComparison;
    project.dpdArchitectureComparisonFigure = dpdArchitectureComparisonFigure;
    project.postPa = postPa;
    project.rxSignal = rxSignal;
    project.channelInfo = channelInfo;
    project.results = results;
end

function matches = blockMetricsMatchMainAclr(metrics, dpdCfg, postPsdInfo)
    matches = false;
    requiredFields = {'monitorSeed', 'monitorOfdmSymbolCount', 'spectrumEstimator', 'spectrumWelchInfo'};
    if ~isstruct(metrics) || ~isscalar(metrics) || isempty(fieldnames(metrics))
        return
    end
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(metrics, requiredFields{fieldIndex})
            return
        end
    end
    if metrics.monitorSeed ~= dpdCfg.adaptationMetricSeed || metrics.monitorOfdmSymbolCount ~= dpdCfg.adaptationMetricOfdmSymbols
        return
    end
    if ~ischar(metrics.spectrumEstimator) || ~strcmpi(strtrim(metrics.spectrumEstimator), 'Welch')
        return
    end
    if ~isstruct(metrics.spectrumWelchInfo) || ~isscalar(metrics.spectrumWelchInfo) || ~isfield(metrics.spectrumWelchInfo, 'segmentLength')
        return
    end
    if ~isstruct(postPsdInfo) || ~isfield(postPsdInfo, 'segmentLength') || metrics.spectrumWelchInfo.segmentLength ~= postPsdInfo.segmentLength
        return
    end
    matches = true;
end
