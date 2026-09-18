function [metrics, candidateCfg] = dpdMonitorCandidateMetrics(trainingInput, candidateCfg, paCfg, monitoringSignal, includeSpectrum)
    if nargin < 4
        error('dpdMonitorCandidateMetrics:MissingInput', 'trainingInput, candidateCfg, paCfg and monitoringSignal are required.');
    end
    if nargin < 5 || isempty(includeSpectrum)
        includeSpectrum = false;
    end
    validateattributes(includeSpectrum, {'logical', 'numeric'}, {'scalar'});
    includeSpectrum = logical(includeSpectrum);
    requiredMonitorFields = {'cfg', 'systemInput', 'txInfo', 'modelInput', 'inputInfo', 'stressInput'};
    for fieldIndex = 1:length(requiredMonitorFields)
        if ~isfield(monitoringSignal, requiredMonitorFields{fieldIndex})
            error('dpdMonitorCandidateMetrics:MissingMonitoringField', 'monitoringSignal must contain %s.', requiredMonitorFields{fieldIndex});
        end
    end
    requiredInputInfoFields = {'inputScale', 'rateFilter'};
    for fieldIndex = 1:length(requiredInputInfoFields)
        if ~isfield(monitoringSignal.inputInfo, requiredInputInfoFields{fieldIndex})
            error('dpdMonitorCandidateMetrics:MissingInputInfoField', 'monitoringSignal.inputInfo must contain %s.', requiredInputInfoFields{fieldIndex});
        end
    end
    requiredSystemFields = {'sampleRate', 'channelBandwidth', 'aclrMeasurementBandwidth', 'fftSize', 'cpLengths', 'numOfdmSymbols', 'numActiveSubcarriers', 'activeIndices'};
    for fieldIndex = 1:length(requiredSystemFields)
        if ~isfield(monitoringSignal.cfg, requiredSystemFields{fieldIndex})
            error('dpdMonitorCandidateMetrics:MissingSystemField', 'monitoringSignal.cfg must contain %s.', requiredSystemFields{fieldIndex});
        end
    end
    if ~isfield(monitoringSignal.txInfo, 'qamSymbols')
        error('dpdMonitorCandidateMetrics:MissingQamReference', 'monitoringSignal.txInfo.qamSymbols is required.');
    end
    trainingInput = trainingInput(:);
    monitorSystemInput = monitoringSignal.systemInput(:);
    monitorModelInput = monitoringSignal.modelInput(:);
    [candidateCfg, normalizationInfo] = dpdNormalizeCandidate(trainingInput, candidateCfg, paCfg);
    stressGuard = dpdStressGuard(candidateCfg, paCfg, monitoringSignal.stressInput);
    [monitorDpdOutput, dpdInfo] = dpdCore(monitorModelInput, candidateCfg);
    monitorPaOutput = gmpCore(monitorDpdOutput, paCfg);
    [monitorSystemOutput, systemOutputInfo] = paSystemOutput(monitorPaOutput, length(monitorSystemInput), monitoringSignal.inputInfo.inputScale, monitoringSignal.inputInfo.rateFilter, monitoringSignal.cfg);
    [evmRmsPercent, evmRmsDb] = qamEvm(monitorSystemOutput, monitoringSignal.txInfo.qamSymbols, monitoringSignal.cfg);
    spectrumSegmentLength = round(configurationValue(candidateCfg, 'adaptationMetricSpectrumSegmentLength', 8192));
    validateattributes(spectrumSegmentLength, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    if spectrumSegmentLength < 2
        error('dpdMonitorCandidateMetrics:InvalidWelchSegmentLength', 'adaptationMetricSpectrumSegmentLength must be at least two.');
    end
    [outputFrequency, outputPsd, outputWelchInfo] = welch(monitorSystemOutput, monitoringSignal.cfg.sampleRate, spectrumSegmentLength);
    [prePaFrequency, prePaPsd, prePaWelchInfo] = welch(monitorSystemInput, monitoringSignal.cfg.sampleRate, spectrumSegmentLength);
    if length(outputFrequency) ~= length(prePaFrequency) || any(outputFrequency ~= prePaFrequency)
        error('dpdMonitorCandidateMetrics:InconsistentSpectrumGrid', 'The pre-PA and DPD-plus-PA Welch spectra must use the same frequency grid.');
    end
    outputAclr = aclr(outputFrequency, outputPsd, monitoringSignal.cfg.channelBandwidth, monitoringSignal.cfg.aclrMeasurementBandwidth);
    prePaAclr = aclr(prePaFrequency, prePaPsd, monitoringSignal.cfg.channelBandwidth, monitoringSignal.cfg.aclrMeasurementBandwidth);
    metrics.normalization = normalizationInfo;
    metrics.coefficientScale = normalizationInfo.coefficientScale;
    metrics.stressGuard = stressGuard;
    metrics.evmRmsPercent = evmRmsPercent;
    metrics.evmRmsDb = evmRmsDb;
    metrics.worstAclr1Db = outputAclr.worstAclr1Db;
    metrics.worstAclr2Db = outputAclr.worstAclr2Db;
    metrics.prePaWorstAclr1Db = prePaAclr.worstAclr1Db;
    metrics.prePaWorstAclr2Db = prePaAclr.worstAclr2Db;
    metrics.prePaAclr1Db = prePaAclr.worstAclr1Db;
    metrics.prePaAclr2Db = prePaAclr.worstAclr2Db;
    metrics.inputWorstAclr1Db = prePaAclr.worstAclr1Db;
    metrics.inputWorstAclr2Db = prePaAclr.worstAclr2Db;
    metrics.aclr1ShortfallDb = prePaAclr.worstAclr1Db - outputAclr.worstAclr1Db;
    metrics.aclr2ShortfallDb = prePaAclr.worstAclr2Db - outputAclr.worstAclr2Db;
    metrics.aclr1GapToPrePaDb = outputAclr.worstAclr1Db - prePaAclr.worstAclr1Db;
    metrics.aclr2GapToPrePaDb = outputAclr.worstAclr2Db - prePaAclr.worstAclr2Db;
    metrics.dpdInputAveragePower = mean(abs(monitorModelInput).^2);
    metrics.dpdOutputAveragePower = dpdInfo.outputAveragePower;
    metrics.paModelOutputAveragePower = mean(abs(monitorPaOutput).^2);
    metrics.paSystemOutputAveragePower = systemOutputInfo.systemOutputAveragePower;
    metrics.dpdOutputPeakMagnitude = dpdInfo.outputPeakMagnitude;
    metrics.spectrumEstimator = 'Welch';
    metrics.spectrumSegmentLength = spectrumSegmentLength;
    metrics.spectrumWelchInfo = outputWelchInfo;
    metrics.prePaSpectrumWelchInfo = prePaWelchInfo;
    metrics.spectrumIncluded = includeSpectrum;
    chainFinite = all(isfinite(monitorDpdOutput)) && all(isfinite(monitorPaOutput)) && all(isfinite(monitorSystemOutput));
    metricFinite = isfinite(evmRmsPercent) && isfinite(evmRmsDb) && isfinite(outputAclr.worstAclr1Db) && isfinite(outputAclr.worstAclr2Db) && isfinite(prePaAclr.worstAclr1Db) && isfinite(prePaAclr.worstAclr2Db);
    metrics.chainFinite = chainFinite;
    metrics.metricFinite = metricFinite;
    metrics.isSafe = normalizationInfo.succeeded && stressGuard.isSafe && chainFinite && metricFinite;
    metrics.evmReference = 'QAM symbols after CP removal, FFT and best complex-gain correction';
    metrics.aclrReference = 'Welch PSD of the PA system output on the fixed monitoring signal';
    metrics.prePaAclrReference = 'Welch PSD of monitoringSignal.systemInput before PA input processing';
    if includeSpectrum
        metrics.spectrumFrequencyHz = outputFrequency;
        metrics.outputPsd = outputPsd;
        metrics.prePaPsd = prePaPsd;
        metrics.outputAclr = outputAclr;
        metrics.prePaAclr = prePaAclr;
    end
end

function [evmRmsPercent, evmRmsDb] = qamEvm(systemOutput, referenceQamSymbols, cfg)
    symbolLengths = cfg.fftSize + cfg.cpLengths;
    requiredSampleCount = sum(symbolLengths);
    if length(systemOutput) < requiredSampleCount
        error('dpdMonitorCandidateMetrics:ShortSystemOutput', 'The PA system output is shorter than the monitoring OFDM record.');
    end
    receivedNoCp = complex(zeros(cfg.fftSize, cfg.numOfdmSymbols));
    readIndex = 1;
    for symbolIndex = 1:cfg.numOfdmSymbols
        usefulStart = readIndex + cfg.cpLengths(symbolIndex);
        usefulStop = usefulStart + cfg.fftSize - 1;
        receivedNoCp(:, symbolIndex) = systemOutput(usefulStart:usefulStop);
        readIndex = usefulStop + 1;
    end
    receivedGrid = fftshift(fft(receivedNoCp, [], 1) * sqrt(cfg.numActiveSubcarriers) / cfg.fftSize, 1);
    receivedQamMatrix = receivedGrid(cfg.activeIndices, :);
    receivedQamSymbols = receivedQamMatrix(:);
    referenceQamSymbols = referenceQamSymbols(:);
    referenceEnergy = sum(abs(referenceQamSymbols).^2);
    bestLinearGain = sum(conj(referenceQamSymbols) .* receivedQamSymbols) / referenceEnergy;
    correctedQamSymbols = receivedQamSymbols / bestLinearGain;
    errorEnergy = sum(abs(correctedQamSymbols - referenceQamSymbols).^2);
    evmRms = sqrt(errorEnergy / referenceEnergy);
    evmRmsPercent = 100 * evmRms;
    evmRmsDb = 20 * log10(max(evmRms, realmin));
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
