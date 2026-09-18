function results = cfrExperimentRun(userCfg)
    if nargin < 1 || isempty(userCfg)
        userCfg = struct();
    end
    if ~isstruct(userCfg) || ~isscalar(userCfg)
        error('cfrExperimentRun:InvalidConfiguration', 'userCfg must be a scalar structure.');
    end

    experimentCfg = defaultExperimentConfiguration();
    experimentCfg = applyUserConfiguration(experimentCfg, userCfg);
    validateExperimentConfiguration(experimentCfg);

    baseCfg = config();
    experimentCfg.numOfdmSymbols = min(experimentCfg.numOfdmSymbols, length(baseCfg.cpLengths));
    baseCfg.numOfdmSymbols = experimentCfg.numOfdmSymbols;
    baseCfg.cpLengths = baseCfg.cpLengths(1:experimentCfg.numOfdmSymbols);
    baseCfg.cfr.cpLengths = baseCfg.cpLengths;
    baseCfg.cfr.enabled = false;
    if isfield(baseCfg, 'dpd')
        baseCfg.dpd.enabled = false;
        baseCfg.dpd.adaptationEnabled = false;
    end

    cfrCfg = baseCfg;
    cfrCfg.cfr.enabled = true;
    if ~isempty(experimentCfg.targetPaprDb)
        cfrCfg.cfr.targetPaprDb = experimentCfg.targetPaprDb;
    end
    if isfield(cfrCfg, 'dpd')
        [cfrCfg.dpd, ~] = dpdRefreshCfrCompatibility(cfrCfg.dpd, cfrCfg.cfr);
        cfrCfg.dpd.enabled = false;
        cfrCfg.dpd.adaptationEnabled = false;
    end

    fprintf('\nCFR peak-cancellation experiment after CP/WOLA.\n');
    fprintf('Seed = %d, OFDM symbols = %d, target PAPR = %.2f dB.\n', experimentCfg.seed, experimentCfg.numOfdmSymbols, cfrCfg.cfr.targetPaprDb);

    [originalSignal, originalTxInfo] = tx(baseCfg, experimentCfg.seed, false);
    [cfrSignal, cfrTxInfo] = tx(cfrCfg, experimentCfg.seed, false);
    if ~isequal(originalTxInfo.bits, cfrTxInfo.bits) || ~isequal(originalTxInfo.qamSymbols, cfrTxInfo.qamSymbols)
        error('cfrExperimentRun:ReferenceMismatch', 'CFR OFF and ON must use identical source bits and QAM symbols.');
    end

    [frequencyOriginal, psdOriginal, welchInfoOriginal] = welch(originalSignal, baseCfg.sampleRate, experimentCfg.welchSegmentLength);
    [frequencyCfr, psdCfr, welchInfoCfr] = welch(cfrSignal, cfrCfg.sampleRate, experimentCfg.welchSegmentLength);
    if ~isequal(frequencyOriginal, frequencyCfr)
        error('cfrExperimentRun:WelchGridMismatch', 'CFR OFF and ON Welch frequency grids must be identical.');
    end
    originalAclr = aclr(frequencyOriginal, psdOriginal, baseCfg.channelBandwidth, baseCfg.aclrMeasurementBandwidth);
    cfrAclr = aclr(frequencyCfr, psdCfr, cfrCfg.channelBandwidth, cfrCfg.aclrMeasurementBandwidth);

    [originalEvmPercent, originalEvmDb] = qamEvm(originalSignal, originalTxInfo.qamSymbols, baseCfg);
    [cfrEvmPercent, cfrEvmDb] = qamEvm(cfrSignal, cfrTxInfo.qamSymbols, cfrCfg);
    originalPaprDb = signalPaprDb(originalSignal);
    cfrPaprDb = signalPaprDb(cfrSignal);

    interpolated = struct();
    if logical(experimentCfg.measureInterpolatedPapr)
        [originalModelInput, originalModelInfo] = paModelInput(originalSignal, baseCfg);
        [cfrModelInput, cfrModelInfo] = paModelInput(cfrSignal, cfrCfg);
        interpolated.available = true;
        interpolated.sampleRate = baseCfg.pa.modelSampleRate;
        interpolated.originalPaprDb = signalPaprDb(originalModelInput);
        interpolated.cfrPaprDb = signalPaprDb(cfrModelInput);
        interpolated.paprReductionDb = interpolated.originalPaprDb - interpolated.cfrPaprDb;
        interpolated.originalPeakMagnitude = max(abs(originalModelInput));
        interpolated.cfrPeakMagnitude = max(abs(cfrModelInput));
        interpolated.originalRms = originalModelInfo.modelInputRms;
        interpolated.cfrRms = cfrModelInfo.modelInputRms;
        clear originalModelInput cfrModelInput
    else
        interpolated.available = false;
    end

    correction = cfrSignal - originalSignal;
    correctionPowerRatio = mean(abs(correction).^2) / max(mean(abs(originalSignal).^2), realmin('double'));
    averagePowerChangeDb = powerRatioDb(mean(abs(cfrSignal).^2), mean(abs(originalSignal).^2));
    worstAclr1OriginalDb = min(originalAclr.aclr1LeftDb, originalAclr.aclr1RightDb);
    worstAclr1CfrDb = min(cfrAclr.aclr1LeftDb, cfrAclr.aclr1RightDb);
    worstAclr2OriginalDb = min(originalAclr.aclr2LeftDb, originalAclr.aclr2RightDb);
    worstAclr2CfrDb = min(cfrAclr.aclr2LeftDb, cfrAclr.aclr2RightDb);

    figureHandles.amplitude = [];
    figureHandles.peakCancellation = [];
    figureHandles.ccdfWelch = [];
    if logical(experimentCfg.showPlots)
        figureHandles.amplitude = plotAmplitudeAndSymbolPapr(cfrTxInfo.cfr, cfrCfg.cfr);
        figureHandles.peakCancellation = cfrPeakCancellationPlot(cfrTxInfo.cfr, cfrCfg.cfr);
        figureHandles.ccdfWelch = plotCcdfAndWelch(cfrTxInfo.cfr, frequencyOriginal, psdOriginal, psdCfr, originalPaprDb, cfrPaprDb, experimentCfg);
    end

    fprintf('Frame PAPR after WOLA, CFR OFF -> ON: %.3f -> %.3f dB (reduction %.3f dB).\n', originalPaprDb, cfrPaprDb, originalPaprDb - cfrPaprDb);
    if interpolated.available
        fprintf('Post-interpolation:   %.3f -> %.3f dB (reduction %.3f dB).\n', interpolated.originalPaprDb, interpolated.cfrPaprDb, interpolated.paprReductionDb);
    end
    fprintf('CFR-only QAM EVM: %.4f %% (baseline %.4f %%).\n', cfrEvmPercent, originalEvmPercent);
    fprintf('Worst ACLR1: %.3f -> %.3f dB; worst ACLR2: %.3f -> %.3f dB.\n', worstAclr1OriginalDb, worstAclr1CfrDb, worstAclr2OriginalDb, worstAclr2CfrDb);
    fprintf('Compensating pulses: %d; correction power = %.3f dB relative to OFDM.\n', cfrTxInfo.cfr.totalPulseCount, 10 * log10(max(correctionPowerRatio, realmin('double'))));

    results.configuration = experimentCfg;
    results.cfrConfiguration = cfrCfg.cfr;
    results.processingPoint = cfrTxInfo.cfr.processingPoint;
    results.spectrumEstimator = 'Welch';
    results.original.signal = originalSignal;
    results.original.txInfo = originalTxInfo;
    results.original.paprDb = originalPaprDb;
    results.original.evmRmsPercent = originalEvmPercent;
    results.original.evmRmsDb = originalEvmDb;
    results.original.frequency = frequencyOriginal;
    results.original.psd = psdOriginal;
    results.original.welchInfo = welchInfoOriginal;
    results.original.aclr = originalAclr;
    results.cfr.signal = cfrSignal;
    results.cfr.txInfo = cfrTxInfo;
    results.cfr.paprDb = cfrPaprDb;
    results.cfr.paprReductionDb = originalPaprDb - cfrPaprDb;
    results.cfr.evmRmsPercent = cfrEvmPercent;
    results.cfr.evmRmsDb = cfrEvmDb;
    results.cfr.frequency = frequencyCfr;
    results.cfr.psd = psdCfr;
    results.cfr.welchInfo = welchInfoCfr;
    results.cfr.aclr = cfrAclr;
    results.cfr.correctionPowerDb = 10 * log10(max(correctionPowerRatio, realmin('double')));
    results.cfr.averagePowerChangeDb = averagePowerChangeDb;
    results.cfr.outputPaprConsistencyErrorDb = cfrPaprDb - cfrTxInfo.cfr.outputPaprDb;
    results.cfr.worstAclr1ChangeDb = worstAclr1CfrDb - worstAclr1OriginalDb;
    results.cfr.worstAclr2ChangeDb = worstAclr2CfrDb - worstAclr2OriginalDb;
    results.interpolated = interpolated;
    results.figures = figureHandles;
end

function cfg = defaultExperimentConfiguration()
    cfg.seed = 7;
    cfg.numOfdmSymbols = 140;
    cfg.targetPaprDb = [];
    cfg.welchSegmentLength = 8192;
    cfg.displayFloorDb = -100;
    cfg.measureInterpolatedPapr = true;
    cfg.showPlots = true;
end

function cfg = applyUserConfiguration(cfg, userCfg)
    userFields = fieldnames(userCfg);
    for fieldIndex = 1:length(userFields)
        fieldName = userFields{fieldIndex};
        if ~isfield(cfg, fieldName)
            error('cfrExperimentRun:UnknownConfigurationField', 'Unknown configuration field: %s.', fieldName);
        end
        cfg.(fieldName) = userCfg.(fieldName);
    end
end

function validateExperimentConfiguration(cfg)
    validateNonnegativeInteger(cfg.seed, 'seed');
    validatePositiveInteger(cfg.numOfdmSymbols, 'numOfdmSymbols');
    validatePositiveInteger(cfg.welchSegmentLength, 'welchSegmentLength');
    if ~isempty(cfg.targetPaprDb) && (~isnumeric(cfg.targetPaprDb) || ~isscalar(cfg.targetPaprDb) || ~isreal(cfg.targetPaprDb) || ~isfinite(cfg.targetPaprDb) || cfg.targetPaprDb <= 0)
        error('cfrExperimentRun:InvalidTargetPapr', 'targetPaprDb must be empty or a positive finite scalar.');
    end
    if ~isnumeric(cfg.displayFloorDb) || ~isscalar(cfg.displayFloorDb) || ~isreal(cfg.displayFloorDb) || ~isfinite(cfg.displayFloorDb) || cfg.displayFloorDb >= 0
        error('cfrExperimentRun:InvalidDisplayFloor', 'displayFloorDb must be a finite negative scalar.');
    end
    validateLogicalScalar(cfg.measureInterpolatedPapr, 'measureInterpolatedPapr');
    validateLogicalScalar(cfg.showPlots, 'showPlots');
end

function figureHandle = plotAmplitudeAndSymbolPapr(cfrInfo, cfrCfg)
    figureHandle = figure('Name', 'CFR amplitude and per-symbol PAPR', 'NumberTitle', 'off', 'Color', 'w');
    sampleIndex = (0:length(cfrInfo.worstSymbolInputEnvelope)-1).';
    inputEnvelope = cfrInfo.worstSymbolInputEnvelope / max(cfrInfo.inputRms, realmin('double'));
    outputEnvelope = cfrInfo.worstSymbolOutputEnvelope / max(cfrInfo.outputRms, realmin('double'));
    thresholdRatio = 10^(cfrCfg.targetPaprDb / 20);

    subplot(2, 1, 1);
    inputLine = plot(sampleIndex, inputEnvelope, 'b', 'LineWidth', 1.0);
    hold on;
    outputLine = plot(sampleIndex, outputEnvelope, 'r', 'LineWidth', 1.0);
    thresholdLine = plot([sampleIndex(1) sampleIndex(end)], [thresholdRatio thresholdRatio], 'k--', 'LineWidth', 1.1);
    grid on;
    xlabel('Sample in OFDM interval including CP');
    ylabel('|x| / RMS');
    title(sprintf('Worst original OFDM interval including CP, index %d', cfrInfo.worstSymbolIndex));
    legend([inputLine outputLine thresholdLine], {'CFR OFF', 'CFR ON', sprintf('Target %.1f dB', cfrCfg.targetPaprDb)}, 'Location', 'best');

    subplot(2, 1, 2);
    symbolIndex = (1:length(cfrInfo.inputPaprPerSymbolDb)).';
    inputPaprLine = plot(symbolIndex, cfrInfo.inputPeakRatioPerSymbolDb, 'bo-', 'LineWidth', 1.0, 'MarkerSize', 3);
    hold on;
    outputPaprLine = plot(symbolIndex, cfrInfo.outputPeakRatioPerSymbolDb, 'ro-', 'LineWidth', 1.0, 'MarkerSize', 3);
    targetLine = plot([symbolIndex(1) symbolIndex(end)], [cfrCfg.targetPaprDb cfrCfg.targetPaprDb], 'k--', 'LineWidth', 1.1);
    grid on;
    xlabel('OFDM symbol');
    ylabel('Peak level relative to frame RMS, dB');
    title(sprintf('Peak level by OFDM interval including CP, pulses = %d', cfrInfo.totalPulseCount));
    legend([inputPaprLine outputPaprLine targetLine], {'CFR OFF', 'CFR ON', 'Target'}, 'Location', 'best');
end

function figureHandle = plotCcdfAndWelch(cfrInfo, frequency, psdOriginal, psdCfr, originalPaprDb, cfrPaprDb, cfg)
    figureHandle = figure('Name', 'CFR PAPR CCDF and Welch spectrum', 'NumberTitle', 'off', 'Color', 'w');
    [inputPapr, inputCcdf] = empiricalCcdf(cfrInfo.inputPaprPerSymbolDb);
    [outputPapr, outputCcdf] = empiricalCcdf(cfrInfo.outputPaprPerSymbolDb);

    subplot(2, 1, 1);
    inputLine = semilogy(inputPapr, inputCcdf, 'b', 'LineWidth', 1.3);
    hold on;
    outputLine = semilogy(outputPapr, outputCcdf, 'r', 'LineWidth', 1.3);
    grid on;
    xlabel('OFDM interval PAPR including CP, dB');
    ylabel('Pr(PAPR > x)');
    title(sprintf('Empirical PAPR CCDF over %d OFDM symbols', length(inputPapr)));
    legend([inputLine outputLine], {'CFR OFF', 'CFR ON'}, 'Location', 'best');

    subplot(2, 1, 2);
    commonReference = max([psdOriginal(:); psdCfr(:)]);
    originalPsdDb = 10 * log10(max(psdOriginal / commonReference, realmin('double')));
    cfrPsdDb = 10 * log10(max(psdCfr / commonReference, realmin('double')));
    originalLine = plot(frequency / 1e6, max(originalPsdDb, cfg.displayFloorDb), 'b', 'LineWidth', 1.1);
    hold on;
    cfrLine = plot(frequency / 1e6, max(cfrPsdDb, cfg.displayFloorDb), 'r', 'LineWidth', 1.1);
    grid on;
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    title(sprintf('Welch spectrum after WOLA/CFR, PAPR %.2f -> %.2f dB', originalPaprDb, cfrPaprDb));
    legend([originalLine cfrLine], {'CFR OFF', 'CFR ON'}, 'Location', 'best');
    ylim([cfg.displayFloorDb 5]);
end

function [sortedPapr, ccdf] = empiricalCcdf(paprValues)
    sortedPapr = sort(paprValues(:));
    sampleCount = length(sortedPapr);
    ccdf = (sampleCount:-1:1).' / sampleCount;
end

function [evmPercent, evmDb] = qamEvm(signal, referenceQamSymbols, cfg)
    signal = signal(:);
    receivedNoCp = complex(zeros(cfg.fftSize, cfg.numOfdmSymbols));
    readIndex = 1;
    for symbolIndex = 1:cfg.numOfdmSymbols
        usefulStart = readIndex + cfg.cpLengths(symbolIndex);
        usefulStop = usefulStart + cfg.fftSize - 1;
        receivedNoCp(:, symbolIndex) = signal(usefulStart:usefulStop);
        readIndex = usefulStop + 1;
    end
    receivedGrid = fftshift(fft(receivedNoCp, [], 1) * sqrt(cfg.numActiveSubcarriers) / cfg.fftSize, 1);
    receivedQamSymbols = receivedGrid(cfg.activeIndices, :);
    receivedQamSymbols = receivedQamSymbols(:);
    referenceQamSymbols = referenceQamSymbols(:);
    referenceEnergy = sum(abs(referenceQamSymbols).^2);
    bestLinearGain = sum(conj(referenceQamSymbols) .* receivedQamSymbols) / max(referenceEnergy, realmin('double'));
    if abs(bestLinearGain) <= realmin('double')
        error('cfrExperimentRun:InvalidQamGain', 'The received QAM signal has zero best-fit linear gain.');
    end
    correctedSymbols = receivedQamSymbols / bestLinearGain;
    errorEnergy = sum(abs(correctedSymbols - referenceQamSymbols).^2);
    evmRatio = sqrt(errorEnergy / max(referenceEnergy, realmin('double')));
    evmPercent = 100 * evmRatio;
    evmDb = 20 * log10(max(evmRatio, realmin('double')));
end

function paprDb = signalPaprDb(signal)
    powerValues = abs(signal(:)).^2;
    averagePower = mean(powerValues);
    if averagePower <= 0
        paprDb = 0;
    else
        paprDb = 10 * log10(max(powerValues) / averagePower);
    end
end

function ratioDb = powerRatioDb(numerator, denominator)
    ratioDb = 10 * log10(max(numerator, realmin('double')) / max(denominator, realmin('double')));
end

function validatePositiveInteger(value, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0 || value ~= round(value)
        error('cfrExperimentRun:InvalidPositiveInteger', '%s must be a positive integer.', fieldName);
    end
end

function validateNonnegativeInteger(value, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value < 0 || value ~= round(value)
        error('cfrExperimentRun:InvalidNonnegativeInteger', '%s must be a nonnegative integer.', fieldName);
    end
end

function validateLogicalScalar(value, fieldName)
    if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || (value ~= 0 && value ~= 1)
        error('cfrExperimentRun:InvalidLogicalScalar', '%s must be a logical scalar.', fieldName);
    end
end
