function results = cfrEvmSweepRun(userCfg)
    if nargin < 1 || isempty(userCfg)
        userCfg = struct();
    end
    if ~isstruct(userCfg) || ~isscalar(userCfg)
        error('cfrEvmSweepRun:InvalidConfiguration', 'userCfg must be a scalar structure.');
    end

    experimentCfg = defaultExperimentConfiguration();
    experimentCfg = applyUserConfiguration(experimentCfg, userCfg);
    validateExperimentConfiguration(experimentCfg);
    targetPaprDb = sort(experimentCfg.targetPaprDbValues(:));
    if any(diff(targetPaprDb) <= 0)
        error('cfrEvmSweepRun:DuplicateTargets', 'targetPaprDbValues must not contain duplicate values.');
    end

    baseCfg = config('indirect', false);
    experimentCfg.numOfdmSymbols = min(experimentCfg.numOfdmSymbols, length(baseCfg.cpLengths));
    baseCfg.numOfdmSymbols = experimentCfg.numOfdmSymbols;
    baseCfg.cpLengths = baseCfg.cpLengths(1:experimentCfg.numOfdmSymbols);
    baseCfg.cfr.cpLengths = baseCfg.cpLengths;
    baseCfg.cfr.enabled = false;
    if isfield(baseCfg, 'dpd')
        baseCfg.dpd.enabled = false;
        baseCfg.dpd.adaptationEnabled = false;
    end

    fprintf('\nCFR EVM versus clipping-level sweep.\n');
    fprintf('Seed = %d, OFDM symbols = %d, points = %d.\n', experimentCfg.seed, experimentCfg.numOfdmSymbols, length(targetPaprDb));
    fprintf('EVM measurement: after WOLA then CFR, before DPD, PA and channel.\n\n');

    [originalSignal, originalTxInfo] = runTxSilently(baseCfg, experimentCfg.seed);
    [baselineEvmPercent, baselineEvmDb] = qamEvmAfterWola(originalSignal, originalTxInfo.qamSymbols, baseCfg);
    originalUsefulPaprDb = originalTxInfo.cfr.usefulInputPaprDb;
    originalPostWolaPaprDb = signalPaprDb(originalSignal);

    pointCount = length(targetPaprDb);
    evmPercent = zeros(pointCount, 1);
    evmDb = zeros(pointCount, 1);
    coreEvmPercent = zeros(pointCount, 1);
    usefulOutputPaprDb = zeros(pointCount, 1);
    postWolaPaprDb = zeros(pointCount, 1);
    totalPulseCount = zeros(pointCount, 1);
    convergedSymbolCount = zeros(pointCount, 1);
    remainingPeakSymbolCount = zeros(pointCount, 1);
    stalledSymbolCount = zeros(pointCount, 1);
    correctionPowerDb = zeros(pointCount, 1);

    for pointIndex = 1:pointCount
        fprintf('Running point %d/%d: target %.2f dB.\n', pointIndex, pointCount, targetPaprDb(pointIndex));
        pointCfg = baseCfg;
        pointCfg.cfr.enabled = true;
        pointCfg.cfr.targetPaprDb = targetPaprDb(pointIndex);

        [pointSignal, pointTxInfo] = runTxSilently(pointCfg, experimentCfg.seed);
        if ~isequal(originalTxInfo.bits, pointTxInfo.bits) || ~isequal(originalTxInfo.qamSymbols, pointTxInfo.qamSymbols)
            error('cfrEvmSweepRun:ReferenceMismatch', 'Every sweep point must use identical source bits and QAM symbols.');
        end

        [evmPercent(pointIndex), evmDb(pointIndex)] = qamEvmAfterWola(pointSignal, pointTxInfo.qamSymbols, pointCfg);
        coreEvmPercent(pointIndex) = pointTxInfo.cfr.cfrOnlyEvmPercent;
        usefulOutputPaprDb(pointIndex) = pointTxInfo.cfr.usefulOutputPaprDb;
        postWolaPaprDb(pointIndex) = signalPaprDb(pointSignal);
        totalPulseCount(pointIndex) = pointTxInfo.cfr.totalPulseCount;
        convergedSymbolCount(pointIndex) = pointTxInfo.cfr.convergedSymbolCount;
        remainingPeakSymbolCount(pointIndex) = pointTxInfo.cfr.remainingPeakSymbolCount;
        stalledSymbolCount(pointIndex) = pointTxInfo.cfr.stalledSymbolCount;
        correctionPowerDb(pointIndex) = pointTxInfo.cfr.correctionPowerDb;
    end

    fprintf('\n Target    A_th/RMS     QAM EVM    Useful PAPR   Frame PAPR   Pulses   Remaining\n');
    fprintf('  (dB)                    (%%)         (dB)         (dB)                 symbols\n');
    for pointIndex = 1:pointCount
        fprintf('%7.2f %11.4f %11.4f %13.4f %12.4f %9d %11d\n', ...
            targetPaprDb(pointIndex), 10^(targetPaprDb(pointIndex) / 20), evmPercent(pointIndex), ...
            usefulOutputPaprDb(pointIndex), postWolaPaprDb(pointIndex), totalPulseCount(pointIndex), remainingPeakSymbolCount(pointIndex));
    end

    targetReached = remainingPeakSymbolCount == 0;
    actualUsefulPaprReductionDb = originalUsefulPaprDb - usefulOutputPaprDb;
    requestedPeakReductionDb = max(0, originalPostWolaPaprDb - targetPaprDb);

    figureHandle = [];
    if logical(experimentCfg.showPlot)
        figureHandle = plotEvmSweep(targetPaprDb, evmPercent);
    end

    fprintf('\nCFR OFF baseline: useful-symbol PAPR = %.4f dB, post-WOLA PAPR = %.4f dB, QAM EVM = %.6f %%.\n', ...
        originalUsefulPaprDb, originalPostWolaPaprDb, baselineEvmPercent);
    defaultTargetIndex = find(abs(targetPaprDb - baseCfg.cfr.targetPaprDb) <= 1e-12, 1, 'first');
    if ~isempty(defaultTargetIndex)
        fprintf('Current %.1f dB point: QAM EVM = %.4f %%, useful-symbol PAPR = %.4f dB, pulses = %d.\n', ...
            baseCfg.cfr.targetPaprDb, evmPercent(defaultTargetIndex), usefulOutputPaprDb(defaultTargetIndex), totalPulseCount(defaultTargetIndex));
    end

    results.configuration = experimentCfg;
    results.processingPoint = 'QAM EVM after WOLA then CFR, before DPD, PA and channel';
    results.reference.seed = experimentCfg.seed;
    results.reference.numOfdmSymbols = experimentCfg.numOfdmSymbols;
    results.reference.usefulPaprDb = originalUsefulPaprDb;
    results.reference.postWolaPaprDb = originalPostWolaPaprDb;
    results.reference.evmRmsPercent = baselineEvmPercent;
    results.reference.evmRmsDb = baselineEvmDb;
    results.targetPaprDb = targetPaprDb;
    results.thresholdAmplitudeOverRms = 10.^(targetPaprDb / 20);
    results.requestedPeakReductionDb = requestedPeakReductionDb;
    results.evmRmsPercent = evmPercent;
    results.evmRmsDb = evmDb;
    results.evmIncreasePercentPoints = evmPercent - baselineEvmPercent;
    results.coreEvmRmsPercent = coreEvmPercent;
    results.usefulOutputPaprDb = usefulOutputPaprDb;
    results.postWolaPaprDb = postWolaPaprDb;
    results.frameOutputPaprDb = postWolaPaprDb;
    results.actualFramePaprReductionDb = originalPostWolaPaprDb - postWolaPaprDb;
    results.actualUsefulPaprReductionDb = actualUsefulPaprReductionDb;
    results.totalPulseCount = totalPulseCount;
    results.convergedSymbolCount = convergedSymbolCount;
    results.remainingPeakSymbolCount = remainingPeakSymbolCount;
    results.stalledSymbolCount = stalledSymbolCount;
    results.targetReached = targetReached;
    results.correctionPowerDb = correctionPowerDb;
    results.figure = figureHandle;
end

function cfg = defaultExperimentConfiguration()
    cfg.seed = 7;
    cfg.numOfdmSymbols = 140;
    cfg.targetPaprDbValues = (6:0.5:11.5).';
    cfg.showPlot = true;
end

function cfg = applyUserConfiguration(cfg, userCfg)
    userFields = fieldnames(userCfg);
    for fieldIndex = 1:length(userFields)
        fieldName = userFields{fieldIndex};
        if ~isfield(cfg, fieldName)
            error('cfrEvmSweepRun:UnknownConfigurationField', 'Unknown configuration field: %s.', fieldName);
        end
        cfg.(fieldName) = userCfg.(fieldName);
    end
end

function validateExperimentConfiguration(cfg)
    validateNonnegativeInteger(cfg.seed, 'seed');
    validatePositiveInteger(cfg.numOfdmSymbols, 'numOfdmSymbols');
    if ~isnumeric(cfg.targetPaprDbValues) || isempty(cfg.targetPaprDbValues) || ~isvector(cfg.targetPaprDbValues) || ...
            ~isreal(cfg.targetPaprDbValues) || any(~isfinite(cfg.targetPaprDbValues(:))) || any(cfg.targetPaprDbValues(:) <= 0)
        error('cfrEvmSweepRun:InvalidTargets', 'targetPaprDbValues must be a nonempty vector of positive finite values.');
    end
    validateLogicalScalar(cfg.showPlot, 'showPlot');
end

function [signal, txInfo] = runTxSilently(cfg, seed)
    evalc('[signal, txInfo] = tx(cfg, seed, false);');
end

function [evmPercent, evmDb] = qamEvmAfterWola(signal, referenceQamSymbols, cfg)
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
        error('cfrEvmSweepRun:InvalidQamGain', 'The recovered QAM signal has zero best-fit linear gain.');
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

function figureHandle = plotEvmSweep(targetPaprDb, evmPercent)
    figureHandle = figure('Name', 'CFR EVM / clipping level', 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 960 600]);
    evmLine = plot(targetPaprDb, evmPercent, 'bo-', 'LineWidth', 1.4, 'MarkerSize', 6, 'MarkerFaceColor', 'b');

    maximumEvm = max(evmPercent);
    if maximumEvm <= 0
        upperEvmLimit = 1;
    else
        upperEvmLimit = 1.12 * maximumEvm;
    end
    ylim([0 upperEvmLimit]);
    if length(targetPaprDb) == 1
        horizontalPadding = max(0.5, 0.05 * abs(targetPaprDb));
        xlim([targetPaprDb - horizontalPadding targetPaprDb + horizontalPadding]);
    else
        xlim([min(targetPaprDb) max(targetPaprDb)]);
    end

    grid on;
    box on;
    xlabel('Target PAPR, dB');
    ylabel('EVM, %');
    title('Target PAPR/EVM');
    legend(evmLine, {' EVM'}, 'Location', 'northeast');
    set(figureHandle, 'PaperPositionMode', 'auto');
end

function validatePositiveInteger(value, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0 || value ~= round(value)
        error('cfrEvmSweepRun:InvalidPositiveInteger', '%s must be a positive integer.', fieldName);
    end
end

function validateNonnegativeInteger(value, fieldName)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value < 0 || value ~= round(value)
        error('cfrEvmSweepRun:InvalidNonnegativeInteger', '%s must be a nonnegative integer.', fieldName);
    end
end

function validateLogicalScalar(value, fieldName)
    if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || (value ~= 0 && value ~= 1)
        error('cfrEvmSweepRun:InvalidLogicalScalar', '%s must be a logical scalar.', fieldName);
    end
end
