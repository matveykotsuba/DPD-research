function results = rateConversionExperimentRun(userCfg)
    if nargin < 1 || isempty(userCfg)
        userCfg = struct();
    end

    projectCfg = config();
    experimentCfg = defaultExperimentConfiguration();
    experimentCfg = applyUserConfiguration(experimentCfg, userCfg);
    validateExperimentConfiguration(experimentCfg, projectCfg);

    projectCfg.numOfdmSymbols = experimentCfg.numOfdmSymbols;
    projectCfg.cpLengths = projectCfg.cpLengths(1:experimentCfg.numOfdmSymbols);

    [txSignal, txInfo] = tx(projectCfg, experimentCfg.seed, false);
    txSignal = txSignal(:);

    inputSampleRate = projectCfg.sampleRate;
    highSampleRate = projectCfg.pa.modelSampleRate;
    rateFactorExact = highSampleRate / inputSampleRate;
    rateFactor = round(rateFactorExact);

    if abs(rateFactorExact - rateFactor) > 1e-12 * max(1, abs(rateFactorExact))
        error('rateConversionExperimentRun:RateFactor', 'The PA model sample rate must be an integer multiple of the TX sample rate.');
    end

    projectCfg.pa.rateFactor = rateFactor;
    [highRateFilter, highRateFilterInfo] = paRateFilter(projectCfg);

    lowRateFilterCfg = projectCfg;
    lowRateFilterCfg.pa.rateFactor = 1;
    lowRateFilterCfg.pa.modelSampleRate = inputSampleRate;
    [lowRateFilter, lowRateFilterInfo] = paRateFilter(lowRateFilterCfg);

    prefilteredSignal = fftFilterSame(txSignal, lowRateFilter);

    interpolatedNoFilter = zeroStuff(txSignal, rateFactor);
    interpolatedPostFilterOnly = paInterpolate(txSignal, rateFactor, highRateFilter);

    filteredBeforeDecimation = fftFilterSame(interpolatedPostFilterOnly, highRateFilter);
    decimatedNoFilter = interpolatedPostFilterOnly(1:rateFactor:end);
    decimatedWithFilter = filteredBeforeDecimation(1:rateFactor:end);

    unfilteredInterpolationFilteredBeforeDecimation = fftFilterSame(interpolatedNoFilter, highRateFilter);
    decimatedFromUnfilteredInterpolationNoFilter = interpolatedNoFilter(1:rateFactor:end);
    decimatedFromUnfilteredInterpolationWithFilter = unfilteredInterpolationFilteredBeforeDecimation(1:rateFactor:end);

    highRateSampleIndex = (0:length(interpolatedPostFilterOnly) - 1).';
    shiftedOfdmCopy = experimentCfg.aliasRelativeAmplitude .* interpolatedPostFilterOnly .* exp(1j * 2 * pi * experimentCfg.aliasOffsetHz .* highRateSampleIndex / highSampleRate);
    aliasingInputHighRate = interpolatedPostFilterOnly + shiftedOfdmCopy;
    aliasingFilteredBeforeDecimation = fftFilterSame(aliasingInputHighRate, highRateFilter);
    aliasingDecimatedNoFilter = aliasingInputHighRate(1:rateFactor:end);
    aliasingDecimatedWithFilter = aliasingFilteredBeforeDecimation(1:rateFactor:end);
    aliasFoldFrequencyHz = foldFrequencyToNyquist(experimentCfg.aliasOffsetHz, inputSampleRate);

    lowWelchLength = min(experimentCfg.welchSegmentLength, length(txSignal));
    highWelchLength = min(experimentCfg.welchSegmentLength * rateFactor, length(interpolatedNoFilter));

    sourceSpectrum = measureSpectrum(txSignal, inputSampleRate, lowWelchLength);
    prefilteredSpectrum = measureSpectrum(prefilteredSignal, inputSampleRate, lowWelchLength);
    interpolatedNoFilterSpectrum = measureSpectrum(interpolatedNoFilter, highSampleRate, highWelchLength);
    interpolatedPostFilterOnlySpectrum = measureSpectrum(interpolatedPostFilterOnly, highSampleRate, highWelchLength);
    decimatedNoFilterSpectrum = measureSpectrum(decimatedNoFilter, inputSampleRate, lowWelchLength);
    decimatedWithFilterSpectrum = measureSpectrum(decimatedWithFilter, inputSampleRate, lowWelchLength);
    decimatedFromUnfilteredInterpolationNoFilterSpectrum = measureSpectrum(decimatedFromUnfilteredInterpolationNoFilter, inputSampleRate, lowWelchLength);
    decimatedFromUnfilteredInterpolationWithFilterSpectrum = measureSpectrum(decimatedFromUnfilteredInterpolationWithFilter, inputSampleRate, lowWelchLength);
    aliasingInputHighRateSpectrum = measureSpectrum(aliasingInputHighRate, highSampleRate, highWelchLength);
    aliasingDecimatedNoFilterSpectrum = measureSpectrum(aliasingDecimatedNoFilter, inputSampleRate, lowWelchLength);
    aliasingDecimatedWithFilterSpectrum = measureSpectrum(aliasingDecimatedWithFilter, inputSampleRate, lowWelchLength);

    sourceReferencePsd = max([sourceSpectrum.psd; prefilteredSpectrum.psd]);
    interpolationReferencePsd = max([interpolatedNoFilterSpectrum.psd; interpolatedPostFilterOnlySpectrum.psd]);
    decimationOutputReferencePsd = max([sourceSpectrum.psd; decimatedNoFilterSpectrum.psd; decimatedWithFilterSpectrum.psd]);
    unfilteredInterpolationDecimationReferencePsd = max([sourceSpectrum.psd; decimatedFromUnfilteredInterpolationNoFilterSpectrum.psd; decimatedFromUnfilteredInterpolationWithFilterSpectrum.psd]);
    aliasingInputReferencePsd = max(aliasingInputHighRateSpectrum.psd);
    aliasingOutputReferencePsd = max([sourceSpectrum.psd; aliasingDecimatedNoFilterSpectrum.psd; aliasingDecimatedWithFilterSpectrum.psd]);

    sourceSpectrum.relativeDb = relativePsdDb(sourceSpectrum.psd, sourceReferencePsd);
    prefilteredSpectrum.relativeDb = relativePsdDb(prefilteredSpectrum.psd, sourceReferencePsd);
    interpolatedNoFilterSpectrum.relativeDb = relativePsdDb(interpolatedNoFilterSpectrum.psd, interpolationReferencePsd);
    interpolatedPostFilterOnlySpectrum.relativeDb = relativePsdDb(interpolatedPostFilterOnlySpectrum.psd, interpolationReferencePsd);
    sourceOutputReferenceSpectrum = sourceSpectrum;
    sourceOutputReferenceSpectrum.relativeDb = relativePsdDb(sourceSpectrum.psd, decimationOutputReferencePsd);
    decimatedNoFilterSpectrum.relativeDb = relativePsdDb(decimatedNoFilterSpectrum.psd, decimationOutputReferencePsd);
    decimatedWithFilterSpectrum.relativeDb = relativePsdDb(decimatedWithFilterSpectrum.psd, decimationOutputReferencePsd);
    sourceUnfilteredInterpolationReferenceSpectrum = sourceSpectrum;
    sourceUnfilteredInterpolationReferenceSpectrum.relativeDb = relativePsdDb(sourceSpectrum.psd, unfilteredInterpolationDecimationReferencePsd);
    decimatedFromUnfilteredInterpolationNoFilterSpectrum.relativeDb = relativePsdDb(decimatedFromUnfilteredInterpolationNoFilterSpectrum.psd, unfilteredInterpolationDecimationReferencePsd);
    decimatedFromUnfilteredInterpolationWithFilterSpectrum.relativeDb = relativePsdDb(decimatedFromUnfilteredInterpolationWithFilterSpectrum.psd, unfilteredInterpolationDecimationReferencePsd);
    aliasingInputHighRateSpectrum.relativeDb = relativePsdDb(aliasingInputHighRateSpectrum.psd, aliasingInputReferencePsd);
    aliasingSourceReferenceSpectrum = sourceSpectrum;
    aliasingSourceReferenceSpectrum.relativeDb = relativePsdDb(sourceSpectrum.psd, aliasingOutputReferencePsd);
    aliasingDecimatedNoFilterSpectrum.relativeDb = relativePsdDb(aliasingDecimatedNoFilterSpectrum.psd, aliasingOutputReferencePsd);
    aliasingDecimatedWithFilterSpectrum.relativeDb = relativePsdDb(aliasingDecimatedWithFilterSpectrum.psd, aliasingOutputReferencePsd);

    filterPassbandHz = projectCfg.pa.resamplerPassbandFrequency;
    sourceFigure = [];
    interpolationFigure = [];
    decimationFigure = [];
    unfilteredInterpolationDecimationFigure = [];
    aliasingInputFigure = [];
    aliasingFigure = [];

    if experimentCfg.showPlots
        sourceFigure = plotSourceExperiment(sourceSpectrum, prefilteredSpectrum, inputSampleRate, filterPassbandHz, experimentCfg);
        interpolationFigure = plotInterpolationExperiment(interpolatedNoFilterSpectrum, interpolatedPostFilterOnlySpectrum, highSampleRate, filterPassbandHz, experimentCfg);
        decimationFigure = plotDecimationExperiment(sourceOutputReferenceSpectrum, decimatedNoFilterSpectrum, decimatedWithFilterSpectrum, inputSampleRate, filterPassbandHz, experimentCfg, 'Decimation after filtered interpolation', 'Filtered interpolation: direct decimation', 'Filtered interpolation: Kaiser FIR before decimation');
        unfilteredInterpolationDecimationFigure = plotDecimationExperiment(sourceUnfilteredInterpolationReferenceSpectrum, decimatedFromUnfilteredInterpolationNoFilterSpectrum, decimatedFromUnfilteredInterpolationWithFilterSpectrum, inputSampleRate, filterPassbandHz, experimentCfg, 'Decimation after zero insertion without filtering', 'Zero insertion x24: direct decimation', 'Zero insertion x24: Kaiser FIR before decimation');
        aliasingInputFigure = plotAliasingInputExperiment(aliasingInputHighRateSpectrum, highSampleRate, inputSampleRate, experimentCfg.aliasOffsetHz, experimentCfg);
        aliasingFigure = plotAliasingExperiment(aliasingSourceReferenceSpectrum, aliasingDecimatedNoFilterSpectrum, aliasingDecimatedWithFilterSpectrum, inputSampleRate, filterPassbandHz, aliasFoldFrequencyHz, rateFactor, experimentCfg);
    end

    nmseNoFilterDb = calculateNmseDb(decimatedNoFilter, txSignal);
    nmseWithFilterDb = calculateNmseDb(decimatedWithFilter, txSignal);
    aliasingNmseNoFilterDb = calculateNmseDb(aliasingDecimatedNoFilter, txSignal);
    aliasingNmseWithFilterDb = calculateNmseDb(aliasingDecimatedWithFilter, txSignal);

    fprintf('\nOFDM rate-conversion experiment:\n');
    fprintf('  OFDM channel/occupied bandwidth = %.3f/%.3f MHz.\n', projectCfg.channelBandwidth / 1e6, projectCfg.occupiedBandwidth / 1e6);
    fprintf('  OFDM symbols = %d, input samples = %d.\n', projectCfg.numOfdmSymbols, length(txSignal));
    fprintf('  Input sample rate = %.3f MHz.\n', inputSampleRate / 1e6);
    fprintf('  Interpolation and decimation factor = %d.\n', rateFactor);
    fprintf('  High sample rate = %.3f MHz.\n', highSampleRate / 1e6);
    fprintf('  Kaiser FIR passband/stopband = %.3f/%.3f MHz.\n', projectCfg.pa.resamplerPassbandFrequency / 1e6, projectCfg.pa.resamplerStopbandFrequency / 1e6);
    fprintf('  Kaiser attenuation = %.1f dB.\n', projectCfg.pa.resamplerAttenuationDb);
    fprintf('  Low-rate FIR order/taps = %d/%d.\n', lowRateFilterInfo.order, length(lowRateFilter));
    fprintf('  High-rate FIR order/taps = %d/%d.\n', highRateFilterInfo.order, length(highRateFilter));
    fprintf('  Low/high FIR group delay = %.6f/%.6f microseconds.\n', lowRateFilterInfo.delaySamples / inputSampleRate * 1e6, highRateFilterInfo.delaySamples / highSampleRate * 1e6);
    fprintf('  Round-trip NMSE after decimation without the Kaiser FIR = %.2f dB.\n', nmseNoFilterDb);
    fprintf('  Round-trip NMSE with the Kaiser FIR before decimation = %.2f dB.\n', nmseWithFilterDb);
    fprintf('  Aliasing test shifted OFDM copy = %+.3f MHz, relative amplitude = %.3f.\n', experimentCfg.aliasOffsetHz / 1e6, experimentCfg.aliasRelativeAmplitude);
    fprintf('  Predicted folded center after decimation = %+.3f MHz.\n', aliasFoldFrequencyHz / 1e6);
    fprintf('  Aliasing-test NMSE without the Kaiser FIR = %.2f dB.\n', aliasingNmseNoFilterDb);
    fprintf('  Aliasing-test NMSE with the Kaiser FIR = %.2f dB.\n', aliasingNmseWithFilterDb);

    results = struct();
    results.configuration = experimentCfg;
    results.txConfiguration = projectCfg;
    results.txInfo = txInfo;
    results.inputSampleRate = inputSampleRate;
    results.highSampleRate = highSampleRate;
    results.outputSampleRate = inputSampleRate;
    results.interpolationFactor = rateFactor;
    results.decimationFactor = rateFactor;
    results.txSignal = txSignal;
    results.prefilteredSignal = prefilteredSignal;
    results.interpolatedNoFilter = interpolatedNoFilter;
    results.interpolatedPostFilterOnly = interpolatedPostFilterOnly;
    results.decimatedNoFilter = decimatedNoFilter;
    results.decimatedWithFilter = decimatedWithFilter;
    results.decimatedFromUnfilteredInterpolationNoFilter = decimatedFromUnfilteredInterpolationNoFilter;
    results.decimatedFromUnfilteredInterpolationWithFilter = decimatedFromUnfilteredInterpolationWithFilter;
    results.shiftedOfdmCopy = shiftedOfdmCopy;
    results.aliasingInputHighRate = aliasingInputHighRate;
    results.aliasingFilteredBeforeDecimation = aliasingFilteredBeforeDecimation;
    results.aliasingDecimatedNoFilter = aliasingDecimatedNoFilter;
    results.aliasingDecimatedWithFilter = aliasingDecimatedWithFilter;
    results.aliasFoldFrequencyHz = aliasFoldFrequencyHz;
    results.lowRateFilter = lowRateFilter;
    results.highRateFilter = highRateFilter;
    results.lowRateFilterInfo = lowRateFilterInfo;
    results.highRateFilterInfo = highRateFilterInfo;
    results.sourceSpectrum = sourceSpectrum;
    results.prefilteredSpectrum = prefilteredSpectrum;
    results.interpolatedNoFilterSpectrum = interpolatedNoFilterSpectrum;
    results.interpolatedPostFilterOnlySpectrum = interpolatedPostFilterOnlySpectrum;
    results.decimatedNoFilterSpectrum = decimatedNoFilterSpectrum;
    results.decimatedWithFilterSpectrum = decimatedWithFilterSpectrum;
    results.decimatedFromUnfilteredInterpolationNoFilterSpectrum = decimatedFromUnfilteredInterpolationNoFilterSpectrum;
    results.decimatedFromUnfilteredInterpolationWithFilterSpectrum = decimatedFromUnfilteredInterpolationWithFilterSpectrum;
    results.aliasingInputHighRateSpectrum = aliasingInputHighRateSpectrum;
    results.aliasingDecimatedNoFilterSpectrum = aliasingDecimatedNoFilterSpectrum;
    results.aliasingDecimatedWithFilterSpectrum = aliasingDecimatedWithFilterSpectrum;
    results.nmseNoFilterDb = nmseNoFilterDb;
    results.nmseWithFilterDb = nmseWithFilterDb;
    results.aliasingNmseNoFilterDb = aliasingNmseNoFilterDb;
    results.aliasingNmseWithFilterDb = aliasingNmseWithFilterDb;
    results.spectrumEstimator = 'Welch';
    results.sourceFigure = sourceFigure;
    results.interpolationFigure = interpolationFigure;
    results.decimationFigure = decimationFigure;
    results.unfilteredInterpolationDecimationFigure = unfilteredInterpolationDecimationFigure;
    results.aliasingInputFigure = aliasingInputFigure;
    results.aliasingFigure = aliasingFigure;
end

function cfg = defaultExperimentConfiguration()
    cfg.numOfdmSymbols = 14;
    cfg.seed = 7;
    cfg.welchSegmentLength = 8192;
    cfg.displayFloorDb = -120;
    cfg.showPlots = true;
    cfg.aliasOffsetHz = 25e6;
    cfg.aliasRelativeAmplitude = 1.0;
end

function cfg = applyUserConfiguration(cfg, userCfg)
    userFields = fieldnames(userCfg);

    for fieldIndex = 1:length(userFields)
        fieldName = userFields{fieldIndex};

        if ~isfield(cfg, fieldName)
            error('rateConversionExperimentRun:UnknownConfigurationField', 'Unknown configuration field: %s.', fieldName);
        end

        cfg.(fieldName) = userCfg.(fieldName);
    end
end

function validateExperimentConfiguration(cfg, projectCfg)
    validatePositiveInteger(cfg.numOfdmSymbols, 'numOfdmSymbols');
    validateNonnegativeInteger(cfg.seed, 'seed');
    validatePositiveInteger(cfg.welchSegmentLength, 'welchSegmentLength');

    if cfg.numOfdmSymbols > length(projectCfg.cpLengths)
        error('rateConversionExperimentRun:OfdmSymbolCount', 'numOfdmSymbols cannot exceed %d for the current TX configuration.', length(projectCfg.cpLengths));
    end

    if projectCfg.pa.modelSampleRate <= projectCfg.sampleRate
        error('rateConversionExperimentRun:SampleRates', 'The PA model sample rate must exceed the TX sample rate.');
    end

    if projectCfg.pa.resamplerPassbandFrequency <= 0 || projectCfg.pa.resamplerStopbandFrequency <= projectCfg.pa.resamplerPassbandFrequency
        error('rateConversionExperimentRun:FilterBand', 'The Kaiser FIR passband and stopband settings are invalid.');
    end

    if projectCfg.pa.resamplerStopbandFrequency > projectCfg.sampleRate / 2
        error('rateConversionExperimentRun:FilterBand', 'The Kaiser FIR stopband must not exceed the low-rate Nyquist frequency.');
    end

    if ~isscalar(cfg.displayFloorDb) || ~isfinite(cfg.displayFloorDb) || cfg.displayFloorDb >= 0
        error('rateConversionExperimentRun:DisplayFloor', 'displayFloorDb must be a finite negative scalar.');
    end

    if ~isscalar(cfg.showPlots)
        error('rateConversionExperimentRun:ShowPlots', 'showPlots must be scalar.');
    end

    if ~isscalar(cfg.aliasOffsetHz) || ~isfinite(cfg.aliasOffsetHz)
        error('rateConversionExperimentRun:AliasOffset', 'aliasOffsetHz must be a finite scalar.');
    end

    if abs(cfg.aliasOffsetHz) - projectCfg.occupiedBandwidth / 2 <= projectCfg.sampleRate / 2
        error('rateConversionExperimentRun:AliasOffset', 'The shifted OFDM copy must lie entirely above the output Nyquist frequency.');
    end

    if abs(cfg.aliasOffsetHz) + projectCfg.occupiedBandwidth / 2 >= projectCfg.pa.modelSampleRate / 2
        error('rateConversionExperimentRun:AliasOffset', 'The shifted OFDM copy must lie inside the high-rate Nyquist interval.');
    end

    if ~isscalar(cfg.aliasRelativeAmplitude) || ~isfinite(cfg.aliasRelativeAmplitude) || cfg.aliasRelativeAmplitude <= 0
        error('rateConversionExperimentRun:AliasAmplitude', 'aliasRelativeAmplitude must be a finite positive scalar.');
    end
end

function validatePositiveInteger(value, name)
    if ~isscalar(value) || ~isfinite(value) || value <= 0 || value ~= round(value)
        error('rateConversionExperimentRun:ConfigurationValue', '%s must be a positive integer.', name);
    end
end

function validateNonnegativeInteger(value, name)
    if ~isscalar(value) || ~isfinite(value) || value < 0 || value ~= round(value)
        error('rateConversionExperimentRun:ConfigurationValue', '%s must be a nonnegative integer.', name);
    end
end

function outputSignal = zeroStuff(inputSignal, factor)
    inputSignal = inputSignal(:);
    outputSignal = complex(zeros(length(inputSignal) * factor, 1));
    outputSignal(1:factor:end) = inputSignal;
end

function spectrum = measureSpectrum(signal, sampleRate, segmentLength)
    [frequency, psd, info] = welch(signal, sampleRate, segmentLength);
    spectrum.frequency = frequency;
    spectrum.psd = psd;
    spectrum.info = info;
end

function relativeDb = relativePsdDb(psd, referencePsd)
    relativeDb = 10 * log10(max(psd / referencePsd, realmin('double')));
end

function foldedFrequencyHz = foldFrequencyToNyquist(frequencyHz, sampleRate)
    foldedFrequencyHz = mod(frequencyHz + sampleRate / 2, sampleRate) - sampleRate / 2;
end

function nmseDb = calculateNmseDb(measuredSignal, referenceSignal)
    measuredSignal = measuredSignal(:);
    referenceSignal = referenceSignal(:);
    sampleCount = min(length(measuredSignal), length(referenceSignal));
    errorSignal = measuredSignal(1:sampleCount) - referenceSignal(1:sampleCount);
    referenceEnergy = sum(abs(referenceSignal(1:sampleCount)).^2);
    nmseDb = 10 * log10(max(sum(abs(errorSignal).^2) / referenceEnergy, realmin('double')));
end

function figureHandle = plotSourceExperiment(sourceSpectrum, prefilteredSpectrum, sampleRate, filterPassbandHz, cfg)
    figureHandle = figure('Name', 'Project OFDM before interpolation', 'NumberTitle', 'off', 'Color', 'w');

    subplot(2, 1, 1);
    plotNormalizedSpectrum(sourceSpectrum, 'b', sampleRate, filterPassbandHz, cfg);
    title('Project OFDM before the low-rate Kaiser FIR');

    subplot(2, 1, 2);
    plotNormalizedSpectrum(prefilteredSpectrum, 'g', sampleRate, filterPassbandHz, cfg);
    title('Project OFDM after the low-rate Kaiser FIR');
end

function figureHandle = plotInterpolationExperiment(noFilterSpectrum, filteredSpectrum, sampleRate, filterPassbandHz, cfg)
    figureHandle = figure('Name', 'Project OFDM interpolation', 'NumberTitle', 'off', 'Color', 'w');

    subplot(2, 1, 1);
    plotNormalizedSpectrum(noFilterSpectrum, 'r', sampleRate, filterPassbandHz, cfg);
    title('Zero insertion x24, no filtering');

    subplot(2, 1, 2);
    plotNormalizedSpectrum(filteredSpectrum, 'g', sampleRate, filterPassbandHz, cfg);
    title('Interpolation x24 with the Kaiser FIR');
end

function figureHandle = plotDecimationExperiment(referenceSpectrum, noFilterSpectrum, filteredSpectrum, outputSampleRate, filterPassbandHz, cfg, figureName, noFilterTitle, filteredTitle)
    figureHandle = figure('Name', figureName, 'NumberTitle', 'off', 'Color', 'w');

    subplot(2, 1, 1);
    referenceLine = plot(referenceSpectrum.frequency / 1e6, max(referenceSpectrum.relativeDb, cfg.displayFloorDb), 'k--', 'LineWidth', 1.0);
    hold on;
    outputLine = plot(noFilterSpectrum.frequency / 1e6, max(noFilterSpectrum.relativeDb, cfg.displayFloorDb), 'r', 'LineWidth', 1.0);
    formatSpectrumAxes(outputSampleRate, filterPassbandHz, cfg);
    title(noFilterTitle);
    legend([referenceLine outputLine], {'Original OFDM', 'Decimated without FIR'}, 'Location', 'best');

    subplot(2, 1, 2);
    referenceLine = plot(referenceSpectrum.frequency / 1e6, max(referenceSpectrum.relativeDb, cfg.displayFloorDb), 'k--', 'LineWidth', 1.0);
    hold on;
    outputLine = plot(filteredSpectrum.frequency / 1e6, max(filteredSpectrum.relativeDb, cfg.displayFloorDb), 'g', 'LineWidth', 1.0);
    formatSpectrumAxes(outputSampleRate, filterPassbandHz, cfg);
    title(filteredTitle);
    legend([referenceLine outputLine], {'Original OFDM', 'Decimated after Kaiser FIR'}, 'Location', 'best');
end

function figureHandle = plotAliasingExperiment(referenceSpectrum, noFilterSpectrum, filteredSpectrum, outputSampleRate, filterPassbandHz, aliasFoldFrequencyHz, rateFactor, cfg)
    figureHandle = figure('Name', 'Aliasing with and without the Kaiser FIR', 'NumberTitle', 'off', 'Color', 'w');

    subplot(2, 1, 1);
    noFilterLine = plot(noFilterSpectrum.frequency / 1e6, max(noFilterSpectrum.relativeDb, cfg.displayFloorDb), 'r', 'LineWidth', 1.0);
    hold on;
    referenceLine = plot(referenceSpectrum.frequency / 1e6, max(referenceSpectrum.relativeDb, cfg.displayFloorDb), 'k--', 'LineWidth', 1.0);
    formatSpectrumAxes(outputSampleRate, filterPassbandHz, cfg);
    plot([aliasFoldFrequencyHz aliasFoldFrequencyHz] / 1e6, [cfg.displayFloorDb 5], 'm:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    title(sprintf('Downsampling x%d without FIR: alias at %+.2f MHz', rateFactor, aliasFoldFrequencyHz / 1e6));
    legend([referenceLine noFilterLine], {'Original OFDM', 'Without FIR before downsampling'}, 'Location', 'best');

    subplot(2, 1, 2);
    filteredLine = plot(filteredSpectrum.frequency / 1e6, max(filteredSpectrum.relativeDb, cfg.displayFloorDb), 'g', 'LineWidth', 1.0);
    hold on;
    referenceLine = plot(referenceSpectrum.frequency / 1e6, max(referenceSpectrum.relativeDb, cfg.displayFloorDb), 'k--', 'LineWidth', 1.0);
    formatSpectrumAxes(outputSampleRate, filterPassbandHz, cfg);
    plot([aliasFoldFrequencyHz aliasFoldFrequencyHz] / 1e6, [cfg.displayFloorDb 5], 'm:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    title(sprintf('Kaiser FIR before downsampling x%d: alias suppressed', rateFactor));
    legend([referenceLine filteredLine], {'Original OFDM', 'Kaiser FIR before downsampling'}, 'Location', 'best');
end

function figureHandle = plotAliasingInputExperiment(highRateSpectrum, highSampleRate, outputSampleRate, aliasOffsetHz, cfg)
    figureHandle = figure('Name', 'Shifted OFDM copy before downsampling', 'NumberTitle', 'off', 'Color', 'w');

    signalLine = plot(highRateSpectrum.frequency / 1e6, max(highRateSpectrum.relativeDb, cfg.displayFloorDb), 'b', 'LineWidth', 1.0);
    [outputNyquistLine, aliasCenterLine] = formatAliasingInputAxes(highSampleRate, outputSampleRate, aliasOffsetHz, cfg);
    title(sprintf('After interpolation FIR, before downsampling: copy at %+.2f MHz', aliasOffsetHz / 1e6));
    legend([signalLine outputNyquistLine aliasCenterLine], {'OFDM and shifted OFDM copy', 'Output Nyquist limits', 'Shifted-copy center'}, 'Location', 'best');
end

function [outputNyquistLine, aliasCenterLine] = formatAliasingInputAxes(highSampleRate, outputSampleRate, aliasOffsetHz, cfg)
    displayLimitHz = min(highSampleRate / 2, max(1.4 * abs(aliasOffsetHz), 1.2 * outputSampleRate / 2));
    grid on;
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    xlim([-displayLimitHz displayLimitHz] / 1e6);
    ylim([cfg.displayFloorDb 5]);
    hold on;
    outputNyquistHz = outputSampleRate / 2;
    outputNyquistLine = plot([outputNyquistHz outputNyquistHz] / 1e6, [cfg.displayFloorDb 5], 'k:', 'LineWidth', 1.0);
    plot([-outputNyquistHz -outputNyquistHz] / 1e6, [cfg.displayFloorDb 5], 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    aliasCenterLine = plot([aliasOffsetHz aliasOffsetHz] / 1e6, [cfg.displayFloorDb 5], 'm:', 'LineWidth', 1.0);
end

function plotNormalizedSpectrum(spectrum, color, sampleRate, filterPassbandHz, cfg)
    plot(spectrum.frequency / 1e6, max(spectrum.relativeDb, cfg.displayFloorDb), color, 'LineWidth', 1.0);
    formatSpectrumAxes(sampleRate, filterPassbandHz, cfg);
end

function formatSpectrumAxes(sampleRate, filterPassbandHz, cfg)
    grid on;
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    xlim([-sampleRate sampleRate] / 2e6);
    ylim([cfg.displayFloorDb 5]);
    hold on;
    plot([filterPassbandHz filterPassbandHz] / 1e6, [cfg.displayFloorDb 5], 'k:', 'HandleVisibility', 'off');
    plot([-filterPassbandHz -filterPassbandHz] / 1e6, [cfg.displayFloorDb 5], 'k:', 'HandleVisibility', 'off');
end
