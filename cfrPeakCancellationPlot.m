function figureHandle = cfrPeakCancellationPlot(cfrInfo, cfrCfg)
    if nargin < 2
        error('cfrPeakCancellationPlot:MissingInput', 'cfrInfo and cfrCfg are required.');
    end
    if ~isstruct(cfrInfo) || ~isscalar(cfrInfo) || ~isfield(cfrInfo, 'enabled') || ~cfrInfo.enabled
        error('cfrPeakCancellationPlot:InvalidCfrInfo', 'cfrInfo must describe an enabled CFR run.');
    end

    figureHandle = figure('Name', 'CFR peak detection and cancellation', 'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1000 600]);

    inputRms = max(cfrInfo.inputRms, realmin('double'));
    sampleIndex = cfrInfo.peakWindowSampleIndices;
    inputEnvelope = cfrInfo.peakWindowInputEnvelope / inputRms;
    outputEnvelope = cfrInfo.peakWindowOutputEnvelope / inputRms;
    thresholdRatio = cfrInfo.thresholdMagnitude / inputRms;

    previousEnvelope = inputEnvelope([1 1:end-1]);
    nextEnvelope = inputEnvelope([2:end end]);
    thresholdWithTolerance = thresholdRatio * (1 + cfrInfo.signature.relativeTolerance);
    peakMask = inputEnvelope > thresholdWithTolerance & inputEnvelope >= previousEnvelope & inputEnvelope > nextEnvelope;
    peakIndices = find(peakMask);
    if isempty(peakIndices)
        [largestEnvelope, largestIndex] = max(inputEnvelope);
        if largestEnvelope > thresholdRatio
            peakIndices = largestIndex;
        end
    end

    [~, strongestPeakIndex] = max(inputEnvelope);
    displayWindowLength = min(512, length(inputEnvelope));
    firstDisplayedIndex = strongestPeakIndex - floor(displayWindowLength / 2);
    firstDisplayedIndex = max(1, min(firstDisplayedIndex, length(inputEnvelope) - displayWindowLength + 1));
    displayedIndices = (firstDisplayedIndex:firstDisplayedIndex + displayWindowLength - 1).';
    displayedPeakIndices = peakIndices(peakIndices >= displayedIndices(1) & peakIndices <= displayedIndices(end));

    aboveThresholdEnvelope = inputEnvelope;
    aboveThresholdEnvelope(inputEnvelope <= thresholdRatio) = NaN;

    outputLine = plot(sampleIndex(displayedIndices), outputEnvelope(displayedIndices), 'r--', 'LineWidth', 1.1);
    hold on;
    inputLine = plot(sampleIndex(displayedIndices), inputEnvelope(displayedIndices), 'b', 'LineWidth', 1.0);
    thresholdLine = plot([sampleIndex(displayedIndices(1)) sampleIndex(displayedIndices(end))], [thresholdRatio thresholdRatio], 'k--', 'LineWidth', 1.1);
    aboveThresholdLine = plot(sampleIndex(displayedIndices), aboveThresholdEnvelope(displayedIndices), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.8);

    if ~isempty(displayedPeakIndices)
        stemX = reshape([sampleIndex(displayedPeakIndices).'; sampleIndex(displayedPeakIndices).'; nan(1, length(displayedPeakIndices))], [], 1);
        stemY = reshape([thresholdRatio * ones(1, length(displayedPeakIndices)); inputEnvelope(displayedPeakIndices).'; nan(1, length(displayedPeakIndices))], [], 1);
        plot(stemX, stemY, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.0);
        peakLine = plot(sampleIndex(displayedPeakIndices), inputEnvelope(displayedPeakIndices), 'o', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', 'w', 'LineWidth', 1.2, 'MarkerSize', 6);
    else
        peakLine = plot(NaN, NaN, 'o', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', 'w', 'LineWidth', 1.2, 'MarkerSize', 6);
    end

    grid on;
    xlabel('Time samples');
    ylabel('|x| /  RMS');
    title(sprintf('CFR '));
    legend([inputLine outputLine thresholdLine aboveThresholdLine peakLine], ...
        {'CFR OFF', 'CFR ON', sprintf('Target %.1f dB', cfrCfg.targetPaprDb), 'Above target', 'Local peaks'}, ...
        'Location', 'northeast');
    xlim([sampleIndex(displayedIndices(1)) sampleIndex(displayedIndices(end))]);
    maximumDisplayedEnvelope = max([inputEnvelope(displayedIndices); outputEnvelope(displayedIndices); thresholdRatio]);
    ylim([0 1.08 * maximumDisplayedEnvelope]);
    set(figureHandle, 'PaperPositionMode', 'auto');
end
