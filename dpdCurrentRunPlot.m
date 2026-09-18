function figures = dpdCurrentRunPlot(comparison, cfg)
    figureNames = {'DPD spectrum comparison for current main signal'; 'DPD ACLR comparison for current main signal'; 'DPD AM-AM and AM-PM comparison for current main signal'};
    for figureIndex = 1:length(figureNames)
        existingFigures = findall(0, 'Type', 'figure', 'Name', figureNames{figureIndex});
        if ~isempty(existingFigures)
            close(existingFigures);
        end
    end

    if ~isfield(comparison, 'enabled') || ~comparison.enabled
        figures = struct();
        return;
    end

    colors = [0.1500 0.1500 0.1500; 0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560; 0.1804 0.4902 0.1961];
    signalNames = comparison.signalNames;

    spectrumFigure = figure('Name', 'DPD spectrum comparison for current main signal', 'NumberTitle', 'off', 'Color', 'w', 'Visible', 'on');
    hold on;
    spectrumLines = gobjects(4, 1);
    for signalIndex = 1:4
        spectrumLines(signalIndex) = plot(comparison.frequency / 1e6, comparison.psdDb(:, signalIndex), 'Color', colors(signalIndex, :), 'LineWidth', 1.4);
    end
    spectrumHalfSpanHz = min(cfg.pa.modelSampleRate / 2, getComparisonOption(cfg.dpd, 'spectrumHalfSpanHz', 20 * cfg.channelBandwidth));
    grid on;
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    xlim(spectrumHalfSpanHz / 1e6 * [-1 1]);
    ylim([-90 5]);
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    title('Current main signal before channel and AWGN');
    legend(spectrumLines, signalNames, 'Location', 'best');

    aclrFigure = figure('Name', 'DPD ACLR comparison for current main signal', 'NumberTitle', 'off', 'Color', 'w', 'Visible', 'on');
    hold on;
    aclrLines = gobjects(4, 1);
    aclrLegend = cell(4, 1);
    for signalIndex = 1:4
        aclrLines(signalIndex) = plot(comparison.frequency / 1e6, comparison.psdDb(:, signalIndex), 'Color', colors(signalIndex, :), 'LineWidth', 1.4);
        aclrLegend{signalIndex} = sprintf('%s: ACLR1 %.2f dB, ACLR2 %.2f dB', signalNames{signalIndex}, comparison.aclr(signalIndex).worstAclr1Db, comparison.aclr(signalIndex).worstAclr2Db);
    end
    bandEdgesMHz = comparison.aclr(1).bandEdges / 1e6;
    bandCentersMHz = comparison.aclr(1).bandCenters / 1e6;
    bandBoundariesMHz = unique(bandEdgesMHz(:));
    aclrYLimits = [-90 5];
    for boundaryIndex = 1:length(bandBoundariesMHz)
        plot([bandBoundariesMHz(boundaryIndex) bandBoundariesMHz(boundaryIndex)], aclrYLimits, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    bandLabels = {'ACLR2 L', 'ACLR1 L', 'Main', 'ACLR1 R', 'ACLR2 R'};
    for bandIndex = 1:length(bandCentersMHz)
        text(bandCentersMHz(bandIndex), aclrYLimits(2) - 2, bandLabels{bandIndex}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    end
    grid on;
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'Layer', 'top');
    xlim([bandEdgesMHz(1, 1) - 0.5 bandEdgesMHz(end, 2) + 0.5]);
    ylim(aclrYLimits);
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    title('Current main signal ACLR before channel and AWGN');
    legend(aclrLines, aclrLegend, 'Location', 'southoutside');

    amFigure = figure('Name', 'DPD AM-AM and AM-PM comparison for current main signal', 'NumberTitle', 'off', 'Color', 'w', 'Visible', 'on');
    subplot(1, 2, 1);
    hold on;
    amLines = gobjects(4, 1);
    for signalIndex = 1:4
        amLines(signalIndex) = plot(comparison.amInputAmplitude, comparison.amOutputAmplitude(:, signalIndex), 'Color', colors(signalIndex, :), 'LineWidth', 1.4);
    end
    maximumAmplitude = max([comparison.amInputAmplitude; comparison.amOutputAmplitude(:)]);
    grid on;
    axis equal;
    xlim([0 maximumAmplitude]);
    ylim([0 maximumAmplitude]);
    xlabel('Original input amplitude');
    ylabel('Output amplitude');
    title('AM-AM');
    legend(amLines, signalNames, 'Location', 'best');

    subplot(1, 2, 2);
    hold on;
    phaseLines = gobjects(4, 1);
    for signalIndex = 1:4
        phaseLines(signalIndex) = plot(comparison.amInputAmplitude, comparison.amOutputPhaseDeg(:, signalIndex), 'Color', colors(signalIndex, :), 'LineWidth', 1.4);
    end
    grid on;
    xlim([0 max(comparison.amInputAmplitude)]);
    ylim([-180 180]);
    xlabel('Original input amplitude');
    ylabel('Phase relative to Original, degrees');
    title('AM-PM');
    legend(phaseLines, signalNames, 'Location', 'best');

    fprintf('\nDPD comparison for the current main signal:\n');
    disp(comparison.aclrSummary);
    fprintf('PA common-gain NMSE:       %.3f dB\n', comparison.paNmseDb);
    fprintf('DPD+PA common-gain NMSE:   %.3f dB\n', comparison.cascadeNmseDb);
    fprintf('PA output power mismatch:  %+.3f dB\n', comparison.paOutputPowerMismatchDb);

    figures.spectrum = spectrumFigure;
    figures.aclr = aclrFigure;
    figures.am = amFigure;
end

function value = getComparisonOption(dpdCfg, fieldName, defaultValue)
    value = defaultValue;
    if isfield(dpdCfg, 'comparison') && isstruct(dpdCfg.comparison) && isfield(dpdCfg.comparison, fieldName)
        value = dpdCfg.comparison.(fieldName);
    end
end
