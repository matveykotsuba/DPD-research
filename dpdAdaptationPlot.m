function figureHandle = dpdAdaptationPlot(metrics, dpdCfg, deployedAclr)
    if nargin < 3
        deployedAclr = struct();
    end
    figureName = 'DPD adaptation metrics by block';
    existingFigures = findall(0, 'Type', 'figure', 'Name', figureName);
    dpdEnabled = isstruct(dpdCfg) && isfield(dpdCfg, 'enabled') && dpdCfg.enabled;

    if ~dpdEnabled || ~isstruct(metrics) || isempty(fieldnames(metrics))
        if ~isempty(existingFigures)
            close(existingFigures);
        end
        figureHandle = [];
        return;
    end

    requiredFields = {'learningArchitecture', 'updateIndex', 'evmRmsPercent', 'worstAclr1Db', 'worstAclr2Db'};
    for fieldIndex = 1:length(requiredFields)
        if ~isfield(metrics, requiredFields{fieldIndex})
            error('dpdAdaptationPlot:MissingMetricField', 'metrics must contain %s.', requiredFields{fieldIndex});
        end
    end

    if isempty(existingFigures)
        figureHandle = figure('Name', figureName, 'NumberTitle', 'off', 'Color', 'w');
    else
        figureHandle = existingFigures(1);
        if length(existingFigures) > 1
            close(existingFigures(2:end));
        end
        figure(figureHandle);
        clf(figureHandle);
    end

    blockIndex = metrics.updateIndex;
    evmRmsPercent = metrics.evmRmsPercent;
    worstAclr1Db = metrics.worstAclr1Db;
    worstAclr2Db = metrics.worstAclr2Db;
    blockIndex = blockIndex(:);
    evmRmsPercent = evmRmsPercent(:);
    worstAclr1Db = worstAclr1Db(:);
    worstAclr2Db = worstAclr2Db(:);
    metricCount = length(blockIndex);
    if metricCount < 2 || length(evmRmsPercent) ~= metricCount || length(worstAclr1Db) ~= metricCount || length(worstAclr2Db) ~= metricCount
        error('dpdAdaptationPlot:InvalidMetricSize', 'Block indices, ACLR values, and RMS EVM values must have the same length of at least two.');
    end
    if any(~isfinite(blockIndex)) || any(~isfinite(evmRmsPercent)) || any(~isfinite(worstAclr1Db)) || any(~isfinite(worstAclr2Db))
        error('dpdAdaptationPlot:InvalidMetricValue', 'Block indices, ACLR values, and RMS EVM values must be finite.');
    end
    if isfield(dpdCfg, 'learningArchitecture') && ischar(dpdCfg.learningArchitecture)
        architecture = lower(strtrim(dpdCfg.learningArchitecture));
    else
        architecture = lower(strtrim(metrics.learningArchitecture));
    end
    if ~strcmp(architecture, 'direct') && ~strcmp(architecture, 'indirect')
        error('dpdAdaptationPlot:InvalidLearningArchitecture', 'The learning architecture must be direct or indirect.');
    end
    architectureName = upper(architecture);
    diagnosticFields = {'diagnosticUpdateIndex', 'diagnosticFrequencyHz', 'diagnosticPsdDb', 'diagnosticAmAmInputAmplitude', 'diagnosticAmAmOutputAmplitude'};
    diagnosticFieldPresent = false(length(diagnosticFields), 1);
    for fieldIndex = 1:length(diagnosticFields)
        diagnosticFieldPresent(fieldIndex) = isfield(metrics, diagnosticFields{fieldIndex});
    end
    if any(diagnosticFieldPresent) && ~all(diagnosticFieldPresent)
        error('dpdAdaptationPlot:IncompleteDiagnostics', 'Spectrum and AM-AM diagnostic fields must be stored together.');
    end
    hasDiagnostics = all(diagnosticFieldPresent);
    if hasDiagnostics && (~isfield(metrics, 'spectrumEstimator') || ~strcmpi(metrics.spectrumEstimator, 'Welch'))
        error('dpdAdaptationPlot:InvalidSpectrumEstimator', 'Block spectrum diagnostics must use Welch PSD.');
    end
    if hasDiagnostics
        set(figureHandle, 'Position', [100 100 1250 760]);
        aclrAxes = subplot(2, 2, 1);
        evmAxes = subplot(2, 2, 3);
    else
        aclrAxes = subplot(2, 1, 1);
        evmAxes = subplot(2, 1, 2);
    end

    axes(aclrAxes);
    aclr1Line = plot(blockIndex, worstAclr1Db, 'bo-', 'LineWidth', 1.4, 'MarkerSize', 6);
    hold on;
    aclr2Line = plot(blockIndex, worstAclr2Db, 'rs-', 'LineWidth', 1.4, 'MarkerSize', 6);
    aclrLegendHandles = gobjects(5, 1);
    aclrLegendLabels = cell(5, 1);
    aclrLegendHandles(1) = aclr1Line;
    aclrLegendHandles(2) = aclr2Line;
    aclrLegendLabels{1} = 'Worst ACLR1 = min(L,R)';
    aclrLegendLabels{2} = 'Worst ACLR2 = min(L,R)';
    aclrLegendCount = 2;
    if isfield(metrics, 'prePaWorstAclr1Db') && isfinite(metrics.prePaWorstAclr1Db)
        prePaAclr1Line = plot([blockIndex(1) blockIndex(end)], [metrics.prePaWorstAclr1Db metrics.prePaWorstAclr1Db], 'k--', 'LineWidth', 1.2);
        aclrLegendCount = aclrLegendCount + 1;
        aclrLegendHandles(aclrLegendCount) = prePaAclr1Line;
        aclrLegendLabels{aclrLegendCount} = 'ACLR1 before PA';
    end
    targetAclr1Db = NaN;
    if isfield(metrics, 'targetAclr1Db')
        targetAclr1Db = metrics.targetAclr1Db;
    elseif isfield(dpdCfg, 'adaptationTargetAclr1Db')
        targetAclr1Db = dpdCfg.adaptationTargetAclr1Db;
    end
    if isnumeric(targetAclr1Db) && isscalar(targetAclr1Db) && isreal(targetAclr1Db) && isfinite(targetAclr1Db)
        targetAclr1Line = plot([blockIndex(1) blockIndex(end)], [targetAclr1Db targetAclr1Db], 'm:', 'LineWidth', 1.4);
        aclrLegendCount = aclrLegendCount + 1;
        aclrLegendHandles(aclrLegendCount) = targetAclr1Line;
        aclrLegendLabels{aclrLegendCount} = sprintf('ACLR1 target %.1f dB', targetAclr1Db);
    end
    [deployedWorstAclr1Db, deployedAclrAvailable] = deployedWorstAclr1(deployedAclr);
    if deployedAclrAvailable
        deployedDifferenceTolerance = 1e-9 * max(1, abs(deployedWorstAclr1Db));
        if abs(deployedWorstAclr1Db - worstAclr1Db(end)) > deployedDifferenceTolerance
            deployedMarker = plot(blockIndex(end), deployedWorstAclr1Db, 'k.', 'MarkerSize', 20, 'Clipping', 'off');
            aclrLegendCount = aclrLegendCount + 1;
            aclrLegendHandles(aclrLegendCount) = deployedMarker;
            aclrLegendLabels{aclrLegendCount} = sprintf('Active deployed worst ACLR1 = %.2f dB', deployedWorstAclr1Db);
        end
    end
    grid on;
    xlim([blockIndex(1) blockIndex(end)]);
    xticks(blockIndex);
    xlabel('Adaptation block');
    ylabel('ACLR, dB');
    title(sprintf('%s DPD: full-frame Welch ACLR, worst = min(left, right)', architectureName));
    legend(aclrLegendHandles(1:aclrLegendCount), aclrLegendLabels(1:aclrLegendCount), 'Location', 'best');

    axes(evmAxes);
    plot(blockIndex, evmRmsPercent, 'ko-', 'LineWidth', 1.4, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    grid on;
    xlim([blockIndex(1) blockIndex(end)]);
    xticks(blockIndex);
    xlabel('Adaptation block');
    ylabel('RMS EVM, %');
    title(sprintf('%s DPD: QAM EVM on the same full OFDM frame', architectureName));
    if hasDiagnostics
        plotSpectrumDiagnostics(metrics, architectureName);
        plotAmAmDiagnostics(metrics, architectureName);
    end
    drawnow;
end

function [value, available] = deployedWorstAclr1(deployedAclr)
    value = NaN;
    available = false;
    if isempty(deployedAclr)
        return
    end
    if isnumeric(deployedAclr) && isscalar(deployedAclr)
        value = deployedAclr;
    elseif isstruct(deployedAclr) && isfield(deployedAclr, 'worstAclr1Db')
        value = deployedAclr.worstAclr1Db;
    elseif isstruct(deployedAclr) && isfield(deployedAclr, 'aclr1LeftDb') && isfield(deployedAclr, 'aclr1RightDb')
        value = min(deployedAclr.aclr1LeftDb, deployedAclr.aclr1RightDb);
    elseif isstruct(deployedAclr) && isempty(fieldnames(deployedAclr))
        return
    else
        error('dpdAdaptationPlot:InvalidDeployedAclr', 'deployedAclr must contain worstAclr1Db or the left and right ACLR1 values.');
    end
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('dpdAdaptationPlot:InvalidDeployedAclr', 'The deployed worst ACLR1 value must be a finite real scalar.');
    end
    available = true;
end

function plotSpectrumDiagnostics(metrics, architectureName)
    diagnosticUpdateIndex = metrics.diagnosticUpdateIndex(:);
    diagnosticFrequencyHz = metrics.diagnosticFrequencyHz(:);
    diagnosticPsdDb = metrics.diagnosticPsdDb;
    diagnosticCount = length(diagnosticUpdateIndex);
    if diagnosticCount < 2 || size(diagnosticPsdDb, 1) ~= length(diagnosticFrequencyHz) || size(diagnosticPsdDb, 2) ~= diagnosticCount
        error('dpdAdaptationPlot:InvalidSpectrumDiagnosticSize', 'Spectrum diagnostics must contain one PSD column per diagnostic update.');
    end
    if any(~isfinite(diagnosticFrequencyHz)) || any(~isfinite(diagnosticPsdDb(:)))
        error('dpdAdaptationPlot:InvalidSpectrumDiagnosticValue', 'Spectrum diagnostics must contain finite values.');
    end
    subplot(2, 2, 2);
    hold on;
    diagnosticColors = lines(diagnosticCount);
    spectrumLines = gobjects(diagnosticCount, 1);
    legendText = cell(diagnosticCount, 1);
    for diagnosticIndex = 1:diagnosticCount
        spectrumLines(diagnosticIndex) = plot(diagnosticFrequencyHz / 1e6, diagnosticPsdDb(:, diagnosticIndex), 'Color', diagnosticColors(diagnosticIndex, :), 'LineWidth', 1.3);
        legendText{diagnosticIndex} = sprintf('Block %d', diagnosticUpdateIndex(diagnosticIndex));
    end
    grid on;
    xlim([diagnosticFrequencyHz(1) diagnosticFrequencyHz(end)] / 1e6);
    ylim([-100 5]);
    xlabel('Frequency, MHz');
    ylabel('Welch PSD, dB');
    title(sprintf('%s DPD: Welch spectrum by adaptation block', architectureName));
    legend(spectrumLines, legendText, 'Location', 'best');
end

function plotAmAmDiagnostics(metrics, architectureName)
    diagnosticUpdateIndex = metrics.diagnosticUpdateIndex(:);
    inputAmplitude = metrics.diagnosticAmAmInputAmplitude(:);
    outputAmplitude = metrics.diagnosticAmAmOutputAmplitude;
    diagnosticCount = length(diagnosticUpdateIndex);
    if diagnosticCount < 2 || size(outputAmplitude, 1) ~= length(inputAmplitude) || size(outputAmplitude, 2) ~= diagnosticCount
        error('dpdAdaptationPlot:InvalidAmAmDiagnosticSize', 'AM-AM diagnostics must contain one output curve per diagnostic update.');
    end
    if any(~isfinite(inputAmplitude)) || any(~isfinite(outputAmplitude(:))) || any(inputAmplitude < 0) || any(outputAmplitude(:) < 0)
        error('dpdAdaptationPlot:InvalidAmAmDiagnosticValue', 'AM-AM diagnostics must contain finite nonnegative amplitudes.');
    end
    maximumAmplitude = max([inputAmplitude; outputAmplitude(:)]);
    if maximumAmplitude <= realmin
        error('dpdAdaptationPlot:ZeroAmAmDiagnostic', 'AM-AM diagnostics must contain a positive amplitude.');
    end
    subplot(2, 2, 4);
    hold on;
    idealLine = plot([0 maximumAmplitude], [0 maximumAmplitude], 'k--', 'LineWidth', 1.1);
    diagnosticColors = lines(diagnosticCount);
    amAmLines = gobjects(diagnosticCount, 1);
    legendText = cell(diagnosticCount + 1, 1);
    legendText{1} = 'Ideal linear';
    for diagnosticIndex = 1:diagnosticCount
        amAmLines(diagnosticIndex) = plot(inputAmplitude, outputAmplitude(:, diagnosticIndex), 'Color', diagnosticColors(diagnosticIndex, :), 'LineWidth', 1.3);
        legendText{diagnosticIndex + 1} = sprintf('Block %d', diagnosticUpdateIndex(diagnosticIndex));
    end
    grid on;
    axis equal;
    xlim([0 maximumAmplitude]);
    ylim([0 maximumAmplitude]);
    xlabel('Input amplitude');
    ylabel('Output amplitude');
    title(sprintf('%s DPD: DPD+PA AM-AM by adaptation block', architectureName));
    legend([idealLine; amAmLines], legendText, 'Location', 'best');
end
