function [figureHandle, comparison] = dpdArchitectureComparisonPlot(dpdCfg, currentMetrics, deployedAclr)
    if nargin < 3
        deployedAclr = struct();
    end
    figureName = 'DPD Direct and Indirect comparison by block';
    existingFigures = findall(0, 'Type', 'figure', 'Name', figureName);
    comparison = struct();
    figureHandle = [];
    enabled = isstruct(dpdCfg) && isfield(dpdCfg, 'enabled') && logical(dpdCfg.enabled);
    plotEnabled = enabled && configurationValue(dpdCfg, 'architectureComparisonPlotEnabled', true);
    if ~plotEnabled
        closeFigures(existingFigures);
        return
    end

    architectures = {'indirect', 'direct'};
    metricsByArchitecture = cell(2, 1);
    structureByArchitecture = cell(2, 1);
    available = false(2, 1);
    architectureStatus = cell(2, 1);
    activeCfrSignature = currentCfrSignature(dpdCfg);
    for architectureIndex = 1:2
        architecture = architectures{architectureIndex};
        modelFile = dpdArchitectureArtifactFile(architecture, activeCfrSignature, 'model');
        [savedMetrics, savedStructure, savedAvailable, savedReason] = savedArchitectureMetrics(modelFile, architecture, dpdCfg);
        [cachedMetrics, cachedStructure, cachedAvailable, cachedReason] = dpdLoadArchitectureComparisonCache(dpdCfg, architecture);
        if cachedAvailable
            [cachedAvailable, cachedValidationReason] = validateMetrics(cachedMetrics, architecture, dpdCfg);
            if ~cachedAvailable
                cachedMetrics = struct();
                cachedStructure = struct();
                cachedReason = ['cached history rejected: ' cachedValidationReason];
            end
        end
        architectureStatus{architectureIndex}.available = false;
        architectureStatus{architectureIndex}.source = 'none';
        architectureStatus{architectureIndex}.reason = [savedReason '; cache: ' cachedReason];
        if cachedAvailable
            metricsByArchitecture{architectureIndex} = cachedMetrics;
            structureByArchitecture{architectureIndex} = cachedStructure;
            available(architectureIndex) = true;
            architectureStatus{architectureIndex}.available = true;
            architectureStatus{architectureIndex}.source = 'cache';
            architectureStatus{architectureIndex}.reason = cachedReason;
        elseif savedAvailable
            metricsByArchitecture{architectureIndex} = savedMetrics;
            structureByArchitecture{architectureIndex} = savedStructure;
            available(architectureIndex) = true;
            architectureStatus{architectureIndex}.available = true;
            architectureStatus{architectureIndex}.source = 'saved';
            architectureStatus{architectureIndex}.reason = savedReason;
        end
    end

    currentMetricsIgnored = false;
    if isstruct(currentMetrics) && ~isempty(fieldnames(currentMetrics)) && isfield(currentMetrics, 'learningArchitecture')
        currentArchitecture = lower(strtrim(currentMetrics.learningArchitecture));
        currentIndex = find(strcmp(architectures, currentArchitecture), 1);
        if ~isempty(currentIndex)
            [currentValid, currentReason] = validateMetrics(currentMetrics, currentArchitecture, dpdCfg);
            if currentValid && ~available(currentIndex)
                metricsByArchitecture{currentIndex} = currentMetrics;
                structureByArchitecture{currentIndex} = dpdStructure(dpdCfg);
                available(currentIndex) = true;
                architectureStatus{currentIndex}.available = true;
                architectureStatus{currentIndex}.source = 'current';
                architectureStatus{currentIndex}.reason = 'compatible current identity-start history loaded';
            elseif ~currentValid
                currentMetricsIgnored = true;
                if available(currentIndex)
                    architectureStatus{currentIndex}.reason = [architectureStatus{currentIndex}.reason '; current history ignored: ' currentReason];
                else
                    architectureStatus{currentIndex}.reason = [architectureStatus{currentIndex}.reason '; current history rejected: ' currentReason];
                end
            end
        end
    end

    if ~all(available)
        closeFigures(existingFigures);
        comparison.available = false;
        comparison.reason = ['one or more compatible identity-start histories are unavailable; Indirect: ' architectureStatus{1}.reason '; Direct: ' architectureStatus{2}.reason];
        comparison.status.indirect = architectureStatus{1};
        comparison.status.direct = architectureStatus{2};
        return
    end

    indirectMetrics = metricsByArchitecture{1};
    directMetrics = metricsByArchitecture{2};
    [comparable, comparisonReason] = comparableMetrics(indirectMetrics, directMetrics, structureByArchitecture{1}, structureByArchitecture{2});
    if ~comparable
        closeFigures(existingFigures);
        comparison.available = false;
        comparison.reason = comparisonReason;
        comparison.status.indirect = architectureStatus{1};
        comparison.status.direct = architectureStatus{2};
        return
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
    set(figureHandle, 'Position', [150 120 1000 720]);

    directBlock = directMetrics.updateIndex(:);
    indirectBlock = indirectMetrics.updateIndex(:);

    subplot(2, 1, 1);
    directAclr1 = plot(directBlock, directMetrics.worstAclr1Db(:), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    hold on;
    directAclr2 = plot(directBlock, directMetrics.worstAclr2Db(:), 'bs--', 'LineWidth', 1.5, 'MarkerSize', 6);
    indirectAclr1 = plot(indirectBlock, indirectMetrics.worstAclr1Db(:), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    indirectAclr2 = plot(indirectBlock, indirectMetrics.worstAclr2Db(:), 'rs--', 'LineWidth', 1.5, 'MarkerSize', 6);
    minimumBlock = min([directBlock; indirectBlock]);
    maximumBlock = max([directBlock; indirectBlock]);
    aclr1ReferenceDb = configurationValue(dpdCfg, 'architectureComparisonAclr1ReferenceDb', 52);
    validateattributes(aclr1ReferenceDb, {'numeric'}, {'scalar', 'real', 'finite'});
    aclr1Reference = plot([minimumBlock maximumBlock], [aclr1ReferenceDb aclr1ReferenceDb], 'k-.', 'LineWidth', 1.5);
    aclrLegendHandles = [directAclr1 directAclr2 indirectAclr1 indirectAclr2 aclr1Reference];
    aclrLegendLabels = {'Direct worst ACLR1', 'Direct worst ACLR2', 'Indirect worst ACLR1', 'Indirect worst ACLR2', sprintf('Pre-PA ACLR1 reference = %.0f dB', aclr1ReferenceDb)};
    [deployedWorstAclr1Db, deployedAclrAvailable] = deployedWorstAclr1(deployedAclr);
    deployedMarkerShown = false;
    activeArchitecture = lower(strtrim(configurationValue(dpdCfg, 'learningArchitecture', '')));
    if deployedAclrAvailable && (strcmp(activeArchitecture, 'direct') || strcmp(activeArchitecture, 'indirect'))
        activeMetrics = metricsByArchitecture{find(strcmp(architectures, activeArchitecture), 1)};
        activeEndpointAclr1Db = activeMetrics.worstAclr1Db(end);
        deployedDifferenceTolerance = 1e-9 * max(1, abs(deployedWorstAclr1Db));
        if abs(deployedWorstAclr1Db - activeEndpointAclr1Db) > deployedDifferenceTolerance
            activeBlock = activeMetrics.updateIndex(:);
            deployedMarker = plot(activeBlock(end), deployedWorstAclr1Db, 'k.', 'MarkerSize', 20, 'Clipping', 'off');
            aclrLegendHandles = [aclrLegendHandles deployedMarker];
            aclrLegendLabels{end + 1} = sprintf('Active %s deployed worst ACLR1 = %.2f dB', upperFirst(activeArchitecture), deployedWorstAclr1Db);
            deployedMarkerShown = true;
        end
    end
    grid on;
    xlim([minimumBlock maximumBlock]);
    xticks(unique([directBlock; indirectBlock]));
    xlabel('Adaptation block');
    ylabel('ACLR, dB');
    title(sprintf('Direct and Indirect DPD: full-frame Welch ACLR (%d symbols), worst = min(L,R)', directMetrics.monitorOfdmSymbolCount));
    legend(aclrLegendHandles, aclrLegendLabels, 'Location', 'best');

    subplot(2, 1, 2);
    directEvm = plot(directBlock, directMetrics.evmRmsPercent(:), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    hold on;
    indirectEvm = plot(indirectBlock, indirectMetrics.evmRmsPercent(:), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    grid on;
    xlim([minimumBlock maximumBlock]);
    xticks(unique([directBlock; indirectBlock]));
    xlabel('Adaptation block');
    ylabel('RMS EVM, %');
    title(sprintf('Direct and Indirect DPD: QAM EVM on the same full %d-symbol OFDM frame', directMetrics.monitorOfdmSymbolCount));
    legend([directEvm indirectEvm], {'Direct', 'Indirect'}, 'Location', 'best');

    comparison.direct = directMetrics;
    comparison.indirect = indirectMetrics;
    comparison.available = true;
    comparison.status.indirect = architectureStatus{1};
    comparison.status.direct = architectureStatus{2};
    if currentMetricsIgnored
        comparison.reason = ['current history ignored: ' currentReason '; comparable saved identity-start histories loaded'];
    else
        comparison.reason = 'comparable identity-start histories loaded';
    end
    comparison.monitorSeed = directMetrics.monitorSeed;
    comparison.monitorOfdmSymbolCount = directMetrics.monitorOfdmSymbolCount;
    comparison.monitorSampleRate = directMetrics.monitorSampleRate;
    comparison.reference = 'same full OFDM frame and Welch ACLR measurement before channel and AWGN';
    comparison.initialPointDefinition = directMetrics.initialPointDefinition;
    comparison.initialPoint.worstAclr1Db = directMetrics.worstAclr1Db(1);
    comparison.initialPoint.worstAclr2Db = directMetrics.worstAclr2Db(1);
    comparison.initialPoint.evmRmsPercent = directMetrics.evmRmsPercent(1);
    comparison.aclr1ReferenceDb = aclr1ReferenceDb;
    comparison.activeArchitecture = activeArchitecture;
    comparison.deployedWorstAclr1Db = deployedWorstAclr1Db;
    comparison.deployedAclrAvailable = deployedAclrAvailable;
    comparison.deployedMarkerShown = deployedMarkerShown;
    drawnow;
end

function [metrics, structure, available, reason] = savedArchitectureMetrics(modelFile, architecture, dpdCfg)
    metrics = struct();
    structure = struct();
    available = false;
    reason = ['model file not found: ' modelFile];
    if exist(modelFile, 'file') ~= 2
        return
    end
    try
        savedModel = load(modelFile, 'fitInfo', 'trainingInfo');
    catch loadError
        reason = ['model file cannot be loaded: ' loadError.message];
        return
    end
    if ~isfield(savedModel, 'fitInfo') || ~isstruct(savedModel.fitInfo) || ~isfield(savedModel.fitInfo, 'learningArchitecture') || ~ischar(savedModel.fitInfo.learningArchitecture)
        reason = 'saved fitInfo or learning architecture is missing or invalid';
        return
    end
    if ~strcmpi(strtrim(savedModel.fitInfo.learningArchitecture), architecture)
        reason = ['saved model architecture does not match ' architecture];
        return
    end
    if ~isfield(savedModel.fitInfo, 'paModelSignature')
        reason = 'saved PA model signature is missing';
        return
    end
    if ~isfield(dpdCfg, 'paModelSignature') || isempty(dpdCfg.paModelSignature)
        reason = 'current PA model signature is unavailable';
        return
    end
    if ~isequaln(savedModel.fitInfo.paModelSignature, dpdCfg.paModelSignature)
        reason = 'saved model was trained for a different PA model or operating point';
        return
    end
    if isfield(savedModel.fitInfo, 'cfrSignature')
        savedCfrSignature = savedModel.fitInfo.cfrSignature;
    else
        savedCfrSignature = cfrSignature(cfrConfig(false));
    end
    if ~isequaln(savedCfrSignature, currentCfrSignature(dpdCfg))
        reason = 'saved model was trained with a different CFR configuration';
        return
    end
    if ~isfield(savedModel, 'trainingInfo') || ~isstruct(savedModel.trainingInfo)
        reason = 'saved trainingInfo is missing or invalid';
        return
    end
    [metrics, available, comparisonMetricsReason] = dpdMakeArchitectureComparisonMetrics(savedModel.trainingInfo);
    if ~available
        reason = comparisonMetricsReason;
        return
    end
    [valid, validationReason] = validateMetrics(metrics, architecture, dpdCfg);
    if ~valid
        available = false;
        metrics = struct();
        reason = validationReason;
        return
    end
    try
        structure = fitInfoStructure(savedModel.fitInfo);
    catch structureError
        available = false;
        metrics = struct();
        structure = struct();
        reason = ['saved DPD structure is invalid: ' structureError.message];
        return
    end
    if ~isequaln(structure, dpdStructure(dpdCfg))
        available = false;
        metrics = struct();
        structure = struct();
        reason = 'saved identity-start history uses a different DPD GMP structure';
        return
    end
    reason = 'compatible saved identity-start history loaded';
end

function [valid, reason] = validateMetrics(metrics, architecture, dpdCfg)
    wrapper.architectureComparisonMetrics = metrics;
    [validatedMetrics, valid, reason] = dpdMakeArchitectureComparisonMetrics(wrapper);
    if ~valid
        return
    end
    if ~strcmpi(strtrim(validatedMetrics.learningArchitecture), architecture)
        valid = false;
        reason = ['metric architecture does not match ' architecture];
        return
    end
    if isfield(dpdCfg, 'adaptationMetricSeed') && validatedMetrics.monitorSeed ~= dpdCfg.adaptationMetricSeed
        valid = false;
        reason = 'metric history uses a different full-frame reporting seed';
        return
    end
    if isfield(dpdCfg, 'adaptationMetricOfdmSymbols') && validatedMetrics.monitorOfdmSymbolCount ~= dpdCfg.adaptationMetricOfdmSymbols
        valid = false;
        reason = 'metric history uses a different full-frame OFDM symbol count';
        return
    end
    if isfield(dpdCfg, 'adaptationMetricSpectrumSegmentLength')
        hasSegmentLength = isfield(validatedMetrics, 'spectrumWelchInfo') && isstruct(validatedMetrics.spectrumWelchInfo) && isscalar(validatedMetrics.spectrumWelchInfo) && isfield(validatedMetrics.spectrumWelchInfo, 'segmentLength');
        if ~hasSegmentLength || validatedMetrics.spectrumWelchInfo.segmentLength ~= dpdCfg.adaptationMetricSpectrumSegmentLength
            valid = false;
            reason = 'metric history uses a different full-frame Welch segment length';
            return
        end
    end
    reason = 'valid identity-start comparison history';
end

function [comparable, reason] = comparableMetrics(indirectMetrics, directMetrics, indirectStructure, directStructure)
    comparable = false;
    if indirectMetrics.monitorSeed ~= directMetrics.monitorSeed || indirectMetrics.monitorOfdmSymbolCount ~= directMetrics.monitorOfdmSymbolCount || indirectMetrics.monitorSampleRate ~= directMetrics.monitorSampleRate
        reason = 'Direct and Indirect histories use different monitoring signals';
        return
    end
    if ~isequaln(indirectStructure, directStructure)
        reason = 'Direct and Indirect histories use different DPD GMP structures';
        return
    end
    if indirectMetrics.trainingSeed ~= directMetrics.trainingSeed || indirectMetrics.trainingSampleCount ~= directMetrics.trainingSampleCount || indirectMetrics.trainingBlockLength ~= directMetrics.trainingBlockLength
        reason = 'Direct and Indirect histories use different training protocols';
        return
    end
    if ~isequal(indirectMetrics.updateIndex(:), directMetrics.updateIndex(:))
        reason = 'Direct and Indirect histories use different adaptation block indices';
        return
    end
    if ~strcmp(indirectMetrics.evmReference, directMetrics.evmReference) || ~strcmp(indirectMetrics.aclrReference, directMetrics.aclrReference) || ~strcmpi(indirectMetrics.spectrumEstimator, directMetrics.spectrumEstimator)
        reason = 'Direct and Indirect histories use different metric definitions';
        return
    end
    initialTolerance = 1e-9;
    initialDifferences = [indirectMetrics.worstAclr1Db(1) - directMetrics.worstAclr1Db(1), indirectMetrics.worstAclr2Db(1) - directMetrics.worstAclr2Db(1), indirectMetrics.evmRmsPercent(1) - directMetrics.evmRmsPercent(1)];
    if any(abs(initialDifferences) > initialTolerance)
        reason = 'Direct and Indirect histories do not share the same identity-start point';
        return
    end
    comparable = true;
    reason = 'comparable identity-start histories loaded';
end

function structure = fitInfoStructure(fitInfo)
    requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount'};
    for fieldIndex = 1:length(requiredFields)
        fieldName = requiredFields{fieldIndex};
        if ~isfield(fitInfo, fieldName)
            error('dpdArchitectureComparisonPlot:MissingStructureField', 'Saved fitInfo must contain %s.', fieldName);
        end
        structure.(fieldName) = fitInfo.(fieldName);
    end
end

function structure = dpdStructure(dpdCfg)
    requiredFields = {'orders', 'signalDelays', 'envelopeDelays', 'diagonalCount'};
    for fieldIndex = 1:length(requiredFields)
        fieldName = requiredFields{fieldIndex};
        if ~isfield(dpdCfg, fieldName)
            error('dpdArchitectureComparisonPlot:MissingConfigurationField', 'dpdCfg must contain %s.', fieldName);
        end
        structure.(fieldName) = dpdCfg.(fieldName);
    end
end

function closeFigures(figures)
    if ~isempty(figures)
        close(figures);
    end
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
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
        error('dpdArchitectureComparisonPlot:InvalidDeployedAclr', 'deployedAclr must contain worstAclr1Db or the left and right ACLR1 values.');
    end
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('dpdArchitectureComparisonPlot:InvalidDeployedAclr', 'The deployed worst ACLR1 value must be a finite real scalar.');
    end
    available = true;
end

function textValue = upperFirst(textValue)
    textValue(1) = upper(textValue(1));
end

function signature = currentCfrSignature(dpdCfg)
    if isfield(dpdCfg, 'cfrSignature') && ~isempty(dpdCfg.cfrSignature)
        signature = dpdCfg.cfrSignature;
    else
        signature = cfrSignature(cfrConfig(false));
    end
end
