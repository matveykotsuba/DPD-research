function monitor = dpdMonitoringSignal(systemCfg, dpdCfg, purpose)

    if nargin < 2 || ~isstruct(systemCfg) || ~isstruct(dpdCfg)
        error('dpdMonitoringSignal:MissingConfiguration', 'systemCfg and dpdCfg are required.');
    end
    if nargin < 3 || isempty(purpose)
        purpose = 'learning';
    end
    if ~ischar(purpose)
        error('dpdMonitoringSignal:InvalidPurpose', 'purpose must be a character vector.');
    end
    purpose = lower(strtrim(purpose));
    if strcmp(purpose, 'reporting')
        purpose = 'adaptation';
    end

    switch purpose
        case 'learning'
            defaultSymbolCount = configurationValue(dpdCfg, 'adaptationMetricOfdmSymbols', systemCfg.symbolsPerSlot);
            symbolCount = round(configurationValue(dpdCfg, 'learningMonitorOfdmSymbols', defaultSymbolCount));
            seed = round(configurationValue(dpdCfg, 'learningMonitorSeed', configurationValue(dpdCfg, 'adaptationMetricSeed', 14)));
        case 'adaptation'
            defaultSymbolCount = configurationValue(systemCfg, 'numOfdmSymbols', systemCfg.symbolsPerSlot);
            symbolCount = round(configurationValue(dpdCfg, 'adaptationMetricOfdmSymbols', defaultSymbolCount));
            seed = round(configurationValue(dpdCfg, 'adaptationMetricSeed', 7));
        otherwise
            error('dpdMonitoringSignal:UnknownPurpose', 'Unsupported monitoring purpose: %s.', purpose);
    end
    symbolCount = min(symbolCount, length(systemCfg.cpLengths));
    learningSampleCount = round(configurationValue(dpdCfg, 'learningMonitorSampleCount', 131072));
    stressHeadroomDb = configurationValue(dpdCfg, 'learningMonitorStressHeadroomDb', 0.5);
    validateattributes(symbolCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(seed, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
    validateattributes(learningSampleCount, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(stressHeadroomDb, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});

    monitorCfg = systemCfg;
    monitorCfg.numOfdmSymbols = symbolCount;
    monitorCfg.cpLengths = systemCfg.cpLengths(1:symbolCount);
    [systemInput, txInfo] = tx(monitorCfg, seed, false);
    [modelInput, inputInfo] = paModelInput(systemInput, monitorCfg);
    learningSampleCount = min(learningSampleCount, length(modelInput));
    [fullInputPeak, peakIndex] = max(abs(modelInput));
    firstLearningSample = peakIndex - floor((learningSampleCount - 1) / 2);
    firstLearningSample = max(1, min(firstLearningSample, length(modelInput) - learningSampleCount + 1));
    lastLearningSample = firstLearningSample + learningSampleCount - 1;
    learningInput = modelInput(firstLearningSample:lastLearningSample);
    stressScale = 10^(stressHeadroomDb / 20);
    stressInput = learningInput * stressScale;

    monitor.cfg = monitorCfg;
    monitor.systemInput = systemInput;
    monitor.txInfo = txInfo;
    monitor.modelInput = modelInput;
    monitor.inputInfo = inputInfo;
    monitor.learningInput = learningInput;
    monitor.stressInput = stressInput;
    monitor.info.purpose = purpose;
    monitor.info.referencePoint = 'after interpolation and PA-input scaling, before DPD';
    monitor.info.seed = seed;
    monitor.info.ofdmSymbolCount = symbolCount;
    monitor.info.systemSampleCount = length(systemInput);
    monitor.info.fullModelSampleCount = length(modelInput);
    monitor.info.fullInputPeakMagnitude = fullInputPeak;
    monitor.info.fullInputPeakIndex = peakIndex;
    monitor.info.learningFirstSample = firstLearningSample;
    monitor.info.learningLastSample = lastLearningSample;
    monitor.info.modelSampleCount = length(learningInput);
    monitor.info.learningInputPeakMagnitude = max(abs(learningInput));
    monitor.info.stressHeadroomDb = stressHeadroomDb;
    monitor.info.stressScale = stressScale;
    monitor.info.stressInputPeakMagnitude = max(abs(stressInput));
end

function value = configurationValue(configuration, fieldName, defaultValue)
    if isfield(configuration, fieldName)
        value = configuration.(fieldName);
    else
        value = defaultValue;
    end
end
