function results = cfrSelfTest()
    fftSize = 256;
    activeSubcarriers = [-48:-1 1:48].';
    systemCfg.fftSize = fftSize;
    systemCfg.activeIndices = fftSize / 2 + 1 + activeSubcarriers;
    systemCfg.cpLengths = [24; 20; 20; 20];
    systemCfg.numOfdmSymbols = 4;
    systemCfg.wolaLength = 8;
    systemCfg.sampleRate = 256 * 15e3;
    systemCfg.normalCpFirst = 24;
    systemCfg.normalCpOther = 20;

    spectrum = complex(zeros(fftSize, 4));
    constellation = [1+1i; 1-1i; -1+1i; -1-1i] / sqrt(2);
    activeCount = length(systemCfg.activeIndices);
    for symbolIndex = 1:4
        patternIndex = mod((0:activeCount-1).' .* (2 * symbolIndex + 1) + symbolIndex, 4) + 1;
        spectrum(systemCfg.activeIndices, symbolIndex) = constellation(patternIndex);
    end
    spectrum(systemCfg.activeIndices(1:64), 1) = (1 + 1i) / sqrt(2);
    inputSymbols = ifft(ifftshift(spectrum, 1), [], 1) * fftSize / sqrt(length(activeSubcarriers));
    [inputFrame, starts] = wola(inputSymbols, systemCfg);

    disabledCfg = cfrConfig(false);
    [disabledOutput, disabledInfo] = cfrCore(inputFrame, disabledCfg, systemCfg);
    assert(isequal(disabledOutput, inputFrame), 'Disabled CFR must be an exact bypass.');
    assert(~disabledInfo.enabled, 'Disabled CFR info must report enabled=false.');

    enabledCfg = cfrConfig(true);
    enabledCfg.targetPaprDb = 7;
    enabledCfg.maximumPulsesPerSymbol = 64;
    [enabledOutput, enabledInfo] = cfrCore(inputFrame, enabledCfg, systemCfg);
    [repeatOutput, repeatInfo] = cfrCore(inputFrame, enabledCfg, systemCfg);
    assert(isequal(size(enabledOutput), size(inputFrame)), 'CFR must preserve the serialized frame size.');
    assert(all(isfinite(enabledOutput(:))), 'CFR output must be finite.');
    assert(isequal(enabledOutput, repeatOutput), 'CFR must be deterministic.');
    assert(enabledInfo.totalPulseCount == repeatInfo.totalPulseCount, 'Repeated CFR runs must use the same pulse count.');
    assert(enabledInfo.totalPulseCount > 0, 'The coherent-carrier test must trigger peak cancellation.');
    assert(enabledInfo.outputPaprDb < enabledInfo.inputPaprDb, 'Enabled CFR must reduce the deterministic test PAPR.');
    [~, worstPeakLinearIndex] = max(abs(inputFrame));
    expectedWorstSymbolIndex = find(starts <= worstPeakLinearIndex, 1, 'last');
    assert(enabledInfo.worstSymbolIndex == expectedWorstSymbolIndex, 'The reported worst symbol must contain the largest frame-relative input peak.');
    assert(enabledInfo.totalPulseCount <= enabledCfg.maximumPulsesPerSymbol * size(inputSymbols, 2), 'The CFR pulse limit was exceeded.');
    assert(abs(enabledInfo.averagePowerChangeDb) < 1e-10, 'Average power must be preserved.');
    assert(enabledInfo.pulseOutOfBandEnergyRatio < 1e-20, 'The cancellation pulse must be confined to active FFT bins.');

    correctionSpectrum = fftshift(fft(enabledOutput / enabledInfo.normalizationScale - inputFrame));
    frameLength = numel(inputFrame);
    f = (-floor(frameLength/2):ceil(frameLength/2)-1).' / frameLength;
    bins = mod(floor(f*fftSize + 0.5) + floor(fftSize/2), fftSize) + 1;
    activeMask = false(fftSize, 1);
    activeMask(systemCfg.activeIndices) = true;
    outsideMask = ~activeMask(bins);
    outsideEnergy = sum(abs(correctionSpectrum(outsideMask)).^2);
    totalEnergy = sum(abs(correctionSpectrum(:)).^2);
    assert(outsideEnergy / max(totalEnergy, realmin('double')) < 1e-20, 'The CFR correction leaked outside active FFT bins.');

    zeroSymbols = complex(zeros(size(inputFrame)));
    [zeroOutput, zeroInfo] = cfrCore(zeroSymbols, enabledCfg, systemCfg);
    assert(isequal(zeroOutput, zeroSymbols), 'A zero signal must remain zero.');
    assert(zeroInfo.totalPulseCount == 0, 'A zero signal must not create cancellation pulses.');

    projectCfg = config();
    projectCfg.numOfdmSymbols = 2;
    projectCfg.cpLengths = projectCfg.cpLengths(1:2);
    projectCfg.cfr.enabled = false;
    [explicitDisabledSignal, explicitDisabledInfo] = tx(projectCfg, 7, false);
    cfgWithoutCfr = rmfield(projectCfg, 'cfr');
    [legacyDisabledSignal, legacyDisabledInfo] = tx(cfgWithoutCfr, 7, false);
    assert(isequal(explicitDisabledSignal, legacyDisabledSignal), 'A missing CFR configuration must behave as CFR disabled.');
    assert(isequal(explicitDisabledInfo.bits, legacyDisabledInfo.bits), 'CFR compatibility must not change source bits.');
    assert(length(explicitDisabledSignal) == sum(projectCfg.fftSize + projectCfg.cpLengths), 'CFR must not change the WOLA frame length.');
    projectCfg.cfr.enabled = true;
    [txCfrSignal, txCfrInfo] = tx(projectCfg, 7, false);
    [expectedSignal, expectedInfo] = cfrCore(explicitDisabledSignal, projectCfg.cfr, projectCfg);
    assert(isequal(txCfrSignal, expectedSignal), 'TX must apply CFR to the final CP/WOLA signal.');
    assert(isequal(txCfrInfo.symbolStartIndices, explicitDisabledInfo.symbolStartIndices), 'CFR must not shift symbol boundaries.');
    assert(abs(txCfrInfo.paprDb - expectedInfo.outputPaprDb) < 1e-12, 'TX and CFR frame PAPR must agree.');
    assert(isequal(txCfrInfo.bits, explicitDisabledInfo.bits), 'CFR must not alter source bits.');

    cpProbe = complex(zeros(size(inputFrame)));
    cpPeak = starts(2) + 1;
    cpProbe(cpPeak) = 10 + 2i;
    probeCfg = enabledCfg;
    probeCfg.preserveAveragePower = false;
    [probeOut, probeInfo] = cfrCore(cpProbe, probeCfg, systemCfg);
    assert(abs(probeOut(cpPeak)) < abs(cpProbe(cpPeak)), 'A CP peak must be included in cancellation.');
    assert(probeInfo.pulsesPerSymbol(2) > 0, 'CP pulse centres must count towards the corresponding interval.');
    assert(any(abs(probeOut(1:starts(2)-1)) > 0), 'Pulse tails must cross OFDM symbol boundaries.');
    cpProbe(:) = 0;
    cpProbe(1) = 10;
    edgeOut = cfrCore(cpProbe, probeCfg, systemCfg);
    assert(abs(edgeOut(end)) > 0, 'Periodic frame correction must cross the frame boundary.');

    newSignature = txCfrInfo.cfr.signature;
    assert(newSignature.contractVersion == 2, 'Frame CFR needs a new model-compatibility contract.');
    assert(strcmp(newSignature.processingPoint, 'after-cp-wola'), 'Wrong CFR processing point.');
    oldConfiguration = projectCfg.cfr;
    oldConfiguration.algorithmVersion = 1;
    oldConfiguration.pulseShape = 'occupied-subcarrier-rectangular-mask';
    oldSignature = cfrSignature(oldConfiguration);
    assert(~isequaln(oldSignature, newSignature), 'Pre-CP CFR models must not match frame CFR.');

    results.passed = true;
    results.disabledExactBypass = true;
    results.inputPaprDb = enabledInfo.inputPaprDb;
    results.outputPaprDb = enabledInfo.outputPaprDb;
    results.paprReductionDb = enabledInfo.paprReductionDb;
    results.cfrOnlyEvmPercent = enabledInfo.cfrOnlyEvmPercent;
    results.totalPulseCount = enabledInfo.totalPulseCount;
    fprintf('CFR self-test passed. PAPR %.3f -> %.3f dB, pulses = %d.\n', enabledInfo.inputPaprDb, enabledInfo.outputPaprDb, enabledInfo.totalPulseCount);
end
