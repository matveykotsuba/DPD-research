function [txSignal, txInfo] = tx(cfg, seed, calculateSpectrum)
    if nargin < 1 || isempty(cfg)
        cfg = config();
    end

    if nargin < 2 || isempty(seed)
        seed = 7;
    end

    if nargin < 3 || isempty(calculateSpectrum)
        calculateSpectrum = true;
    end

    rng(seed, 'twister');

    numQamSymbols = cfg.numActiveSubcarriers * cfg.numOfdmSymbols;
    numBits = numQamSymbols * cfg.bitsPerQamSymbol;
    bits = randi([0 1], numBits, 1);
    qamSymbols = mapping(bits);

    qamMatrix = reshape(qamSymbols, cfg.numActiveSubcarriers, cfg.numOfdmSymbols);

    shiftedGrid = complex(zeros(cfg.fftSize, cfg.numOfdmSymbols));
    shiftedGrid(cfg.activeIndices, :) = qamMatrix;

    timeNoCp = ifft(ifftshift(shiftedGrid, 1), [], 1);

    timeNoCp = timeNoCp * cfg.fftSize / sqrt(cfg.numActiveSubcarriers);

    if isfield(cfg, 'cfr') && isstruct(cfg.cfr)
        activeCfrCfg = cfg.cfr;
    else
        activeCfrCfg = cfrConfig(false);
    end
    [txSignal, symbolStartIndices, wolaInfo] = wola(timeNoCp, cfg);
    [txSignal, cfrInfo] = cfrCore(txSignal, activeCfrCfg, cfg);
    activeCfrCfg.fftSize = cfg.fftSize;
    activeCfrCfg.activeIndices = cfg.activeIndices;
    activeCfrCfg.cpLengths = cfg.cpLengths;
    activeCfrCfg.wolaLength = cfg.wolaLength;
    cfg.cfr = activeCfrCfg;

    if calculateSpectrum
        [frequency, psd, psdInfo] = welch(txSignal, cfg.sampleRate);
        txAclr = aclr(frequency, psd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);
    else
        frequency = [];
        psd = [];
        psdInfo = struct();
        txAclr = struct();
    end

    txInfo = struct();
    txInfo.cfg = cfg;
    txInfo.seed = seed;
    txInfo.bits = bits;
    txInfo.qamSymbols = qamSymbols;
    txInfo.txSignal = txSignal;
    txInfo.symbolStartIndices = symbolStartIndices;
    txInfo.wola = wolaInfo;
    txInfo.cfr = cfrInfo;
    txInfo.frequency = frequency;
    txInfo.psd = psd;
    txInfo.psdInfo = psdInfo;
    txInfo.aclr = txAclr;
    txInfo.averagePower = mean(abs(txSignal).^2);
    txInfo.peakPower = max(abs(txSignal).^2);
    txInfo.paprDb = 10 * log10(txInfo.peakPower / txInfo.averagePower);
    fprintf('TX:  channel bandwidth = %.1f MHz.\n', cfg.channelBandwidth / 1e6);
    fprintf('     %d RB, %d active subcarriers, 16-QAM.\n', cfg.numResourceBlocks, cfg.numActiveSubcarriers);
    fprintf('     FFT size = %d, Fs = %.2f MHz, deltaF = %.0f kHz.\n', cfg.fftSize, cfg.sampleRate / 1e6, cfg.subcarrierSpacing / 1e3);
    fprintf('     normal CP = %d/%d samples, frame = %d samples (%.1f ms).\n', cfg.normalCpFirst, cfg.normalCpOther, length(txSignal), length(txSignal) / cfg.sampleRate * 1e3);
    fprintf('     OFDM symbols = %d, PAPR = %.2f dB.\n', cfg.numOfdmSymbols, txInfo.paprDb);
    if cfrInfo.enabled
        fprintf('     CFR after CP/WOLA: target %.2f dB, frame PAPR %.2f -> %.2f dB, pulses = %d.\n', cfrInfo.targetPaprDb, cfrInfo.inputPaprDb, cfrInfo.outputPaprDb, cfrInfo.totalPulseCount);
        fprintf('     CFR-only waveform EVM = %.3f %%, TX output PAPR = %.2f dB.\n', cfrInfo.cfrOnlyEvmPercent, txInfo.paprDb);
    else
        fprintf('     CFR after CP/WOLA: disabled (exact bypass).\n');
    end
    fprintf('     WOLA length = %d samples (%.3f microseconds).\n', wolaInfo.length, wolaInfo.duration * 1e6);
    fprintf('     remaining CP = %d/%d samples.\n', wolaInfo.remainingCpFirst, wolaInfo.remainingCpOther);
    if calculateSpectrum
        fprintf('     ACLR measurement: %d records, %.0f Hz bin spacing.\n', psdInfo.recordCount, psdInfo.frequencyStep);
        fprintf('     ACLR1 L/R/worst = %.2f/%.2f/%.2f dB.\n', txAclr.aclr1LeftDb, txAclr.aclr1RightDb, txAclr.worstAclr1Db);
        fprintf('     ACLR2 L/R/worst = %.2f/%.2f/%.2f dB.\n', txAclr.aclr2LeftDb, txAclr.aclr2RightDb, txAclr.worstAclr2Db);
    end
end
