function results = snrRun()
    [txSignal, txInfo] = tx();
    cfg = txInfo.cfg;
    [paSignal, ~] = pa(txSignal, cfg);
    [paFrequency, paPsd] = welch(paSignal, cfg.sampleRate);
    paAclr = aclr(paFrequency, paPsd, cfg.channelBandwidth, cfg.aclrMeasurementBandwidth);

    snrValues = (-5:1:60).';
    numberOfPoints = length(snrValues);

    evmRmsPercent = zeros(numberOfPoints, 1);
    evmRmsDb = zeros(numberOfPoints, 1);
    rxAclr1Db = zeros(numberOfPoints, 1);
    rxAclr2Db = zeros(numberOfPoints, 1);
    referenceEnergy = sum(abs(txInfo.qamSymbols).^2);

    rng(11, 'twister');

    for snrIndex = 1:numberOfPoints
        sweepCfg = cfg;
        sweepCfg.snrDb = snrValues(snrIndex);

        [rxSignal, ~] = channel(paSignal, sweepCfg);
        rxQamSymbols = extractQamSymbols(rxSignal, sweepCfg);

        bestLinearGain = sum( conj(txInfo.qamSymbols) .* rxQamSymbols) / referenceEnergy;
        rxQamSymbols = rxQamSymbols / bestLinearGain;
        symbolError = rxQamSymbols - txInfo.qamSymbols;
        errorEnergy = sum(abs(symbolError).^2);
        evmRms = sqrt(errorEnergy / referenceEnergy);

        evmRmsPercent(snrIndex) = 100 * evmRms;
        evmRmsDb(snrIndex) = 20 * log10(evmRms);

        [frequency, rxPsd] = welch(rxSignal, sweepCfg.sampleRate);

        rxAclr = aclr(frequency, rxPsd, sweepCfg.channelBandwidth, sweepCfg.aclrMeasurementBandwidth);

        rxAclr1Db(snrIndex) = rxAclr.worstAclr1Db;
        rxAclr2Db(snrIndex) = rxAclr.worstAclr2Db;
    end

    txAclr1Db = paAclr.worstAclr1Db * ones(numberOfPoints, 1);

    txAclr2Db = paAclr.worstAclr2Db * ones(numberOfPoints, 1);

    figure('Name', 'EVM / SNR', 'Color', 'w');
    plot(snrValues, evmRmsPercent, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);

    grid on;
    xlabel('SNR, dB');
    ylabel(' EVM, %');
    title(' EVM / SNR');

    figure('Name', 'EVM dB / SNR', 'Color', 'w');
    plot(snrValues, evmRmsDb, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);

    grid on;
    xlabel('SNR, dB');
    ylabel(' EVM, dB');
    title(' EVM dB / SNR');

    figure('Name', 'ACLR / SNR', 'Color', 'w');

    subplot(2, 1, 1);
    plot(snrValues, txAclr1Db, 'b-', 'LineWidth', 1.2);
    hold on;
    plot(snrValues, txAclr2Db, 'r-', 'LineWidth', 1.2);
    grid on;
    xlabel('SNR, dB');
    ylabel(' ACLR, dB');
    title(' ACLR before AWGN');
    legend(' ACLR1BW', ' ACLR2BW', 'Location', 'best');

    subplot(2, 1, 2);
    plot(snrValues, rxAclr1Db, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);
    hold on;
    plot(snrValues, rxAclr2Db, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 5);
    grid on;
    xlabel('SNR, dB');
    ylabel(' ACLR , dB');
    title(' ACLR after AWGN');
    legend(' ACLR1BW', ' ACLR2BW', 'Location', 'northwest');

    results = struct();
    results.snrDb = snrValues;
    results.evmRmsPercent = evmRmsPercent;
    results.evmRmsDb = evmRmsDb;

    results.txAclr1Db = txAclr1Db;
    results.txAclr2Db = txAclr2Db;
    results.rxAclr1Db = rxAclr1Db;
    results.rxAclr2Db = rxAclr2Db;
end

function rxQamSymbols = extractQamSymbols(rxSignal, cfg)
    rxSignal = rxSignal(:);

    rxNoCp = complex(zeros(cfg.fftSize, cfg.numOfdmSymbols));

    readIndex = 1;

    for symbolIndex = 1:cfg.numOfdmSymbols
        cpLength = cfg.cpLengths(symbolIndex);
        usefulStart = readIndex + cpLength;
        usefulStop = usefulStart + cfg.fftSize - 1;

        rxNoCp(:, symbolIndex) = rxSignal(usefulStart:usefulStop);

        readIndex = usefulStop + 1;
    end

    receivedGrid = fft(rxNoCp, [], 1) * sqrt(cfg.numActiveSubcarriers) / cfg.fftSize;

    receivedGrid = fftshift(receivedGrid, 1);

    rxQamMatrix = receivedGrid(cfg.activeIndices, :);

    rxQamSymbols = rxQamMatrix(:);
end
