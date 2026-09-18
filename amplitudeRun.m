function results = amplitudeRun()
    [txSignal, txInfo] = tx();
    cfg = txInfo.cfg;

    inputDriveDb = cfg.pa.sweepInputLevelDb(:);
    paInputRms = cfg.pa.referenceInputRms .* 10.^(inputDriveDb / 20);
    inputBackoffDb = -inputDriveDb;
    numberOfPoints = length(inputDriveDb);

    paOutputRms = zeros(numberOfPoints, 1);
    paOutputPower = zeros(numberOfPoints, 1);
    paOutputPowerDb = zeros(numberOfPoints, 1);
    evmRmsPercent = zeros(numberOfPoints, 1);
    nmseDb = zeros(numberOfPoints, 1);
    aclr1Db = zeros(numberOfPoints, 1);
    aclr2Db = zeros(numberOfPoints, 1);

    referenceEnergy = sum(abs(txInfo.qamSymbols).^2);
    referencePower = cfg.pa.referenceInputRms^2;

    for amplitudeIndex = 1:numberOfPoints
        sweepCfg = cfg;
        sweepCfg.pa.enabled = true;
        sweepCfg.pa.inputBackoffDb = inputBackoffDb(amplitudeIndex);

        [paSignal, paInfo] = pa(txSignal, sweepCfg);

        paOutputPower(amplitudeIndex) = paInfo.modelOutputAveragePower;
        paOutputRms(amplitudeIndex) = sqrt(paOutputPower(amplitudeIndex));
        paOutputPowerDb(amplitudeIndex) = 10 * log10(paOutputPower(amplitudeIndex) / referencePower);

        rxQamSymbols = extractQamSymbols(paSignal, sweepCfg);
        bestLinearGain = sum( conj(txInfo.qamSymbols) .* rxQamSymbols) / referenceEnergy;
        rxQamSymbols = rxQamSymbols / bestLinearGain;
        symbolError = rxQamSymbols - txInfo.qamSymbols;
        errorEnergy = sum(abs(symbolError).^2);

        evmRmsPercent(amplitudeIndex) = 100 * sqrt(errorEnergy / referenceEnergy);
        nmseDb(amplitudeIndex) = 10 * log10(errorEnergy / referenceEnergy);

        [frequency, psd] = welch(paSignal, sweepCfg.sampleRate);
        aclrInfo = aclr(frequency, psd, sweepCfg.channelBandwidth, sweepCfg.aclrMeasurementBandwidth);

        aclr1Db(amplitudeIndex) = aclrInfo.worstAclr1Db;
        aclr2Db(amplitudeIndex) = aclrInfo.worstAclr2Db;
    end

    results = table(inputDriveDb, paInputRms, paOutputRms, paOutputPower, paOutputPowerDb, inputBackoffDb, evmRmsPercent, nmseDb, aclr1Db, aclr2Db, 'VariableNames', {'PAInputLevel_dB', 'PAInputRms', 'PAOutputRms', 'PAOutputPower', 'PAOutputPower_dB', 'InputBackoffDb', 'EVM_percent', 'NMSE_dB', 'ACLR1_dB', 'ACLR2_dB'});

    figure('Name', 'PA output power ', 'Color', 'w');

    subplot(3, 1, 1);
    plot(paOutputPowerDb, evmRmsPercent, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);
    grid on;
    xlabel('PA output power , dB');
    ylabel(' EVM, %');

    subplot(3, 1, 2);
    plot(paOutputPowerDb, nmseDb, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);
    grid on;
    xlabel('PA output power , dB');
    ylabel(' NMSE, dB');

    subplot(3, 1, 3);
    plot(paOutputPowerDb, aclr1Db, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);
    hold on;
    plot(paOutputPowerDb, aclr2Db, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 5);
    grid on;
    xlabel('PA output power , dB');
    ylabel('ACLR, dB');
    legend('ACLR1', 'ACLR2', 'Location', 'best');
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
