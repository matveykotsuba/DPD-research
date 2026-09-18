function results = fftRun(cfg, seed, showPlot)
    if nargin < 1 || isempty(cfg)
        cfg = config();
    end

    if nargin < 2 || isempty(seed)
        seed = 7;
    end

    if nargin < 3 || isempty(showPlot)
        showPlot = true;
    end

    calculateWelchSpectrum = false;
    [txSignal, txInfo] = tx(cfg, seed, calculateWelchSpectrum);
    [rxSignal, channelInfo] = channel(txSignal, cfg);

    sampleCount = min(length(txSignal), length(rxSignal));
    txSignal = txSignal(1:sampleCount);
    rxSignal = rxSignal(1:sampleCount);

    [frequency, txPsd, txWelchInfo] = welch(txSignal, cfg.sampleRate);
    [rxFrequency, rxPsd, rxWelchInfo] = welch(rxSignal, cfg.sampleRate);
    if length(rxFrequency) ~= length(frequency) || any(rxFrequency ~= frequency)
        error('fftRun:InconsistentWelchGrid', 'TX and RX Welch spectra must use the same frequency grid.');
    end

    referencePsd = max([txPsd; rxPsd]);
    txPsdRelativeDb = 10 * log10(max(txPsd / referencePsd, realmin('double')));
    rxPsdRelativeDb = 10 * log10(max(rxPsd / referencePsd, realmin('double')));
    displayFloorDb = -160;

    if showPlot
        figure('Name', 'TX RX Welch spectrum', 'NumberTitle', 'off', 'Color', 'w');
        plot(frequency / 1e6, max(txPsdRelativeDb, displayFloorDb), 'b-', 'LineWidth', 1.0);
        hold on;
        plot(frequency / 1e6, max(rxPsdRelativeDb, displayFloorDb), 'r-', 'LineWidth', 1.0);
        plot([0 0], [displayFloorDb 5], 'k:', 'HandleVisibility', 'off');
        grid on;
        xlabel('Frequency, MHz');
        ylabel('Normalized PSD, dB');
        title('Welch spectrum of TX and RX signals');
        legend('TX', 'RX', 'Location', 'best');
        xlim([frequency(1) frequency(end)] / 1e6);
        ylim([displayFloorDb 5]);
    end

    frequencyStep = txWelchInfo.frequencyStep;
    fftLength = txWelchInfo.fftLength;
    [~, dcIndex] = min(abs(frequency));

    results.frequency = frequency;
    results.txPsd = txPsd;
    results.rxPsd = rxPsd;
    results.txPsdRelativeDb = txPsdRelativeDb;
    results.rxPsdRelativeDb = rxPsdRelativeDb;
    results.referencePsd = referencePsd;
    results.frequencyStep = frequencyStep;
    results.fftLength = fftLength;
    results.sampleCount = sampleCount;
    results.dcIndex = dcIndex;
    results.txDcPsd = txPsd(dcIndex);
    results.rxDcPsd = rxPsd(dcIndex);
    results.txDcRelativeDb = txPsdRelativeDb(dcIndex);
    results.rxDcRelativeDb = rxPsdRelativeDb(dcIndex);
    results.txTimePower = mean(abs(txSignal).^2);
    results.rxTimePower = mean(abs(rxSignal).^2);
    results.txPowerFromPsd = sum(txPsd) * frequencyStep;
    results.rxPowerFromPsd = sum(rxPsd) * frequencyStep;
    results.spectrumEstimator = 'Welch';
    results.txWelchInfo = txWelchInfo;
    results.rxWelchInfo = rxWelchInfo;
    results.txInfo = txInfo;
    results.channelInfo = channelInfo;

    fprintf('Welch spectrum: %d TX records, %d RX records, FFT length %d, frequency step %.0f Hz.\n', txWelchInfo.recordCount, rxWelchInfo.recordCount, fftLength, frequencyStep);
    fprintf('TX DC = %.3e (%.2f dB relative to the common maximum).\n', results.txDcPsd, results.txDcRelativeDb);
    fprintf('RX DC = %.3e (%.2f dB relative to the common maximum).\n', results.rxDcPsd, results.rxDcRelativeDb);
end
