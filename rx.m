function results = rx(rxSignal, txInfo, channelInfo)
    cfg = txInfo.cfg;
    rxSignal = rxSignal(:);
    symbolLengths = cfg.fftSize + cfg.cpLengths;

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

    referenceQamSymbols = txInfo.qamSymbols;
    referenceEnergy = sum(abs(referenceQamSymbols).^2);
    bestLinearGain = sum( conj(referenceQamSymbols) .* rxQamSymbols) / referenceEnergy;
    rxQamSymbols = rxQamSymbols / bestLinearGain;
    rxBits = demapping(rxQamSymbols);

    symbolError = rxQamSymbols - referenceQamSymbols;
    errorEnergy = sum(abs(symbolError).^2);
    evmRms = sqrt(errorEnergy / referenceEnergy);
    evmRmsPercent = 100 * evmRms;
    evmRmsDb = 20 * log10(evmRms);
    nmse = errorEnergy / referenceEnergy;
    nmseDb = 10 * log10(nmse);

    prePa = txInfo.prePa;
    postPa = txInfo.postPa;

    frequency = postPa.frequency;
    psdPrePa = prePa.psd;
    psdPostPa = postPa.psd;
    [~, psdRx] = welch(rxSignal, cfg.sampleRate);

    prePaAclr = prePa.aclr;
    postPaAclr = postPa.aclr;

    results = struct();
    results.cfg = cfg;
    results.evmRmsPercent = evmRmsPercent;
    results.evmRmsDb = evmRmsDb;
    results.nmse = nmse;
    results.nmseDb = nmseDb;

    results.rxBits = rxBits;
    results.rxQamSymbols = rxQamSymbols;
    results.frequency = frequency;
    results.psdPrePa = psdPrePa;
    results.psdPostPa = psdPostPa;
    results.psdRx = psdRx;
    results.prePaAclr = prePaAclr;
    results.postPaAclr = postPaAclr;
    results.psdTx = psdPostPa;
    results.txAclr = postPaAclr;
    results.channelInfo = channelInfo;
    results.paInfo = txInfo.paInfo;

    fprintf('\nRX metrics:\n');
    fprintf('		  Recovered bits = %d\n', length(rxBits));
    fprintf('		  Corrected RMS EVM = %.3f %%\n', evmRmsPercent);

    fprintf('\t\t  Corrected RMS EVM = %.3f dB\n', evmRmsDb);
    fprintf('\t\t  Corrected NMSE    = %.3f dB\n', nmseDb);

    numSymbolsToPlot = min(2, cfg.numOfdmSymbols);
    samplesToPlot = sum(symbolLengths(1:numSymbolsToPlot));
    timeUs = (0:samplesToPlot-1).' / cfg.sampleRate * 1e6;

    figure('Name', ' waveforms', 'Color', 'w');
    subplot(1, 1, 1);
    plot(timeUs, real(txInfo.txSignal(1:samplesToPlot)), 'b');
    hold on;
    plot(timeUs, real(rxSignal(1:samplesToPlot)), 'r');
    grid on;
    xlabel('time, microseconds');
    ylabel('');
    title('waveform');
    legend('TX after PA', 'RX after AWGN');

    ofdmSignal = txInfo.prePa.signal;
    instantaneousPower = abs(ofdmSignal).^2;
    averageOfdmPower = mean(instantaneousPower);
    powerRelativeDb = 10 * log10( instantaneousPower / averageOfdmPower + eps);
    sortedPowerRelativeDb = sort(powerRelativeDb);
    numberOfPowerSamples = length(sortedPowerRelativeDb);
    powerCdf = (1:numberOfPowerSamples).' /  numberOfPowerSamples;
    paprDb = sortedPowerRelativeDb(end);

    figure('Name', ' CDF', 'Color', 'w');
    stairs(sortedPowerRelativeDb, powerCdf, 'b', 'LineWidth', 1.2);
    grid on;
    xlabel('Instantaneous power relative to average, dB');
    ylabel('Cumulative probability');
    title(sprintf('OFDM  CDF, PAPR = %.2f dB', paprDb));
    ylim([0 1]);

    psdReference = max([psdPrePa; psdPostPa; psdRx]);
    psdPrePaDb = 10 * log10(psdPrePa / psdReference + eps);
    psdPostPaDb = 10 * log10(psdPostPa / psdReference + eps);
    psdRxDb = 10 * log10(psdRx / psdReference + eps);

    figure('Name', 'PSD', 'Color', 'w');
    prePaPsdLine = plot(frequency / 1e6, psdPrePaDb, 'b', 'LineWidth', 1.0);
    hold on;
    postPaPsdLine = plot(frequency / 1e6, psdPostPaDb, 'm', 'LineWidth', 1.0);
    rxPsdLine = plot(frequency / 1e6, psdRxDb, 'r', 'LineWidth', 1.0);
    grid on;
    set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
    xlabel('Frequency, MHz');
    ylabel(' PSD, dB');
    title('PSD');
    legend([prePaPsdLine postPaPsdLine rxPsdLine],  {'before PA', 'after PA',  'tx->pa->channel->rx'}, 'Location', 'northeast');
    xlim([frequency(1) frequency(end)] / 1e6);
    ylim([-80 0]);

    aclrFrequencyMHz = postPaAclr.frequency / 1e6;
    prePaAclrPsdDb = 10 * log10(  prePaAclr.psd / max(prePaAclr.psd) + eps);
    postPaAclrPsdDb = 10 * log10(  postPaAclr.psd / max(postPaAclr.psd) + eps);
    aclrYLimits = [-80 0];
    aclrBandEdgesMHz = postPaAclr.bandEdges / 1e6;
    aclrBandCentersMHz = postPaAclr.bandCenters / 1e6;
    aclrBoundariesMHz = unique(aclrBandEdgesMHz(:));
    aclrLabelCentersMHz = aclrBandCentersMHz([1 2 4 5]);
    prePaValuesDb = [prePaAclr.aclr2LeftDb  prePaAclr.aclr1LeftDb prePaAclr.aclr1RightDb  prePaAclr.aclr2RightDb];
    postPaValuesDb = [postPaAclr.aclr2LeftDb postPaAclr.aclr1LeftDb postPaAclr.aclr1RightDb postPaAclr.aclr2RightDb];
    aclrLabels = {'ACLR -2BW', 'ACLR -BW', 'ACLR BW', 'ACLR 2BW'};

    figure('Name', ' ACLR', 'Color', 'w');
    prePaLine = plot(aclrFrequencyMHz, prePaAclrPsdDb,  'b', 'LineWidth', 1.1);
    hold on;
    postPaLine = plot(aclrFrequencyMHz, postPaAclrPsdDb,  'm', 'LineWidth', 1.1);

    for boundaryIndex = 1:length(aclrBoundariesMHz)
        plot([aclrBoundariesMHz(boundaryIndex)  aclrBoundariesMHz(boundaryIndex)], aclrYLimits, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end

    for labelIndex = 1:4
        text(aclrLabelCentersMHz(labelIndex), aclrYLimits(2)-2, sprintf('%s\n%.2f dB\n%.2f dB',  aclrLabels{labelIndex}, prePaValuesDb(labelIndex), postPaValuesDb(labelIndex)),  'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    end

    grid on;
    xlabel('Frequency, MHz');
    ylabel(' PSD, dB');
    title('ACLR');
    legend([prePaLine postPaLine],  {' before PA', 'after PA'},  'Location', 'northeast');
    xlim([aclrBandEdgesMHz(1, 1)-0.5 aclrBandEdgesMHz(5, 2)+0.5]);
    ylim(aclrYLimits);
    set(gca, 'Layer', 'top');

    pointsToPlot = min(4000, length(rxQamSymbols));
    figure('Name', ' constellation', 'Color', 'w');
    plot(real(rxQamSymbols(1:pointsToPlot)), imag(rxQamSymbols(1:pointsToPlot)), 'r.', 'MarkerSize', 6);
    hold on;
    plot(real(txInfo.qamSymbols(1:pointsToPlot)), imag(txInfo.qamSymbols(1:pointsToPlot)), 'bo', 'MarkerSize', 6, 'LineWidth', 1.2);
    grid on;
    axis equal;
    xlim([-1.6 1.6]);
    ylim([-1.6 1.6]);
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf([' EVM = %.2f %% ' '(%.2f dB), NMSE = %.2f dB'], evmRmsPercent, evmRmsDb, nmseDb));
    legend('RX symbols', 'TX reference');

    amPointsToPlot = min(6000, length(rxSignal));
    amInputComplex = txInfo.prePa.signal(1:amPointsToPlot);
    amOutputComplex = rxSignal(1:amPointsToPlot);
    amInput = abs(amInputComplex);
    amOutput = abs(amOutputComplex);
    amPhase = angle(amOutputComplex .* conj(amInputComplex)) * 180 / pi;
    amMaximum = max([amInput; amOutput]);

    figure('Name', 'AM/AM and AM/PM', 'Color', 'w');
    subplot(1, 2, 1);
    plot(amInput,  amOutput,  'b.', 'LineStyle', 'none', 'MarkerSize', 4);
    hold on;
    plot([0 amMaximum], [0 amMaximum], 'k--', 'LineWidth', 1.0);
    grid on;
    axis equal;
    xlim([0 amMaximum]);
    ylim([0 amMaximum]);
    xlabel(' input amplitude');
    ylabel('output amplitude');
    title('TX-RX AM/AM');

    subplot(1, 2, 2);
    plot(amInput, amPhase, 'r.', 'LineStyle', 'none',  'MarkerSize', 4);
    hold on;
    plot([0 max(amInput)], [0 0], 'k--', 'LineWidth', 1.0);
    grid on;
    xlim([0 max(amInput)]);
    ylim([-180 180]);
    xlabel('input amplitude');
    ylabel('phase shift, degrees');
    title('TX-RX AM/PM');

    paAmPlot(txInfo.paInfo);

    systemRateSignalLength = min(length(txInfo.paInputSignal), length(txInfo.paSignal));
    systemRatePointsToPlot = min(24000, systemRateSignalLength);
    systemRatePointIndices = round(linspace(1, systemRateSignalLength, systemRatePointsToPlot)).';
    systemRateInputAmplitude = abs(txInfo.paInputSignal(systemRatePointIndices));
    systemRateOutputAmplitude = abs(txInfo.paSignal(systemRatePointIndices));
    systemRateMaximum = max([systemRateInputAmplitude; systemRateOutputAmplitude]);

    figure('Name', 'PA AM/AM before interpolation and after decimation', 'Color', 'w');
    plot(systemRateInputAmplitude, systemRateOutputAmplitude, 'b.', 'LineStyle', 'none', 'MarkerSize', 5);
    hold on;
    plot([0 systemRateMaximum], [0 systemRateMaximum], 'k--', 'LineWidth', 1.0);
    grid on;
    axis equal;
    xlim([0 systemRateMaximum]);
    ylim([0 systemRateMaximum]);
    xlabel('input amplitude before interpolation');
    ylabel('output amplitude after decimation');
    title('PA AM/AM before interpolation and after decimation');
end
