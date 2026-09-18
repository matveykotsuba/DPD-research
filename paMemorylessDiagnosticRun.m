function results = paMemorylessDiagnosticRun(fileName, snrDb, seed, numberOfOfdmSymbols)
    if nargin < 1 || isempty(fileName)
        fileName = fullfile(fileparts(mfilename('fullpath')), '122M_dpd_775.mat');
    end

    if nargin < 2 || isempty(snrDb)
        snrDb = 25;
    end

    if nargin < 3 || isempty(seed)
        seed = 7;
    end

    if nargin < 4 || isempty(numberOfOfdmSymbols)
        numberOfOfdmSymbols = 14;
    end

    if ~isscalar(numberOfOfdmSymbols) || ~isfinite(numberOfOfdmSymbols) || numberOfOfdmSymbols < 1 || numberOfOfdmSymbols ~= round(numberOfOfdmSymbols)
        error('paMemorylessDiagnosticRun:OfdmSymbolCount', 'numberOfOfdmSymbols must be a positive integer.');
    end

    orders = [1 3 5 7];
    alignmentMemoryOrder = 30;
    lambdaValues = logspace(-12, 2, 300);

    measuredData = paMeasuredCaptures(fileName);
    [fitBasis, fitSampleIndices] = paBuildBasis(measuredData.inputBlock, 'memoryless', orders, 0, 0, alignmentMemoryOrder);
    fitOutputBlocks = measuredData.outputBlocks(fitSampleIndices, :);
    fitScore = paCrossValidate(fitBasis, fitOutputBlocks, lambdaValues, true);
    coefficients = reshape(fitScore.finalCoefficients, 1, []);

    cfg = config();
    cfg.snrDb = snrDb;
    numberOfOfdmSymbols = min(numberOfOfdmSymbols, cfg.numOfdmSymbols);
    cfg.numOfdmSymbols = numberOfOfdmSymbols;
    cfg.cpLengths = cfg.cpLengths(1:numberOfOfdmSymbols);
    cfg.dpd.enabled = false;
    cfg.dpd.adaptationEnabled = false;
    cfg.pa.enabled = true;
    cfg.pa.orders = orders;
    cfg.pa.signalDelays = 0;
    cfg.pa.envelopeDelays = 0;
    cfg.pa.diagonalCount = 0;
    cfg.pa.coefficients = coefficients;
    cfg.pa.referenceInputRms = sqrt(mean(abs(measuredData.inputBlock).^2));
    cfg.pa.maximumInputMagnitude = max(abs(measuredData.inputBlock));
    cfg.pa.modelSampleRate = 737.28e6;
    cfg.pa.systemSampleRate = cfg.sampleRate;
    cfg.pa.rateFactor = cfg.pa.modelSampleRate / cfg.sampleRate;

    [txSignal, txInfo] = tx(cfg, seed, false);
    [paSignal, paInfo] = pa(txSignal, cfg);
    rng(seed + 100, 'twister');
    [rxSignal, channelInfo] = channel(paSignal, cfg);

    referenceSymbols = txInfo.qamSymbols;
    rxSymbolsRaw = extractQamSymbols(rxSignal, cfg);
    [rxSymbolsCorrected, rxGain] = correctCommonGain(referenceSymbols, rxSymbolsRaw);

    rxEvmPercent = calculateEvmPercent(referenceSymbols, rxSymbolsCorrected);
    rxEvmDb = 20 * log10(rxEvmPercent / 100);
    rxPhaseRotationDeg = angle(rxGain) * 180 / pi;

    [frequency, txPsd] = welch(txSignal, cfg.sampleRate);
    [~, paPsd] = welch(paSignal, cfg.sampleRate);
    [~, rxPsd] = welch(rxSignal, cfg.sampleRate);
    psdReference = max([txPsd; paPsd; rxPsd]);
    txPsdDb = 10 * log10(txPsd / psdReference + eps);
    paPsdDb = 10 * log10(paPsd / psdReference + eps);
    rxPsdDb = 10 * log10(rxPsd / psdReference + eps);

    spectrumFigure = figure('Name', 'Memoryless PA spectrum', 'NumberTitle', 'off', 'Color', 'w');
    plot(frequency / 1e6, txPsdDb, 'b', 'LineWidth', 1.0);
    hold on;
    plot(frequency / 1e6, paPsdDb, 'm', 'LineWidth', 1.0);
    plot(frequency / 1e6, rxPsdDb, 'r', 'LineWidth', 1.0);
    grid on;
    xlabel('Frequency, MHz');
    ylabel('Normalized PSD, dB');
    title(sprintf('Memoryless PA, P = 7, SNR = %.1f dB', snrDb));
    legend('TX before PA', 'TX after memoryless PA', 'RX after AWGN', 'Location', 'northeast');
    xlim([frequency(1) frequency(end)] / 1e6);
    ylim([-80 0]);

    pointCount = min(4000, length(referenceSymbols));
    referencePlot = referenceSymbols(1:pointCount);
    constellationFigure = figure('Name', 'Memoryless PA RX constellation', 'NumberTitle', 'off', 'Color', 'w');
    plotConstellation(rxSymbolsCorrected(1:pointCount), referencePlot);
    title(sprintf('RX corrected, SNR %.1f dB, EVM %.2f %%', snrDb, rxEvmPercent));

    amFigure = paAmPlot(paInfo);

    fprintf('\nMemoryless PA diagnostic:\n');
    fprintf('  Orders = [1 3 5 7].\n');
    fprintf('  Tikhonov validation NMSE = %.6f dB.\n', fitScore.validationNmseDb);
    fprintf('  Best lambda = %.9e.\n', fitScore.bestLambda);
    fprintf('  First-order coefficient phase = %.6f degrees.\n', angle(coefficients(1)) * 180 / pi);
    fprintf('  Effective RX constellation rotation = %.6f degrees.\n', rxPhaseRotationDeg);
    fprintf('  Corrected RX EVM at %.1f dB SNR = %.6f %% (%.6f dB).\n', snrDb, rxEvmPercent, rxEvmDb);
    disp('  Tikhonov coefficients:');
    disp(coefficients.');

    results = struct();
    results.fileName = fileName;
    results.orders = orders;
    results.coefficients = coefficients;
    results.bestLambda = fitScore.bestLambda;
    results.validationNmseDb = fitScore.validationNmseDb;
    results.fitScore = fitScore;
    results.cfg = cfg;
    results.txInfo = txInfo;
    results.paInfo = paInfo;
    results.channelInfo = channelInfo;
    results.txSignal = txSignal;
    results.paSignal = paSignal;
    results.rxSignal = rxSignal;
    results.rxGain = rxGain;
    results.rxPhaseRotationDeg = rxPhaseRotationDeg;
    results.rxEvmPercent = rxEvmPercent;
    results.rxEvmDb = rxEvmDb;
    results.snrDb = snrDb;
    results.numberOfOfdmSymbols = numberOfOfdmSymbols;
    results.frequency = frequency;
    results.txPsdDb = txPsdDb;
    results.paPsdDb = paPsdDb;
    results.rxPsdDb = rxPsdDb;
    results.spectrumFigure = spectrumFigure;
    results.constellationFigure = constellationFigure;
    results.amFigure = amFigure;
end

function qamSymbols = extractQamSymbols(signal, cfg)
    signal = signal(:);
    requiredSampleCount = sum(cfg.fftSize + cfg.cpLengths);

    if length(signal) < requiredSampleCount
        error('paMemorylessDiagnosticRun:SignalLength', 'The signal is shorter than one complete OFDM frame.');
    end

    noCp = complex(zeros(cfg.fftSize, cfg.numOfdmSymbols));
    readIndex = 1;

    for symbolIndex = 1:cfg.numOfdmSymbols
        cpLength = cfg.cpLengths(symbolIndex);
        usefulStart = readIndex + cpLength;
        usefulStop = usefulStart + cfg.fftSize - 1;
        noCp(:, symbolIndex) = signal(usefulStart:usefulStop);
        readIndex = usefulStop + 1;
    end

    receivedGrid = fft(noCp, [], 1) * sqrt(cfg.numActiveSubcarriers) / cfg.fftSize;
    receivedGrid = fftshift(receivedGrid, 1);
    qamMatrix = receivedGrid(cfg.activeIndices, :);
    qamSymbols = qamMatrix(:);
end

function [correctedSymbols, commonGain] = correctCommonGain(referenceSymbols, receivedSymbols)
    referenceEnergy = sum(abs(referenceSymbols).^2);
    commonGain = sum(conj(referenceSymbols) .* receivedSymbols) / referenceEnergy;

    if abs(commonGain) == 0
        error('paMemorylessDiagnosticRun:ZeroGain', 'The estimated common complex gain is zero.');
    end

    correctedSymbols = receivedSymbols / commonGain;
end

function evmPercent = calculateEvmPercent(referenceSymbols, receivedSymbols)
    errorEnergy = sum(abs(receivedSymbols - referenceSymbols).^2);
    referenceEnergy = sum(abs(referenceSymbols).^2);
    evmPercent = 100 * sqrt(errorEnergy / referenceEnergy);
end

function plotConstellation(receivedSymbols, referenceSymbols)
    plot(real(receivedSymbols), imag(receivedSymbols), 'r.', 'MarkerSize', 6);
    hold on;
    plot(real(referenceSymbols), imag(referenceSymbols), 'bo', 'MarkerSize', 5, 'LineWidth', 1.0);
    grid on;
    axis equal;
    xlim([-1.6 1.6]);
    ylim([-1.6 1.6]);
    xlabel('In-phase');
    ylabel('Quadrature');
    legend('Output symbols', 'TX reference', 'Location', 'best');
end
