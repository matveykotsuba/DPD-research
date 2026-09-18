function [txSignal, symbolStartIndices, info] = wola(timeNoCp, cfg)
    windowLength = cfg.wolaLength;

    windowIndex = (0:windowLength-1).';
    riseWindow = 0.5 - 0.5 * cos(pi * (windowIndex + 0.5) / windowLength);
    fallWindow = flipud(riseWindow);

    symbolLengths = cfg.fftSize + cfg.cpLengths;
    totalSamples = sum(symbolLengths);

    txSignal = complex(zeros(totalSamples, 1));
    symbolStartIndices = zeros(cfg.numOfdmSymbols, 1);
    writeIndex = 1;

    for symbolIndex = 1:cfg.numOfdmSymbols
        cpLength = cfg.cpLengths(symbolIndex);
        cyclicPrefix = timeNoCp(end-cpLength+1:end, symbolIndex);
        cyclicSuffix = timeNoCp(1:windowLength, symbolIndex);
        extendedSymbol = [cyclicPrefix; timeNoCp(:, symbolIndex); cyclicSuffix];
        symbolWindow = [riseWindow; ones(length(extendedSymbol) - 2 * windowLength, 1); fallWindow];
        windowedSymbol = extendedSymbol .* symbolWindow;

        symbolStartIndices(symbolIndex) = writeIndex;
        outputIndices = (writeIndex:writeIndex+length(windowedSymbol)-1).';
        outputIndices = mod(outputIndices - 1, totalSamples) + 1;
        txSignal(outputIndices) = txSignal(outputIndices) + windowedSymbol;
        writeIndex = writeIndex + symbolLengths(symbolIndex);
    end

    info.length = windowLength;
    info.duration = windowLength / cfg.sampleRate;
    info.remainingCpFirst = cfg.normalCpFirst - windowLength;
    info.remainingCpOther = cfg.normalCpOther - windowLength;
end
