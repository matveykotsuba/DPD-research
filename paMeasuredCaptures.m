function data = paMeasuredCaptures(fileName)
    loadedData = load(fileName);

    if isfield(loadedData, 'signal_in')
        inputSignal = loadedData.signal_in(:);
        outputSignal = loadedData.signal_out(:);
    else
        inputSignal = loadedData.txSignal(:);
        outputSignal = loadedData.rxSignal(:);
    end

    blockCount = 5;
    sampleCount = min(length(inputSignal), length(outputSignal));
    blockLength = floor(sampleCount / blockCount);
    sampleCount = blockCount * blockLength;
    inputBlocks = reshape(inputSignal(1:sampleCount), blockLength, blockCount);
    outputBlocks = reshape(outputSignal(1:sampleCount), blockLength, blockCount);

    data.fileName = fileName;
    data.blockCount = blockCount;
    data.blockLength = blockLength;
    data.inputBlock = inputBlocks(:, 1);
    data.inputBlocks = inputBlocks;
    data.outputBlocks = outputBlocks;
end
