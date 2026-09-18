function result = aclr(frequency, psd, channelBandwidth, measurementBandwidth)
    frequency = frequency(:);
    psd = psd(:);
    frequencyStep = frequency(2) - frequency(1);

    bandCenters = channelBandwidth * [-2 -1 0 1 2].';
    halfMeasurementBandwidth = measurementBandwidth / 2;
    bandEdges = [bandCenters-halfMeasurementBandwidth, bandCenters+halfMeasurementBandwidth];

    binLowerEdges = frequency - frequencyStep / 2;
    binUpperEdges = frequency + frequencyStep / 2;
    bandWeights = zeros(length(frequency), length(bandCenters));

    for bandIndex = 1:length(bandCenters)
        overlapHz = max(0, min(binUpperEdges, bandEdges(bandIndex, 2)) - max(binLowerEdges, bandEdges(bandIndex, 1)));
        bandWeights(:, bandIndex) = overlapHz / frequencyStep;
    end

    weightedPsd = psd .* bandWeights;
    bandPowers = sum(weightedPsd, 1) * frequencyStep;
    mainPower = bandPowers(3);

    aclr1LeftDb = 10 * log10(mainPower / bandPowers(2));
    aclr1RightDb = 10 * log10(mainPower / bandPowers(4));
    aclr2LeftDb = 10 * log10(mainPower / bandPowers(1));
    aclr2RightDb = 10 * log10(mainPower / bandPowers(5));

    result = struct();
    result.aclr1LeftDb = aclr1LeftDb;
    result.aclr1RightDb = aclr1RightDb;
    result.aclr2LeftDb = aclr2LeftDb;
    result.aclr2RightDb = aclr2RightDb;
    result.worstAclr1Db = min(aclr1LeftDb, aclr1RightDb);
    result.worstAclr2Db = min(aclr2LeftDb, aclr2RightDb);
    result.frequency = frequency;
    result.psd = psd;
    result.bandCenters = bandCenters;
    result.bandEdges = bandEdges;
    result.bandPowers = bandPowers;
end
