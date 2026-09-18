function [candidateCfg, normalizationInfo] = dpdNormalizeCandidate(modelInput, candidateCfg, paCfg)
    [nominalDpdOutput, ~] = dpdCore(modelInput, candidateCfg);
    firstValidSample = maximumModelDelay(candidateCfg) + maximumModelDelay(paCfg) + 1;
    validRows = firstValidSample:length(modelInput);
    inputPower = mean(abs(modelInput(validRows)).^2);
    dpdPower = mean(abs(nominalDpdOutput(validRows)).^2);
    normalizationSucceeded = isfinite(inputPower) && inputPower > realmin && isfinite(dpdPower) && dpdPower > realmin && all(isfinite(nominalDpdOutput(validRows)));
    if normalizationSucceeded
        coefficientScale = sqrt(inputPower / dpdPower);
        candidateCfg.coefficients = candidateCfg.coefficients * coefficientScale;
    else
        coefficientScale = 1;
    end
    normalizationInfo.reference = 'training signal after interpolation, before validation';
    normalizationInfo.inputAveragePower = inputPower;
    normalizationInfo.nominalDpdAveragePower = dpdPower;
    normalizationInfo.coefficientScale = coefficientScale;
    normalizationInfo.succeeded = logical(normalizationSucceeded);
end

function maximumDelay = maximumModelDelay(modelCfg)
    maximumDelay = max([modelCfg.signalDelays(:); modelCfg.envelopeDelays(:)]);
end
