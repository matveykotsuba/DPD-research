function [rxSignal, channelInfo] = channel(txSignal, cfg)

    txSignal = txSignal(:);
    signalPower = mean(abs(txSignal).^2);
    snrLinear = 10^(cfg.snrDb / 10);
    bandwidthRatio = cfg.sampleRate / cfg.occupiedBandwidth;
    targetNoisePower = signalPower / snrLinear * bandwidthRatio;

    noise = sqrt(targetNoisePower / 2) * (randn(size(txSignal)) + 1i * randn(size(txSignal)));

    rxSignal = txSignal + noise;

    channelInfo = struct();
    channelInfo.type = 'awgn';
    channelInfo.snrDb = cfg.snrDb;
    channelInfo.signalPower = signalPower;
    channelInfo.targetNoisePower = targetNoisePower;

    fprintf('Channel: AWGN noise, in-band SNR = %.2f dB.\n', cfg.snrDb);
end
