function cfg = config(dpdLearningArchitecture, cfrEnabled)
    if nargin < 1
        dpdLearningArchitecture = [];
    end

    if nargin < 2
        cfrEnabled = [];
    end

    cfg.channelBandwidth = 5e6;
    cfg.numResourceBlocks = 25;
    cfg.subcarriersPerResourceBlock = 12;
    cfg.numActiveSubcarriers = cfg.numResourceBlocks * cfg.subcarriersPerResourceBlock;

    cfg.subcarrierSpacing = 15e3;
    cfg.fftSize = 2048;
    cfg.sampleRate = cfg.fftSize * cfg.subcarrierSpacing;
    cfg.samplePeriod = 1 / cfg.sampleRate;
    cfg.occupiedBandwidth = cfg.numActiveSubcarriers * cfg.subcarrierSpacing;

    cfg.symbolsPerSlot = 7;
    cfg.slotsPerSubframe = 2;
    cfg.subframesPerFrame = 10;
    cfg.slotsPerFrame = cfg.slotsPerSubframe * cfg.subframesPerFrame;
    cfg.numOfdmSymbols = cfg.symbolsPerSlot * cfg.slotsPerFrame;

    cfg.normalCpFirst = 160;
    cfg.normalCpOther = 144;
    cpPerSlot = [cfg.normalCpFirst; repmat(cfg.normalCpOther, cfg.symbolsPerSlot - 1, 1)];
    cfg.cpLengths = repmat(cpPerSlot, cfg.slotsPerFrame, 1);

    cfg.usefulSymbolDuration = 1 / cfg.subcarrierSpacing;
    cfg.cpDurations = cfg.cpLengths / cfg.sampleRate;
    cfg.slotDuration = 0.5e-3;
    cfg.subframeDuration = 1e-3;
    cfg.frameDuration = 10e-3;
    cfg.samplesPerSlot = 15360;
    cfg.samplesPerSubframe = 30720;
    cfg.samplesPerFrame = 307200;

    cfg.qamOrder = 16;
    cfg.bitsPerQamSymbol = 4;
    cfg.snrDb = 25;

    cfg.wolaLength = 36;

    cfg.aclrMeasurementBandwidth = cfg.occupiedBandwidth;

    halfActive = cfg.numActiveSubcarriers / 2;
    cfg.activeSubcarriers = [-halfActive:-1 1:halfActive].';
    cfg.dcIndex = cfg.fftSize / 2 + 1;
    cfg.activeIndices = cfg.dcIndex + cfg.activeSubcarriers;

    cfg.cfr = cfrConfig(cfrEnabled);
    cfg.cfr.fftSize = cfg.fftSize;
    cfg.cfr.activeIndices = cfg.activeIndices;
    cfg.cfr.cpLengths = cfg.cpLengths;
    cfg.cfr.wolaLength = cfg.wolaLength;

    cfg.pa = paConfig();
    cfg.pa.systemSampleRate = cfg.sampleRate;
    cfg.pa.rateFactor = cfg.pa.modelSampleRate / cfg.sampleRate;
    cfg.pa.resamplerPassbandFrequency =  2 * cfg.channelBandwidth + cfg.aclrMeasurementBandwidth / 2;
    cfg.pa.resamplerStopbandFrequency = cfg.sampleRate / 2;
    cfg.dpd = dpdConfig(cfg.pa, true, dpdLearningArchitecture, cfg.cfr);
    cfg.dpdLearningArchitecture = cfg.dpd.learningArchitecture;
end
