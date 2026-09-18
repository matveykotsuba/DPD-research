function results = paDataPsd()

sampleRate1 = 737.28e6;
sampleRate2 = 737.28e6;

data122 = load('122M_dpd_775.mat');
dataSeed0 = load('seed_0.mat');
dataSeed1 = load('seed_1.mat');

signal122Input = data122.txSignal(:);
signal122Output = data122.rxSignal(:);

seed0Input = dataSeed0.signal_in(:);
seed0Output = dataSeed0.signal_out(:);

seed1Input = dataSeed1.signal_in(:);
seed1Output = dataSeed1.signal_out(:);

[frequency122, psd122Input, info122] = welch(signal122Input, sampleRate1);

[~, psd122Output] = welch(signal122Output, sampleRate1);

[frequencySeed0, psdSeed0Input, infoSeed0] = welch(seed0Input, sampleRate2);

[~, psdSeed0Output] =  welch(seed0Output, sampleRate2);

[frequencySeed1, psdSeed1Input, infoSeed1] = welch(seed1Input, sampleRate2);

[~, psdSeed1Output] = welch(seed1Output, sampleRate2);

reference122 = max([psd122Input; psd122Output]);
referenceSeed0 = max([psdSeed0Input; psdSeed0Output]);
referenceSeed1 = max([psdSeed1Input; psdSeed1Output]);

psd122InputDb = 10 * log10(psd122Input / reference122 + eps);

psd122OutputDb = 10 * log10(psd122Output / reference122 + eps);

psdSeed0InputDb = 10 * log10(psdSeed0Input / referenceSeed0 + eps);

psdSeed0OutputDb = 10 * log10(psdSeed0Output / referenceSeed0 + eps);

psdSeed1InputDb = 10 * log10(psdSeed1Input / referenceSeed1 + eps);

psdSeed1OutputDb = 10 * log10(psdSeed1Output / referenceSeed1 + eps);

figure('Name', 'data PSD', 'Color', 'w');

subplot(3, 1, 1);

plot(frequency122 / 1e6, psd122InputDb, 'b', 'LineWidth', 1.0);

hold on;

plot(frequency122 / 1e6, psd122OutputDb, 'r', 'LineWidth', 1.0);

grid on;
xlabel('Frequency, MHz');
ylabel(' PSD, dB');
title('122M\_dpd\_775.mat');
legend('txSignal', 'rxSignal', 'Location', 'best');
xlim([frequency122(1) frequency122(end)] / 1e6);
ylim([-100 5]);

subplot(3, 1, 2);

plot(frequencySeed0 / 1e6, psdSeed0InputDb, 'b', 'LineWidth', 1.0);

hold on;

plot(frequencySeed0 / 1e6, psdSeed0OutputDb, 'r', 'LineWidth', 1.0);

grid on;
xlabel('Frequency, MHz');
ylabel(' PSD, dB');
title('seed\_0.mat');
legend('signal\_in', 'signal\_out', 'Location', 'best');
xlim([frequencySeed0(1) frequencySeed0(end)] / 1e6);
ylim([-100 5]);

subplot(3, 1, 3);

plot(frequencySeed1 / 1e6, psdSeed1InputDb, 'b', 'LineWidth', 1.0);

hold on;

plot(frequencySeed1 / 1e6, psdSeed1OutputDb, 'r', 'LineWidth', 1.0);

grid on;
xlabel('Frequency, MHz');
ylabel(' PSD, dB');
title('seed\_1.mat');
legend('signal\_in', 'signal\_out', 'Location', 'best');
xlim([frequencySeed1(1) frequencySeed1(end)] / 1e6);
ylim([-100 5]);

results = struct();

results.capture122.frequency = frequency122;
results.capture122.inputPsd = psd122Input;
results.capture122.outputPsd = psd122Output;
results.capture122.info = info122;

results.seed0.frequency = frequencySeed0;
results.seed0.inputPsd = psdSeed0Input;
results.seed0.outputPsd = psdSeed0Output;
results.seed0.info = infoSeed0;

results.seed1.frequency = frequencySeed1;
results.seed1.inputPsd = psdSeed1Input;
results.seed1.outputPsd = psdSeed1Output;
results.seed1.info = infoSeed1;

end
