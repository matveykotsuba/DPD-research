function figureHandle = paAmPlot(paInfo)
    dpdEnabled = isfield(paInfo, 'dpdEnabled') && paInfo.dpdEnabled;

    if dpdEnabled
        pointCount = min([24000 length(paInfo.amReferenceInput) length(paInfo.amPaOutput) length(paInfo.amDpdOutput) length(paInfo.amDpdPaOutput)]);
        referenceComplex = paInfo.amReferenceInput(1:pointCount);
        paOutputComplex = paInfo.amPaOutput(1:pointCount);
        dpdOutputComplex = paInfo.amDpdOutput(1:pointCount);
        dpdPaOutputComplex = paInfo.amDpdPaOutput(1:pointCount);
        referenceAmplitude = abs(referenceComplex);
        paOutputAmplitude = abs(paOutputComplex);
        dpdOutputAmplitude = abs(dpdOutputComplex);
        dpdPaOutputAmplitude = abs(dpdPaOutputComplex);
        paPhase = angle(paOutputComplex .* conj(referenceComplex)) * 180 / pi;
        dpdPhase = angle(dpdOutputComplex .* conj(referenceComplex)) * 180 / pi;
        dpdPaPhase = angle(dpdPaOutputComplex .* conj(referenceComplex)) * 180 / pi;
        paMaximum = max([referenceAmplitude; paOutputAmplitude; dpdOutputAmplitude; dpdPaOutputAmplitude]);

        figureHandle = figure('Name', 'PA AM/AM and AM/PM', 'Color', 'w');
        subplot(1, 2, 1);
        paAmLine = plot(referenceAmplitude, paOutputAmplitude, 'b.', 'LineStyle', 'none', 'MarkerSize', 4);
        hold on;
        dpdAmLine = plot(referenceAmplitude, dpdOutputAmplitude, 'g.', 'LineStyle', 'none', 'MarkerSize', 4);
        dpdPaAmLine = plot(referenceAmplitude, dpdPaOutputAmplitude, 'm.', 'LineStyle', 'none', 'MarkerSize', 4);
        plot([0 paMaximum], [0 paMaximum], 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        grid on;
        axis equal;
        xlim([0 paMaximum]);
        ylim([0 paMaximum]);
        xlabel(' input amplitude');
        ylabel(' output amplitude');
        title('PA AM/AM');
        legend([paAmLine dpdAmLine dpdPaAmLine], {'PA', 'DPD', 'DPD+PA'}, 'Location', 'best');

        subplot(1, 2, 2);
        paPmLine = plot(referenceAmplitude, paPhase, 'r.', 'LineStyle', 'none', 'MarkerSize', 4);
        hold on;
        dpdPmLine = plot(referenceAmplitude, dpdPhase, 'g.', 'LineStyle', 'none', 'MarkerSize', 4);
        dpdPaPmLine = plot(referenceAmplitude, dpdPaPhase, 'm.', 'LineStyle', 'none', 'MarkerSize', 4);
        plot([0 max(referenceAmplitude)], [0 0], 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        grid on;
        xlim([0 max(referenceAmplitude)]);
        ylim([-180 180]);
        xlabel('input amplitude');
        ylabel('phase shift, degrees');
        title('PA AM/PM');
        legend([paPmLine dpdPmLine dpdPaPmLine], {'PA', 'DPD', 'DPD+PA'}, 'Location', 'best');
    else
        pointCount = min([24000 length(paInfo.amInput) length(paInfo.amOutput)]);
        paInputComplex = paInfo.amInput(1:pointCount);
        paOutputComplex = paInfo.amOutput(1:pointCount);
        paInputAmplitude = abs(paInputComplex);
        paOutputAmplitude = abs(paOutputComplex);
        paPhase = angle(paOutputComplex .* conj(paInputComplex)) * 180 / pi;
        paMaximum = max([paInputAmplitude; paOutputAmplitude]);

        figureHandle = figure('Name', 'PA AM/AM and AM/PM', 'Color', 'w');
        subplot(1, 2, 1);
        plot(paInputAmplitude, paOutputAmplitude, 'b.', 'LineStyle', 'none', 'MarkerSize', 4);
        hold on;
        plot([0 paMaximum], [0 paMaximum], 'k--', 'LineWidth', 1.0);
        grid on;
        axis equal;
        xlim([0 paMaximum]);
        ylim([0 paMaximum]);
        xlabel(' input amplitude');
        ylabel(' output amplitude');
        title('PA AM/AM');

        subplot(1, 2, 2);
        plot(paInputAmplitude, paPhase, 'r.', 'LineStyle', 'none', 'MarkerSize', 4);
        hold on;
        plot([0 max(paInputAmplitude)], [0 0], 'k--', 'LineWidth', 1.0);
        grid on;
        xlim([0 max(paInputAmplitude)]);
        ylim([-180 180]);
        xlabel('input amplitude');
        ylabel('phase shift, degrees');
        title('PA AM/PM');
    end
end
