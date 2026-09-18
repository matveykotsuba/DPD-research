function outputSignal = paInterpolate(inputSignal, factor, coefficients)
    inputSignal = inputSignal(:);
    factor = round(factor);

    expandedSignal = complex(zeros(length(inputSignal) * factor, 1));
    expandedSignal(1:factor:end) = inputSignal;

    outputSignal = fftFilterSame(expandedSignal, factor * coefficients);

end
