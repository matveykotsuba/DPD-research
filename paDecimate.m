function outputSignal = paDecimate(inputSignal, factor, coefficients)
    inputSignal = inputSignal(:);
    factor = round(factor);

    filteredSignal = fftFilterSame(inputSignal , coefficients);

    outputSignal = filteredSignal(1:factor:end);
end
