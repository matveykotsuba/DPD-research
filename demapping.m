function bits = demapping(symbols)
    symbols = symbols(:);
    scaledI = real(symbols) * sqrt(10);
    scaledQ = imag(symbols) * sqrt(10);

    bitGroups = zeros(length(symbols), 4);
    bitGroups(:, 1) = scaledI < 0;
    bitGroups(:, 2) = scaledQ < 0;
    bitGroups(:, 3) = abs(scaledI) >= 2;
    bitGroups(:, 4) = abs(scaledQ) >= 2;

    bits = reshape(bitGroups.', [], 1);
end
