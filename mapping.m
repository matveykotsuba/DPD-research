function symbols = mapping(bits)
    bits = bits(:);

    bitGroups = reshape(bits, 4, []).';
    b0 = bitGroups(:, 1);
    b1 = bitGroups(:, 2);
    b2 = bitGroups(:, 3);
    b3 = bitGroups(:, 4);

    inPhase = (1 - 2 * b0) .* (2 - (1 - 2 * b2));
    quadrature = (1 - 2 * b1) .* (2 - (1 - 2 * b3));

    symbols = (inPhase + 1i * quadrature) / sqrt(10);
    symbols = symbols(:);
end
