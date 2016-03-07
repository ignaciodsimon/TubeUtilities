function H = H_triode_classA(frequencies)%, Vgk)

    % Tube characteristics: Miller capacitances
    Cga = 1.7 * 10^-12;     % Increasing Cgk and Cak one order of magnitude
    Cgk = 1.8 * 10^-12;     % the reduction on hi-freq is around 0.5 dB,
    Cak = 1.9 * 10^-12;     % while that increase on Cga produces a loss in
                            % hi-freq of around 14 dB, so it introduces
                            % around 10^(14.5/20) = 5.3 times more losses.

    % Bias data
%     Ri = (352.5 - 220.0) / (1.625 - 0.15) * 1000; % Ri = ~89.8 k?
%     Gm = (1.25 - 0.15) / (-2 - -3) / 1000; % Gm = ~1.1 mA/V
    Ibias = 0.505 * 10^-3;
    Vbias = -2.5;

    % Use the Vgk value to obtain the corresponding Gm and Ri
    Ri = (352.5 - 220.0) / (1.625 - 0.15) * 1000;
    Gm = (1.25 - 0.15) / (-2 - -3) / 1000;

    % Power supply characteristics
    V_supply = 300;
    Rps = 250;
    Cps = 2 * 47 * 10^-6;

    % Circuit characteristics
    Ra = 300 / 3.5 * 1000;
    Co = 1 * 470 * 10^-9;
    Rl = 500000;
    Cin = 1 * 470 * 10^-9;
    Rin = 5000;     % Added this for the output Z of source
    Rg = 500000;
    Rb = 10000;
    Rk = abs(Vbias) / Ibias;
    Ck = 2.2 * 100 * 10^-6;

    % Computation of capacitors impedances
    ZCin = Rin + 1 ./ (1i * 2 * pi .* frequencies * Cin);
    ZCgk = 1 ./ (1i * 2 * pi .* frequencies * Cgk);
    ZCo = 1 ./ (1i * 2 * pi .* frequencies * Co);
    ZCga = 1 ./ (1i * 2 * pi .* frequencies * Cga);
    ZCak = 1 ./ (1i * 2 * pi .* frequencies * Cak);
    ZCk = 1 ./ (1i * 2 * pi .* frequencies * Ck);
    ZCathode = (ZCk * Rk) ./ (ZCk + Rk);
    ZCps = 1 ./ (1i * 2 * pi .* frequencies * Cps);

    % Computation of intermediate parts
    alpha = 1 ./ (Ra * ( (1 / Rps) + (1 ./ ZCps) + (1 / Ra) ));
    beta = (1 ./ ZCga) + (1 / Ri) + (1 ./ ZCak) + (1 ./ Ra);
    gamma = -Gm + (-1 / Ri) + (-1 ./ ZCak);
    delta = (-1 ./ ZCga) + Gm;
    epsilon = (1 ./ ZCgk) + Gm;
    eta = (-1 / Ri) + (-1 ./ ZCak);
    xsi = (-1 ./ ZCgk) -Gm + (-1 / Ri) + (-1 ./ ZCak) + (-1 ./ ZCathode);
    theta = (-1 ./ ZCin) + (-1 / Rg) + (-1 / Rb);
    kappa =  (-1 ./ ZCga) + (-1 / Rb) + (-1 ./ ZCgk);
    lambda = 1 + (ZCo / Rl);
    sigma = beta - (alpha / Ra) + (eta .* delta ./ epsilon);
    mhu = lambda .* ( (1 ./ ZCga) + (-eta./(theta .* epsilon * Rb * Rb)) + (eta .* kappa ./epsilon));
    phi = (xsi ./ (theta .* epsilon .* Rb * Rb)) + (-xsi .* kappa ./ epsilon) + (1 ./ ZCgk);
    rho = 1 ./ (theta .* ZCin .* Rb);
    tau = ( (-1/Rl) - (lambda .* sigma) ) ./ ( ( -xsi .* delta ./ epsilon) + gamma);

    % Final response computation
    H = rho ./ ((tau .* phi) + mhu);

    % TODO: check why the DC component is NaN
    if isnan(H(1))
        disp('[H_triode_classA()] Warning: the DC coefficient was not a number! ');
        H(1) = 0;
    end
end
