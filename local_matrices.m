function [K, M] = local_matrices(P1, P2, P3, kappa, q)
    x1 = P1(1); x2 = P2(1); x3 = P3(1);
    y1 = P1(2); y2 = P2(2); y3 = P3(2);

    T = ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2;
    Dx = [y2 - y3; y3 - y1; y1 - y2] / (2 * T);    Dy = [x3 - x2; x1 - x3; x2 - x1] / (2 * T);

    K = 1 / (4 * T) * (kappa(1,1) * (Dx * Dx') + kappa(2,2) * (Dy * Dy'));
    M = q * T / 12 * (ones(3,3) + eye(3));
end