function C = spm_drift(N, K)



n = 0:(N - 1);
C = zeros(N,K);
C(:, 1) = 1/sqrt(N);
for k = 2:K
    C(:, k) = sqrt(2/N) * 10 * cos(pi * (2 * n + 1) * (k - 1)/(2 * N));
end


end

