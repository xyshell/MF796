function approxPDFs = blrnd(prc, k, r, T)
% Compute Risk neutral density(pdf) by Breeden-Litzenberger 
% prc: option price (1*n vec)
% k: strike (1*n vec)

% Estimate the implied densities by approximating derivatives.
% Approximate the first derivatives.
dK = diff(k);
prc_dash = diff(prc) ./ repmat(dK, 1, size(prc, 2));
% Approximate the second derivatives.
d2K = dK(2:end);
prc_ddash = diff(prc_dash) ./ repmat(d2K, 1, size(prc_dash, 2));

% Use the approximate derivatives to estimate each density function.
approxPDFs = NaN(size(prc_ddash));
% Set any negative values to zero.
prc_ddash(prc_ddash < 0) = 0;
for i = 1:size(prc_ddash, 2)
    approxPDFs(:, i) = exp(r * T(i)) * prc_ddash(:, i);    
end

end

