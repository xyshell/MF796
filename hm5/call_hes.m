function mse = call_hes(x, r, s, t, q, k, call_bid, call_ask)
%CAL_ Summary of this function goes here
% calibrate κ, θ, σ, ρ and ν0
% given rf, S0, maturity, q, call price.
kappa=x(1);
theta=x(2);
sigma=x(3);
rho=x(4);
v0=x(5);
call = (call_bid + call_ask) / 2;
par = struct(...
    'rf', r, 'q', q, 'x0', log(s), 'v0', v0, 'kappa', kappa, ...
    'theta', theta, 'sigma', sigma, 'rho', rho);
aux.x0 = par.x0;
aux.N = 1024;
aux.K = k;
aux.damp = 1.5;
cfHes = @(u) cflib(u, t, par, 'Heston');
[C, K] = cf2call(cfHes, aux);
% find index of closest strike
mse = sum((C - call).^2);
end

