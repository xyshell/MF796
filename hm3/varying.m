function varying(list, type, par, t, t_list)

rf = 0.025; s0 = 150; sigma = 0.4; v0 = 0.09; 
kappa = 0.5; rho = 0.25; theta = 0.12; k = 150;

clear aux;
aux.x0 = par.x0;

f1 = subplot(1,2,1);
for i = 1:length(list)
    par.(type) = list(i);
    cfHes = @(u) cflib(u, t, par, 'Heston');
    [C K] = cf2call(cfHes,aux);
    vol = blsimpv(s0, K, rf, t, C);
    plot(K, vol)
    hold on;
end 
hold off;
title(f1, 'Implied Vol Skews')
xlabel(f1, 'K')
ylabel(f1, 'Volatility')
legend(f1)

f2 = subplot(1,2,2);
aux.K = k;
for i = 1:length(list)
    par.(type) = list(i);
    C = zeros(length(t_list),1);
    vol = zeros(length(t_list),1);
    for j = 1:length(t_list)
        cfHes = @(u) cflib(u, t_list(j), par, 'Heston');
        C(j) = cf2call(cfHes,aux);
        vol(j) = blsimpv(s0, k, rf, t_list(j), C(j));
    end
    plot(t_list, vol)
    hold on;
end
hold off;
title(f2, 'Vol Term Structure')
xlabel(f2, 'Tau')
ylabel(f2, 'Volatility')
legend(f2)
suptitle(['Varing ',type])

end

