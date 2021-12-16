function V = dn_Carandini1997(param, rcprm, t, stim)

% EXMAPLES ----------------------------------------------------------------
% dt = 0.001;
% t  = dt : dt : 1;
% 
% stim = zeros(1, length(t));
% stim(t > 0.2 & t <= 0.7) = 1;
% 
% param(1) = 0.05; % tau
% param(2) = 0.5;  % w
% rcprm.C  = 2; % capacitance
% rcprm.g0 = 0.1;
% rcprm.k  = 0.8; % determines the effectiveness of the normalization pool

%% COMPUTE IRF

irf1 = t.*exp(-t./param(1)); 
irf1 = irf1./sum(irf1);

irf2 = t.*exp(-t./(param(1)*1.5)); 
irf2 = irf2./sum(irf2);

irf = irf1 - param(2)*irf2;

%% COMPUTE MODEL PREDICTIONS

compute_g = @(g0, k, V) g0./sqrt(1 - k .* V);

delta_V = 0;
V    = zeros(1, length(t));
g    = zeros(1, length(t));

I = convCut(stim, irf, length(t));

% compute the change in potential
for it = 1 : length(t) - 1
    g(it) = compute_g(rcprm.g0, rcprm.k, V(it));
    delta_V   = (I(it) - g(it) * V(it))/rcprm.C;
    V(it + 1) = V(it) + delta_V;
    V = max(V, 0);
end

%% VISUALIZE RESPONSE

% figure (100), clf
% 
% plot(t, stim), hold on
% plot(t, V)
% plot(t, irf)
end