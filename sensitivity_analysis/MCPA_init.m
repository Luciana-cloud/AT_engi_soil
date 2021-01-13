%% SOLUTION MCPA - INITIAL CONDITIONS %%

% C_L = MCPA in the solution phase, K_F*C_L^n_F = MCPA in the sorbed phase
% C_T = Total MCPA = C_L*(psy/roh) + K_F*C_L^n_F

function out = MCPA_init(C_L,C_tot,K_F,n_F,th_V,rho_B)
 out = K_F*C_L^n_F + C_L*(th_V/rho_B) - C_tot;
end
