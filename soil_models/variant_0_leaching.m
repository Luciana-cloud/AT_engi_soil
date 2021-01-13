function vf_ = variant_0_leaching(t,x,p1,q)

%% STOCKS

ATol     = x(1);  % Lighter AT in the system in solution [µg C/L]
AToh     = x(2);  % Heavier AT in the system in solution [µg C/L]
HYo      = x(3);  % HY in the system in solution [µg C/L]

%% PARAMETERS %%

% SORPTION PARAMETERS

K_F_AT  = p1(9);   % Freundlich coeff of AT
n_F_AT  = p1(10);  % Freundlich exponent of AT
K_F_HY  = p1(11);  % Freundlich coeff of HY
n_F_HY  = p1(12);  % Freundlich exponent of HY

%% VALUES OF CONSTANT %%

th_V    = q(10);            % Average volumetric soil water content [cm^3 cm^-3]
rho_B   = q(11);            % Bulk density of soils [g cm^-3]
lv      = q(12);

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(3,1);
   
vf_(1) = (-lv*ATol)*...
         (1+th_V^(-1)*rho_B*(max(ATol,eps^2))^(n_F_AT-1)*K_F_AT*n_F_AT)^(-1);        % Lighter AT outside the cell [µg C/L]

vf_(2) = (-lv*AToh)*...
         (1+th_V^(-1)*rho_B*(max(AToh,eps^2))^(n_F_AT-1)*K_F_AT*n_F_AT)^(-1);        % Heavier AT outside the cell [µg C/L]

vf_(3) = (-lv*HYo)*...
         (1+th_V^(-1)*rho_B*(max(HYo,eps^2))^(n_F_HY-1)*K_F_HY*n_F_HY)^(-1);         % HY outside the cell [µg C/L]

end