% Thermodynamic growth constrains model + Leaching flux %

function vf_ = variant_2_leaching_s(t,x,p1,q)

%% STOCKS

CB       = x(1);  % Active bacteria A biomass [µg C/L]
ATil     = x(2);  % Lighter AT inside the cell [µg C/L]
ATih     = x(3);  % Heavier AT inside the cell [µg C/L]
ATol     = x(4);  % Lighter AT in the system in solution [µg C/L]
AToh     = x(5);  % Heavier AT in the system in solution [µg C/L]
HYi      = x(6);  % HY inside the cell [µg C/L]
HYo      = x(7);  % HY in the system in solution [µg C/L]

%% PARAMETERS %%

k_m_a   = p1(1);  % Maximum degradation rate micro [1/d]
K_m_a   = p1(2);  % Half saturation constant for degradation [µg C/L]
k_m_h   = p1(3);  % Maximum degradation rate micro [1/d]
K_m_h   = p1(4);  % Half saturation constant for degradation [µg C/L]
m       = p1(5);  % Maintenance coeficient [1/d]
Y       = p1(6);  % Substrate uptake efficiency of pesticide degraders bacteria A [1]
fcell   = p1(7);  % Fraction of the biomass that is filter out from the chemostat [-]
reAT    = p1(8);  % Rate of transfer of AT from in and out of the cell [L/d]
reHY    = p1(8);  % Rate of transfer of AT from in and out of the cell [L/d]
M       = p1(13); % Min bacterial concentration in soils [µg C/Kg]

% SORPTION PARAMETERS

K_F_AT  = p1(9);   % Freundlich coeff of AT
n_F_AT  = p1(10);  % Freundlich exponent of AT
K_F_HY  = p1(11);  % Freundlich coeff of HY
n_F_HY  = p1(12);  % Freundlich exponent of HY

%% VALUES OF CONSTANT %%

r_D     = q(1);             % Dilution rate coeficient [1/d]
C_Il    = q(2);             % Concentration of light AT in the inlet [µg C/L]
C_Io    = q(3);             % Concentration of heavier AT in the inlet [µg C/L]
beta    = q(4);             % Factor for chemostat [-]
mAT_C   = q(5);             % µg C/µg AT
mHY_C   = q(6);             % µg C/µg HY
V_u     = q(7);             % Volume of a single bacterium [L]
alpha   = q(8);             % Isotopic fractionation factor (?)
th_V    = q(10);            % Average volumetric soil water content [cm^3 cm^-3]
rho_B   = q(11);            % Bulk density of soils [g cm^-3]
lv      = q(12);

%% BIOKINETIC FUNCTIONS %%

% Rate of degradation of lighter AT [1/d]
r_ATl   = k_m_a*ATil*(K_m_a+ATil+ATih)^(-1);

% Rate of degradation of heavier AT [1/d]
r_ATh   = alpha*k_m_a*ATih*(K_m_a+ATil+ATih)^(-1);

% Rate of degradation of HY [1/d]
r_HY    = k_m_h*exp(-K_m_h/HYi);

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(7,1);

vf_(1) = (CB)*(r_HY*Y)-m*Y*(CB-M);               % Bacteria Biomass [[µg C/L]]

vf_(2) = -r_ATl*fcell/V_u + reAT*((ATol-ATil))*(fcell/V_u) ...
         -ATil*((CB)*(r_HY*Y)-m*Y*(CB-M))/CB;            % Lighter AT inside the cell [µg C/L]

vf_(3) = -r_ATh*fcell/V_u + reAT*((AToh-ATih))*(fcell/V_u) ...
         -ATih*((CB)*(r_HY*Y)-m*Y*(CB-M))/CB;            % Heavier AT inside the cell [µg C/L]
     
vf_(4) = (-reAT*CB*(ATol-ATil)-lv*ATol)*...
         (1+th_V^(-1)*rho_B*(max(ATol,eps^2))^(n_F_AT-1)*K_F_AT*n_F_AT)^(-1);        % Lighter AT outside the cell [µg C/L]

vf_(5) = (-reAT*CB*(AToh-ATih)-lv*AToh)*...
         (1+th_V^(-1)*rho_B*(max(AToh,eps^2))^(n_F_AT-1)*K_F_AT*n_F_AT)^(-1);        % Heavier AT outside the cell [µg C/L]

vf_(6) = r_ATl*fcell/V_u + r_ATh*fcell/V_u - r_HY*fcell/V_u...
       + reHY*(HYo-HYi)*fcell/V_u - HYi*((CB)*(r_HY*Y)-m*Y*(CB-M))/CB;  % HY inside the cell [µg C/L]

vf_(7) = (-reHY*CB*(HYo-HYi)-lv*HYo)*...
         (1+th_V^(-1)*rho_B*(max(HYo,eps^2))^(n_F_HY-1)*K_F_HY*n_F_HY)^(-1);                 % HY outside the cell [µg C/L]

end