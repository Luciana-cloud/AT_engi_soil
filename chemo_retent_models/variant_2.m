function vf_ = variant_2(t,x,p1,q)

% Model with thermodynamic growth constraint %

%% STOCKS

CB       = x(1);  % Active bacteria A biomass [µg C/L]
ATil     = x(2);  % Lighter AT inside the cell [µg C/L]
ATih     = x(3);  % Heavier AT inside the cell [µg C/L]
ATol     = x(4);  % Lighter AT in the system in solution [µg C/L]
AToh     = x(5);  % Heavier AT in the system in solution [µg C/L]
HYi      = x(6);  % HY inside the cell [µg C/L]
HYo      = x(7);  % HY in the system in solution [µg C/L]

%% PARAMETERS %%

k_m_a   = 10^p1(1);  % Maximum degradation rate micro [1/d]
K_m_a   = 10^p1(2);  % Half saturation constant for degradation [µg C/L]
k_m_h   = 10^p1(3);  % Maximum degradation rate micro [1/d]
K_m_h   = 10^p1(4);  % Half saturation constant for degradation [µg C/L]
m       = 10^p1(5);  % Maintenance coeficient [1/d]
Y       = 10^p1(6);  % Substrate uptake efficiency of pesticide degraders bacteria A [1]
fcell   = 10^p1(7);  % Fraction of the biomass that is filter out from the chemostat [-]
reAT    = 10^p1(8);  % Rate of transfer of AT from in and out of the cell [L/d]
reHY    = 10^p1(8);  % Rate of transfer of AT from in and out of the cell [L/d]

%% VALUES OF CONSTANT %%

r_D     = q(1);             % Dilution rate coeficient [1/d]
C_Il    = q(2);             % Concentration of light AT in the inlet [µg C/L]
C_Io    = q(3);             % Concentration of heavier AT in the inlet [µg C/L]
beta    = q(4);             % Factor for chemostat [-]
mAT_C   = q(5);             % µg C/µg AT
mHY_C   = q(6);             % µg C/µg HY
V_u     = q(7);             % Volume of a single bacterium [L]
alpha   = q(8);             % Isotopic fractionation factor (?)

%% BIOKINETIC FUNCTIONS %%

% Rate of degradation of lighter AT [1/d]
r_ATl   = k_m_a*ATil*(K_m_a+ATil+ATih)^(-1);

% Rate of degradation of heavier AT [1/d]
r_ATh   = alpha*k_m_a*ATih*(K_m_a+ATil+ATih)^(-1);

% Rate of degradation of HY [1/d]
r_HY    = k_m_h*exp(-K_m_h/HYi);

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(7,1);

vf_(1) = CB*(r_HY*Y - m*Y - r_D*beta);               % Bacteria Biomass [[µg C/L]]

vf_(2) = -r_ATl*fcell/V_u + reAT*((ATol-ATil))*(fcell/V_u) ...
         -ATil*(r_HY*Y - m*Y - r_D*beta);            % Lighter AT inside the cell [µg C/L]

vf_(3) = -r_ATh*fcell/V_u + reAT*((AToh-ATih))*(fcell/V_u) ...
         -ATih*(r_HY*Y - m*Y - r_D*beta);            % Heavier AT inside the cell [µg C/L]
     
vf_(4) = r_D*(C_Il-ATol)-reAT*CB*(ATol-ATil);        % Lighter AT outside the cell [µg C/L]

vf_(5) = r_D*(C_Io-AToh)-reAT*CB*(AToh-ATih);        % Heavier AT outside the cell [µg C/L]

vf_(6) = r_ATl*fcell/V_u + r_ATh*fcell/V_u - r_HY*fcell/V_u...
       + reHY*(HYo-HYi)*fcell/V_u - HYi*(r_HY*Y - m*Y - r_D*beta);  % HY inside the cell [µg C/L]

vf_(7) = -reHY*CB*(HYo-HYi)-r_D*HYo;                 % HY outside the cell [µg C/L]

end