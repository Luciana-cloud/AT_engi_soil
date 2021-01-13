function vf_ = variant_complex(t,x,p1,q)

%% STOCKS

CB         = x(1);  % Active growing bacteria A biomass [µg C/L]
CBn        = x(2);  % Active non-growing bacteria A biomass [µg C/L]
ATilg      = x(3);  % Lighter AT inside the growing cell [µg C/L]
ATihg      = x(4);  % Heavier AT inside the growing cell [µg C/L]
ATilng     = x(5);  % Lighter AT inside the non-growing cell [µg C/L]
ATihng     = x(6);  % Heavier AT inside the non-growing cell [µg C/L]
ATol       = x(7);  % Lighter AT in the system in solution [µg C/L]
AToh       = x(8);  % Heavier AT in the system in solution [µg C/L]
HYig       = x(9);  % HY inside the growing cell [µg C/L]
HYing      = x(10); % HY inside the non-growing cell [µg C/L]
HYo        = x(11); % HY in the system in solution [µg C/L]

%% PARAMETERS %%

k_m_a   = 10^p1(1);  % Maximum degradation rate micro [1/d]
K_m_a   = 10^p1(2);  % Half saturation constant for degradation [µg C/L]
k_m_h   = 10^p1(3);  % Maximum degradation rate micro [1/d]
K_m_h   = 10^p1(4);  % Half saturation constant for degradation [µg C/L]
m       = 10^p1(5);  % Maintenance coeficient [1/d]
Y       = 10^p1(6);  % Substrate uptake efficiency of pesticide degraders bacteria A [1]
fcell   = 10^p1(7);  % Conversion factor from cells to C [-]
reAT    = 10^p1(8);  % Rate of transfer of AT from in and out of the cell [L/d]
reHY    = 10^p1(8);  % Rate of transfer of AT from in and out of the cell [L/d]
D_cb    = 10^p1(9);  % Rate of transfer of AT from in and out of the cell [L/d]
C_T     = 10^p1(10); % Threshold concentration [µg C/L]
k_r     = 10^p1(11); % Coefficient rate of reactivation [d-1]
k_d     = 10^p1(12); % Coefficient rate of deactivation [d-1]

%% VALUES OF CONSTANT %%

r_D     = q(1);             % Dilution rate coeficient [1/d]
C_Il    = q(2);             % Concentration of light AT in the inlet [µg C/L]
C_Io    = q(3);             % Concentration of heavier AT in the inlet [µg C/L]
beta    = q(4);             % Factor for chemostat [-]
mAT_C   = q(5);             % µg C/µg AT
mHY_C   = q(6);             % µg C/µg HY
V_u     = q(7);             % Volume of a single bacterium [L]
alpha   = q(8);             % Isotopic fractionation factor (?)
n_S     = q(10);            % Switch parameter

%% BIOKINETIC FUNCTIONS %%

% Rate of degradation of lighter AT growing cells [1/d]
r_ATlg    = k_m_a*ATilg*(K_m_a+ATilg+ATihg)^(-1);

% Rate of degradation of lighter AT non-growing cells [1/d]
r_ATlng   = k_m_a*ATilng*(K_m_a+ATilng+ATihng)^(-1);

% Rate of degradation of heavier AT growing cells [1/d]
r_AThg    = alpha*k_m_a*ATihg*(K_m_a+ATilg+ATihg)^(-1);

% Rate of degradation of heavier AT non-growing cells [1/d]
r_AThng   = alpha*k_m_a*ATihng*(K_m_a+ATilng+ATihng)^(-1);

% Density depended growth rate:
k_m_h_T = k_m_h*(CB+CBn)/(D_cb*fcell+CB+CBn);

% Rate of degradation of HY growing cells [1/d]
r_HYg   = k_m_h_T*HYig*(K_m_h+HYig)^(-1);

% Rate of degradation of HY non-growing cells [1/d]
r_HYng  = k_m_h_T*HYing*(K_m_h+HYing)^(-1);

% Switch function for activation/inactivation of growing bacteria pools: []
taug    = (1+exp(C_T^(-1)*n_S^(-1)*(C_T-(ATilg+ATihg))))^(-1);

% Switch function for activation/inactivation of non-growing bacteria pools: []
taung   = (1+exp(C_T^(-1)*n_S^(-1)*(C_T-(ATilng+ATihng))))^(-1);

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(11,1);

% Bacteria Biomass [[µg C/L]]
vf_(1) = CB*(r_HYg*Y - (1 - HYig*(K_m_h+HYig)^(-1))*m*Y - r_D*beta...
       -(1-taug)*k_d)+taung*k_r*CBn;                   

% Non growing Bacteria Biomass [[µg C/L]]
vf_(2) = CBn*(- (1 - HYing*(K_m_h+HYing)^(-1))*m*Y- r_D*beta-taung*k_r)...
       +(1-taug)*k_d*CB;     
   
% Lighter AT inside the growing cell [µg C/L]
vf_(3) = -r_ATlg*fcell/V_u + reAT*(ATol-ATilg)*fcell/V_u ...
         -ATilg*(CB*(r_HYg*Y - (1 - HYig*(K_m_h+HYig)^(-1))*m*Y - r_D*beta...
         -(1-taug)*k_d)+taung*k_r*CBn)/CB-(1-taug)*k_d*ATilg+taung*k_r*ATilng*CBn/CB;             

% Heavier AT inside the growing cell [µg C/L]     
vf_(4) = -r_AThg*fcell/V_u + reAT*(AToh-ATihg)*fcell/V_u ...
         -ATihg*(CB*(r_HYg*Y - (1 - HYig*(K_m_h+HYig)^(-1))*m*Y - r_D*beta...
         -(1-taug)*k_d)+taung*k_r*CBn)/CB-(1-taug)*k_d*ATihg+taung*k_r*ATihng*CBn/CB;     
   
% Lighter AT inside the non-growing cell [µg C/L]
vf_(5) = -r_ATlng*fcell/V_u + reAT*(ATol-ATilng)*fcell/V_u ...
         -ATilng*(CBn*(- (1 - HYing*(K_m_h+HYing)^(-1))*m*Y- r_D*beta-taung*k_r)...
         +(1-taug)*k_d*CB)/CBn+(1-taug)*k_d*ATilg*CB/CBn-taung*k_r*ATilng;      
     
% Heavier AT inside the non-growing cell [µg C/L]
vf_(6) = -r_AThng*fcell/V_u + reAT*(AToh-ATihng)*fcell/V_u ...
         -ATihng*(CBn*(- (1 - HYing*(K_m_h+HYing)^(-1))*m*Y- r_D*beta-taung*k_r)...
         +(1-taug)*k_d*CB)/CBn+(1-taug)*k_d*ATihg*CB/CBn-taung*k_r*ATihng;       
     
% Lighter AT outside the cell [µg C/L] 
vf_(7) = r_D*(C_Il-ATol)-reAT*CB*(ATol-ATilg)-reAT*CBn*(ATol-ATilng);     
     
% Heavier AT outside the cell [µg C/L]
vf_(8) = r_D*(C_Io-AToh)-reAT*CB*(AToh-ATihg)-reAT*CBn*(AToh-ATihng);

% HY inside growing cell [µg C/L]
vf_(9) = r_ATlg*fcell/V_u + r_AThg*fcell/V_u - r_HYg*fcell/V_u...
       + reHY*(HYo-HYig)*fcell/V_u -(m*HYig*(K_m_h+HYig)^(-1))*fcell/V_u...
       - HYig*(CB*(r_HYg*Y - (1 - HYig*(K_m_h+HYig)^(-1))*m*Y - r_D*beta...
       -(1-taug)*k_d)+taung*k_r*CBn)/CB-(1-taug)*k_d*HYig+taung*k_r*HYing*CBn/CB;

% HY inside non-growing cell [µg C/L]
vf_(10)= r_ATlng*fcell/V_u + r_AThng*fcell/V_u - r_HYng*fcell/V_u...
       + reHY*(HYo-HYing)*fcell/V_u -(m*HYing*(K_m_h+HYing)^(-1))*fcell/V_u...
       - HYing*(CBn*(- (1 - HYing*(K_m_h+HYing)^(-1))*m*Y- r_D*beta-taung*k_r)...
       +(1-taug)*k_d*CB)/CBn+(1-taug)*k_d*HYig*CB/CBn-taung*k_r*HYing;

% HY outside the cell [µg C/L]
vf_(11)= -reHY*CB*(HYo-HYig)-reHY*CBn*(HYo-HYing)-r_D*HYo;

end