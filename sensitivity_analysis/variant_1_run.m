%% MANUAL CALIBRATION %%
function [output] = variant_1_run(p1,mode,out)
%% Fixed Information %%

delta13_I = (-29/1000);
AT_tot    =  0*(8*12/216);                  % Input of AT, including light and heavy isotopologues [mg C/L]
frac_coef = -0.0054;                        % Initial fraction Coefficient for the AT applied in the experiment "?" initial
SN        = (delta13_I+1)*0.011180;         % Initial C13/C12. -29 is taken from Beno´s paper (-29 per mil)

%% Define fixed parameters %% 

q(2) = AT_tot;                                  % Concentration of light AT in the inlet [microg C/L]
q(3) = AT_tot*SN;                               % Concentration of heavier AT in the inlet [microg C/L]
q(4) = 1;                                       % Factor for chemostat (1) and retentostat (0).
q(5) = 8*12/216;                                % mg C/mg AT
q(6) = 8*12/197;                                % mg C/mg HY
q(7) = 1e-15;                                   % Volume of a single bacterium [L]
q(8) = 1+frac_coef;                             % Isotopic fractionation factor "alpha" that comes from the enzymatic activity [-]
q(9) = 0.011180;                                % reference isotope ratio of VPDB (Vienna Pee Dee Belemnite) [-]
q(10)= 0.35;                                    % Average volumetric soil water content [cm^3 cm^-3]
q(11)= 1.2;                                     % Bulk density of soils [g cm^-3]
q(12)= (0.5+0.63)/2/300;

%% Calling Model %%

abstol = 1e-9;
reltol = 1e-9;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:7); %,'Jacobian',@MCPA1_jac);%,'NonNegative',1:7);

T_bioC   = 1.5e10;  % Initial Biomass [cells/Kg]
C_AT     = q(5);     % µg C/µg AT
C_HY     = q(6);     % µg C/µg HY
isf      = q(8);     % Isotopic fractionation factor (epsilon)
Ref      = q(9);     % reference isotope ratio of VPDB
C_tot    = p1(14);   % Total initial MCPA conc. micro g-1 soil
C_L0     = 0;
th_V     = q(10);     % Average volumetric soil water content (cm^3 cm^-3)
rho_B    = q(11);     % Bulk density of soils (g cm^-3)
K_F      = p1(9);
n_F      = p1(10);
K_F_H    = p1(11);
n_F_H    = p1(12);
M        = p1(13); %

fun=@(C_L0) MCPA_init(C_L0,C_tot,K_F,n_F,th_V,rho_B);
        
%  switch off solver progress information on solver progress of fsolve
f_opts = optimoptions('fsolve','Display','none');
        
% calculate initial MCPA in solution phase
C_L0 = fsolve(fun,C_L0,f_opts); % MCPA in solution

ATinitC  = C_L0;  % µg C/ L OUTSIDE

time = linspace(0,10950,10950); % simulation period and vector of output times [d]

        c(1)    = M;                  % Active bacteria A biomass [µg C/L]
        c(2)    = 0;                       % Lighter AT inside the cell in the system in solution [µg C/L]
        c(3)    = 0;                    % Heavier AT inside the cell in the system in solution [µg C/L]
        c(4)    = ATinitC;                       % Lighter AT outside the cell in the system in solution [µg C/L]
        c(5)    = ATinitC*SN;                    % Heavier AT outside the cell in the system in solution [µg C/L]
        c(6)    = 0;                       % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = 0;                       % Hydroxyatrazine outisde the cell [µg C/L]

%% First dilution rate Chemostat %%

tic
if mode == 1  % For model without leaching
    [ty,cu] = ode15s(@variant_1_s,time,c,o_opts,p1',q);
else
    [ty,cu] = ode15s(@variant_1_leaching_s,time,c,o_opts,p1',q);
end

%% Desired outputs for Soil with sorption %%

        AT_outlr     = (cu(:,4).*th_V/rho_B+K_F.*cu(:,4).^n_F)/C_AT;     % AT lighter
        AT_outhr     = (cu(:,5).*th_V/rho_B+K_F.*cu(:,5).^n_F)/C_AT;     % AT heavier 
        AT_outTS     = AT_outlr + AT_outhr; % Total AT
        
        HY_outrTS      = (cu(:,7).*th_V/rho_B+K_F_H.*cu(:,7).^n_F_H)/C_HY; % HY
        HY_outrSS     = cu(:,7)/C_HY; % HY
        
        TimeorS       = (ty);               % Time

%% Target variable %%

if out == 1 % For residual AT in soil as target variable
    output = [AT_outTS(end)];
    [AT_outTS(end)]
else
    output = [HY_outrTS(end)]; % For residual HY in soil as target variable
    [HY_outrTS(end)]
end