%% MANUAL CALIBRATION %%
function [SSE] = retentostat_calibration(p1)
%% Fixed Information %%

delta13_I = (-29/1000);
AT_tot    =  30000*(8*12/216);                  % Input of AT, including light and heavy isotopologues [micro g C/L]
frac_coef = -0.0054;                            % Initial fraction Coefficient for the AT applied in the experiment "?" initial
SN        = (delta13_I+1)*0.011180;             % Initial C13/C12. -29 is taken from Beno´s paper (-29 per mil)

%% Define fixed parameters %% 

q(2) = AT_tot;                                  % Concentration of light AT in the inlet [micro g C/L]
q(3) = AT_tot*SN;                               % Concentration of heavier AT in the inlet [micro g C/L]
q(4) = 1;                                       % Factor for chemostat (1) and retentostat (0).
q(5) = 8*12/216;                                % mg C/mg AT
q(6) = 8*12/197;                                % mg C/mg HY
q(7) = 1e-15;                                   % Volume of a single bacterium [L]
q(8) = 1+frac_coef;                             % Isotopic fractionation factor "alpha" that comes from the enzymatic activity [-]
q(9) = 0.011180;                                % reference isotope ratio of VPDB (Vienna Pee Dee Belemnite) [-]

%% Calling Model %%

abstol = 1e-9;
reltol = 1e-7;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:7); 

T_bioC   = 3.66e11;  % Initial Biomass [cells/L]
fcell    = 10^p1(7);
C_AT     = q(5);     % µg C/µg AT
C_HY     = q(6);     % µg C/µg HY
ATinitC  = 87*C_AT;  % µg C/ L OUTSIDE
HYinitC  = 18*C_HY;  % µg C/ L OUTSIDE
isf      = q(8);     % Isotopic fractionation factor (?)
Ref      = q(9);     % reference isotope ratio of VPDB
T_bioR   = 5.4e10;   % Initial Biomass [cells/L]
ATinitR  = 100*C_AT; % microgr/ L OUTSIDE
HYinitR  = 250*C_HY; % microgr/ L OUTSIDE 

%% First dilution rate Retentostat %%

q(1) = 0.02*24;       % Dilution rate coeficient (1/d)        
t5   = [0 1000];
q(4) = 0;                                       % Factor for chemostat (1) and retentostat (0).

        c(1)    = T_bioR*fcell;             % Active bacteria A biomass [µg C/L]
        c(2)    = ATinitR;                  % Lighter AT inside the cell in the system in solution [µg C/L]
        c(3)    = ATinitR*SN;               % Heavier AT inside the cell in the system in solution [µg C/L]
        c(4)    = ATinitR;                  % Lighter AT outside the cell in the system in solution [µg C/L]
        c(5)    = ATinitR*SN;               % Heavier AT outside the cell in the system in solution [µg C/L]
        c(6)    = HYinitR;                  % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = HYinitR;                  % Hydroxyatrazine outisde the cell [µg C/L]

try
     warning off
tic
[ty5,cu5] = ode15s(@variant_1,t5,c,o_opts,p1',q);
% [ty5,cu5] = ode15s(@variant_2,t5,c,o_opts,p1',q); % activated for variant
% 2
 catch ME
     warning off
end
if length(cu5) < length(t5)
    cu5 = ones(length(t5),length(c))*1e+99;
end

if isreal(cu5)==0
    cu5 = ones(length(t5),length(c))*1e+99;    
end

%% Desired outputs for Retentostat %%

        Lcellr       = (cu5(end,1));    % Biomass 
        AT_outlr     = (cu5(end,4))/C_AT;     % AT lighter
        AT_outhr     = (cu5(end,5))/C_AT;     % AT lighter
        AT_Tot2      = AT_outlr + AT_outhr;
        HY_outr      = (cu5(end,7))/C_HY;     % AT lighter
        OUT_R        = vertcat(AT_Tot2,HY_outr,Lcellr);


%% Enrichment Factor Retentostat %% 

d_not_i     = delta13_I;  % epsilon-notation for the inlet.
d_not_oR    = (cu5(end,5)/cu5(end,4)/Ref)-1;  % epsilon-notation for the outlet in the chemostat.
EF_Rt       = (d_not_i - d_not_oR)*1e3;

%% Calling Data %

% Retentostat %

load('retentostat_data.txt')

DatamR   = retentostat_data(:,1);
DatasR   = retentostat_data(:,2);

DAT_R    = vertcat(DatamR(1:2),DatamR(3)*fcell);
DAT_Rsd  = vertcat(DatasR(1:2),DatasR(3)*fcell);
EF_R     = -0.45;
sd_R     = 0.36;

%% SSE %%

SSE  = sum(vertcat((OUT_R-DAT_R)./DAT_Rsd,(EF_R-EF_Rt)/sd_R).^2)

end