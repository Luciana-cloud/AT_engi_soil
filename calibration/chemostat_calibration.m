%% MANUAL CALIBRATION %%
function [SSE] = chemostat_calibration(p1)
%% Fixed Information %%

delta13_I = (-29/1000);
AT_tot    =  30000*(8*12/216);                  % Input of AT, including light and heavy isotopologues [micro g C/L]
frac_coef = -0.0054;                            % Initial fraction Coefficient for the AT applied in the experiment "?" initial
SN        = (delta13_I+1)*0.011180;             % Initial C13/C12. -29 is taken from Beno´s paper (-29 per mil)

%% Define fixed parameters %% 

q(2) = AT_tot;                                  % Concentration of light AT in the inlet [microg C/L]
q(3) = AT_tot*SN;                               % Concentration of heavier AT in the inlet [microg C/L]
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

        c(1)    = T_bioC*fcell;                  % Active bacteria A biomass [µg C/L]
        c(2)    = ATinitC;                       % Lighter AT inside the cell in the system in solution [µg C/L]
        c(3)    = ATinitC*SN;                    % Heavier AT inside the cell in the system in solution [µg C/L]
        c(4)    = ATinitC;                       % Lighter AT outside the cell in the system in solution [µg C/L]
        c(5)    = ATinitC*SN;                    % Heavier AT outside the cell in the system in solution [µg C/L]
        c(6)    = HYinitC;                       % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = HYinitC;                       % Hydroxyatrazine outisde the cell [µg C/L]

%% First dilution rate Chemostat %%

q(1) = 0.023*24;    % Dilution rate coeficient (1/d)
t    = [0 20];

try
     warning off
tic
    [ty,cu] = ode15s(@variant_1,t,c,o_opts,p1',q); % select for variant 1
%     [ty,cu] = ode15s(@variant_2,t,c,o_opts,p1',q); % select for variant 2
 catch ME
     warning off
end
if length(cu) < length(t)
    cu = ones(length(t),length(c))*1e+99;
end

if isreal(cu)==0
    cu = ones(length(t),length(c))*1e+99;    
end

%% Second dilution rate Chemostat %%

q(1) = 0.032*24;       % Dilution rate coeficient (1/d)        
t1   = [20 31];

        c(1)    = cu(end,1);                        % Active bacteria A biomass [µg C/L]
        c(2)    = cu(end,2);                % Lighter AT inside the cell in the system in solution [µg C/L]
        c(3)    = cu(end,3);                % Heavier AT inside the cell in the system in solution [µg C/L]
        c(4)    = cu(end,4);                % Lighter AT outside the cell in the system in solution [µg C/L]
        c(5)    = cu(end,5);                % Heavier AT outside the cell in the system in solution [µg C/L]
        c(6)    = cu(end,6);                % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = cu(end,7);                % Hydroxyatrazine outisde the cell [µg C/L]

try
     warning off

tic
    [ty1,cu1] = ode15s(@variant_1,t,c,o_opts,p1',q);
%     [ty1,cu1] = ode15s(@variant_2,t,c,o_opts,p1',q);

 catch ME
     warning off
end

if length(cu1) < length(t1)
    cu1 = ones(length(t1),length(c))*1e+99;
end

if isreal(cu1)==0
    cu1 = ones(length(t1),length(c))*1e+99;    
end

%% Third dilution rate %%

q(1) = 0.048*24;       % Dilution rate coeficient (1/d)        
t2   = [31 37];

        c(1)    = cu1(end,1);                % Active bacteria A biomass [µg C/L]
        c(2)    = cu1(end,2);                % Lighter AT inside the cell in the system in solution [µg C/L]
        c(3)    = cu1(end,3);                % Heavier AT inside the cell in the system in solution [µg C/L]
        c(4)    = cu1(end,4);                % Lighter AT outside the cell in the system in solution [µg C/L]
        c(5)    = cu1(end,5);                % Heavier AT outside the cell in the system in solution [µg C/L]
        c(6)    = cu1(end,6);                % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = cu1(end,7);                % Hydroxyatrazine outisde the cell [µg C/L]

try
     warning off
tic
    [ty2,cu2] = ode15s(@variant_1,t,c,o_opts,p1',q);
%     [ty2,cu2] = ode15s(@variant_2,t,c,o_opts,p1',q);
 catch ME
     warning off
end
if length(cu2) < length(t2)
    cu2 = ones(length(t2),length(c))*1e+99;
end

if isreal(cu2)==0
    cu2 = ones(length(t2),length(c))*1e+99;    
end

%% Fourth dilution rate %%

q(1) = 0.056*24;       % Dilution rate coeficient (1/d)        
t3   = [37 42];

        c(1)    = cu2(end,1);                % Active bacteria A biomass [µg C/L]
        c(2)    = cu2(end,2);                % Lighter AT inside the cell in the system in solution [µg C/L]
        c(3)    = cu2(end,3);                % Heavier AT inside the cell in the system in solution [µg C/L]
        c(4)    = cu2(end,4);                % Lighter AT outside the cell in the system in solution [µg C/L]
        c(5)    = cu2(end,5);                % Heavier AT outside the cell in the system in solution [µg C/L]
        c(6)    = cu2(end,6);                % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = cu2(end,7);                % Hydroxyatrazine outisde the cell [µg C/L]
            
try
     warning off
tic
    [ty3,cu3] = ode15s(@variant_1,t,c,o_opts,p1',q);
%     [ty3,cu3] = ode15s(@variant_2,t,c,o_opts,p1',q);
 catch ME
     warning off
end
if length(cu3) < length(t3)
    cu3 = ones(length(t3),length(c))*1e+99;
end

if isreal(cu3)==0
    cu3 = ones(length(t3),length(c))*1e+99;    
end    

%% Fifth dilution rate %%

q(1) = 0.068*24;       % Dilution rate coeficient (1/d)        
t4   = [42 45];

        c(1)    = cu3(end,1);                % Active bacteria A biomass [µg C/L]
        c(2)    = cu3(end,2);                % Heavier AT inside the cell in the system in solution [µg C/L]
        c(3)    = cu3(end,3);                % Lighter AT inside the cell in the system in solution [µg C/L]
        c(4)    = cu3(end,4);                % Heavier AT outside the cell in the system in solution [µg C/L]
        c(5)    = cu3(end,5);                % Lighter AT outside the cell in the system in solution [µg C/L]
        c(6)    = cu3(end,6);                % Hydroxyatrazine INSIDE the cell [µg C/L]
        c(7)    = cu3(end,7);                % Hydroxyatrazine outisde the cell [µg C/L]     

try
     warning off
tic
    [ty4,cu4] = ode15s(@variant_1,t,c,o_opts,p1',q);
%     [ty4,cu4] = ode15s(@variant_2,t,c,o_opts,p1',q);
 catch ME
     warning off
end
if length(cu4) < length(t4)
    cu4 = ones(length(t4),length(c))*1e+99;
end

if isreal(cu4)==0
    cu4 = ones(length(t4),length(c))*1e+99;    
end

%% Desired outputs for Chemostat %%
        
        Lcellc       = vertcat(cu(end,1),cu1(end,1),cu2(end,1),cu3(end,1),cu4(end,1));              % Biomass 
        AT_outlc     = vertcat(cu(end,4),cu1(end,4),cu2(end,4),cu3(end,4),cu4(end,4))/C_AT;         % AT lighter
        AT_outhc     = vertcat(cu(end,5),cu1(end,5),cu2(end,5),cu3(end,5),cu4(end,5))/C_AT;         % AT lighter
        AT_Tot       = AT_outlc + AT_outhc;
        HY_outc      = vertcat(cu(end,7),cu1(end,7),cu2(end,7),cu3(end,7),cu4(end,6))/C_HY;         % AT lighter
        OUT_C        = vertcat(AT_Tot,HY_outc,Lcellc);                                              % Time

%% Enrichment Factor Chemostat %% 

d_not_i = delta13_I;                        % epsilon-notation for the inlet.
d_not_o = ((cu(end,5)/cu(end,4))/Ref)-1;    % epsilon-notation for the outlet in the chemostat.
EF_Ch   = (d_not_i - d_not_o)*1e3;

%% Calling Data %

% Chemostat %

load('chemostat_data.txt')
EF_C  = -5.36;
sd_C  = 0.2;
Datam = chemostat_data(:,1);
Datas = chemostat_data(:,2);

DAT_C   = vertcat(Datam(1:10),Datam(11:15)*fcell);
DAT_Csd = vertcat(Datas(1:10),Datas(11:15)*fcell);

%% SSE %%

SSE     = sum(vertcat((OUT_C-DAT_C)./DAT_Csd,(EF_C-EF_Ch)/sd_C).^2)  

end