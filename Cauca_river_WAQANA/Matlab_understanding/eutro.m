function [RES1,RES2,RES3,RES4,RES5,RES6,RES7,RES8] = eutro(params,init,tf,dt)
%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Waqtel modules
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Eutro Module of Waqtel in OD
%
% Entries
%
% % params = [U,H,T,alfa1,kpe,beta,BEN,M1,M2,k520,k120,KN,KP,n,f,IK,dtn,dtp,fn,fp,k620...
%           k320,RP,WNOR,WLOR,WPOR,Io,Cmax]
% Initial conditions
% init = [PHYo,PO4o,PORo,NO3o,NORo,NH4o,Lo,O2o]
%
% tf, final time in seconds
%
% dt, time step in seconds
%
% Output
%
%         I1,I2,I3,....;P1,P2,P3,......;T1,T2,T3,T4,.....,FT
%         ___________________________________________________
%   PHYo |
%   PO4  |
%   POR  |
%   NO3  |
%   NOR  |
%   NH4  |
%   L    |
%   O2   |
%
%   I1,...; INITIAL CONDITIONS
%   P1,...; PARAMETERS
%   T1,...,TF; RESULTS FOR TIME T
%
% Typical values of velocity(m/s), water depth(m) and temperature(°C)
U = params(1);
H = params(2);
T = params(3);

%ALGAL GROWTH, CP (d^-1)

% Cmax = algal growth maximum rate at 20°C; typical vale of 2
Cmax = params(28);

% ke, extinction coefficient of solar rays in water (Aktins formula)
% Zs = 0.15;      %Secchi depth
% ke = 1.7/Zs;
kpe = params(5);
beta = params(6);

% IK,     Calibration parameter (W/m^2)
IK = params(16);

%   Io,     The flux density of solar radiation on the surface (W/m^2)
Io = params(27);

% maxI = 160;
%En Hernandez-Jimenez se muestra la variación a lo largo de un día de la 
%intensidad lumínica
% Io = [0;0; 0; 0; 0; 0; 0.01*maxI; 0.09*maxI; 0.18*maxI; 0.3*maxI; 0.41*maxI;...
%     0.48*maxI; 0.5*maxI; 0.49*maxI; 0.45*maxI; 0.41*maxI; 0.36*maxI; 0.25*maxI;...
%     0.03*maxI;0;0;0;0;0;0];

%   Ih,     The flux density of solar radiation at the bed bottom (W/m^2)
% Ih = Io*exp(-ke*H);

%RAY effect
% RAY = (1/(ke*H))*log((Io + sqrt(IK^2 + Io^2))/(Ih + sqrt(IK^2 + Ih^2)));

% g1, effect of temperature on algal growth
g1 = T/20;

%   KP,     Phoshate half-saturation constant (mg/L), 0.005 aprox., 0.006 Cienaga Santa Marta paper
KP = params(13);

%   KN,     Nitrate half-saturation constant (mg/L), 0.03 aprox., 0.025 Cienaga Santa Marta paper
KN = params(12);

% α1, water toxicity coefficient for algae (1.0 means no toxicity)
alfa1 = params(4);


%ALGAL DISAPPEARANCE, DP (d^-1)
% On the order of 0.5 (d^-1)

% RP = algal biomass respiration rate at 20°C (d^-1)
RP = params(23);

% Algal mortality coefficients M1 and M2
% MP = M1 + M2*PHY + α2
M1 = params(8);
M2 = params(9);

% α2, water toxicity coefficient for algae (0.0 means no toxicity)
alfa2 = (1.0 - alfa1);

% g2, effect of temperature on algal disappearance
g2 = (1.05)^(T-20);

%FOR NITRIC AND PHOSPHORIC NUTRIENTS
%
% fp, average proportion of phosphorus in the cells of living phytoplankton
% (mgP/micgChl"a"); 18 (mgC/mgChl"a") according with Cienaga Grande Paper
fp = params(20);

% dtp, proportion of directly assimilable phosphorus in dead phytoplankton
% (%); 0.05 according with Cienaga Grande Paper phytoplankton diet
dtp = params(18);

% k320, transformation rate of POR into PO4 through bacterial mineralization
% (d^-1)
k320 = params(22);

% k620, transformation rate of NOR into NO3 through heterotrophic and autotrophic bacterial mineralization
% (d^-1)
k620 = params(21);

% fn, average proportion of nitrogen in living phtyoplankton
% (mgN/micgChl"a"); 18 (mgC/mgChl"a") according with Cienaga Grande Paper
fn = params(19);

% dtn, proportion of directly assimilable nitrogen in dead phytoplankton
% (%); 0.95 according with Cienaga Grande Paper phytoplankton diet
dtn = params(17);

% n, quantity of oxygen consumed by nitrification
% (mgO2/mgNH4), paper Ciénaga Grande
n = params(14);

% k520, kinetics of nitrification at 20°C (d^-1) 
k520 = params(10);

% WPOR, sedimentation velocity of non-algal organic phosphorus (m/s)
WPOR = params(26);
% WPOR = (1.0/86400);

% WNOR, sedimentation velocity of non-algal organic nitrogen (m/s)
WNOR = params(24);
% WNOR = (1.0/86400);

% ORGANIC LOAD
% k120, kinetic degradation constant for the Organic Load at 20°C
% (d^-1)
k120 = params(11);

% g3, effect of temperature on the degradation of Organic Load
g3 = (1.047)^(T-20);

% WLOR, sedimentation velocity of Organic Load (m/s)
WLOR = params(25);

% DISOLVED OXYGEN
% f, Oxygen quantity produced by photosynthesis
% (mgO2/micgChl"a")
f = params(15);

% BEN, Bhentic demand (gO2/m2/d) %
BEN = params(7);
BEN = BEN*(1.065)^(T-20);

% k2 (d^-1) % O'Connor and Dobbins formula
k2 = 3.9*U^0.5*H^(-1.5);

% g4, effect of temperature on natural reaeration
g4 = (1.025)^(T-20);

% Cs, Oxygen concentration at saturation (mgO2/L) % Elmore and Hayes formula
Cs = 14.652 - 0.41022*T + 0.007991*(T^2) - 7.7774e-5*(T^3);

%%%%%%%
%
%%%%%%%
%%%% INITIAL VALUES %%%%
PHYo = init(1);     % Phytoplankton Biomass (micgClh"a"/L)
PO4o = init(2);    % Dissolved mineral phosphorus assimilable by phytoplankton(mgP/L)
PORo = init(3);    % Degradable phosphorus not assimilable by phytoplankton(mgP/L)
NO3o = init(4);    % Dissolved mineral nitrogen assimilable by phytoplankton(mgN/L)
NORo = init(5);    % Degradable nitrogen not assimilable by phytoplankton(mgN/L)
NH4o = init(6);     %Nitrógeno amoniacal (mgH4/L)
Lo = init(7);        %Carga orgáica (mgO2/L)
O2o = init(8);      %Oxígeno disuelto (mgO2/L)

%Tiempo
to = 0;               %segundos
%tf = 86400*60;         %segundos
%dt = 60;              %segundos
time = to:dt:tf;
tsteps = length(time);

%Storage of results
% RES1 = zeros(1,tsteps-1 + length(params) + length(init));
% RES2 = zeros(1,tsteps-1 + length(params) + length(init));
% RES3 = zeros(1,tsteps-1 + length(params) + length(init));
% RES4 = zeros(1,tsteps-1 + length(params) + length(init));
% RES5 = zeros(1,tsteps-1 + length(params) + length(init));
% RES6 = zeros(1,tsteps-1 + length(params) + length(init));
% RES7 = zeros(1,tsteps-1 + length(params) + length(init));
% RES8 = zeros(1,tsteps-1 + length(params) + length(init));

RES1 = zeros(1,2 + length(params) + length(init));
RES2 = zeros(1,2 + length(params) + length(init));
RES3 = zeros(1,2 + length(params) + length(init));
RES4 = zeros(1,2 + length(params) + length(init));
RES5 = zeros(1,2 + length(params) + length(init));
RES6 = zeros(1,2 + length(params) + length(init));
RES7 = zeros(1,2 + length(params) + length(init));
RES8 = zeros(1,2 + length(params) + length(init));

RES1(1:length(init)) = init;
RES2(1:length(init)) = init;
RES3(1:length(init)) = init;
RES4(1:length(init)) = init;
RES5(1:length(init)) = init;
RES6(1:length(init)) = init;
RES7(1:length(init)) = init;
RES8(1:length(init)) = init;

RES1(length(init)+1:length(init)+length(params)) = params;
RES2(length(init)+1:length(init)+length(params)) = params;
RES3(length(init)+1:length(init)+length(params)) = params;
RES4(length(init)+1:length(init)+length(params)) = params;
RES5(length(init)+1:length(init)+length(params)) = params;
RES6(length(init)+1:length(init)+length(params)) = params;
RES7(length(init)+1:length(init)+length(params)) = params;
RES8(length(init)+1:length(init)+length(params)) = params;

j = length(init)+length(params);
k = length(init)+length(params);
jmax = tsteps - 1 + length(init)+length(params);

%Loop temporal
for t = to+dt:dt:tf
    
    %RIGHT HAND SIDES
    
    % Phytoplankton
    
    % Algae growth
    % ke, extinction coefficient of solar rays in water (Cienaga Grande Paper)
    ke = kpe + beta*PHYo;
    
    %   Ih,     The flux density of solar radiation at the bed bottom (W/m^2)
    Ih = Io*exp(-ke*H);

    %RAY effect
    RAY = (1/(ke*H))*log((Io + sqrt(IK^2 + Io^2))/(Ih + sqrt(IK^2 + Ih^2)));
    
    %LNUT
    LNUT = min(PO4o/(KP + PO4o) , (NO3o + NH4o)/(KN + NO3o +  NH4o));
    
    %CP
    CP = Cmax*RAY*g1*LNUT*alfa1;
    % lim = min([RAY,g1,LNUT,alfa1]);
    % CP = Cmax*lim;
    
    % Algae Dissapearance
    % MP = algal biomass disappearance rate at 20°C (d^-1)
    MP = M1 + M2*PHYo + alfa2;
    
    %DP
    DP = (RP + MP)*g2;
    %CP - DP
    
    % RHS
    PHYrhs = (1/86400)*(CP - DP)*PHYo;    
    
    % Phosphorus
    % PO4 RHS
    PO4rhs = (1/86400)*(fp*(dtp*DP - CP)*PHYo + k320*g2*PORo);
    
    % POR RHS
    FPOR = WPOR*PORo;
    PORrhs = (1/86400)*(fp*(1-dtp)*DP*PHYo - k320*g2*PORo) - FPOR/H;
    
    %Nitrogen
    %NO3 RHS
    Rn = (NH4o/(NH4o + NO3o));
    NO3rhs = (1/86400)*(-fn*(1-Rn)*CP*PHYo + k520*g2*NH4o);
    
    %NOR RHS
    FNOR = WNOR*NORo;
    NORrhs = (1/86400)*(fn*(1-dtn)*DP*PHYo - k620*g2*NORo) - FNOR/H;
    
    %NH4 RHS
    NH4rhs = (1/86400)*(fn*(dtn*DP - Rn*CP)*PHYo + k620*g2*NORo - k520*g2*NH4o);
    
    %Organic Load
    FLOR = WLOR*Lo;
    Lrhs = (1/86400)*(f*MP*PHYo - k120*g3*Lo) - FLOR/H;
    
    %Disolved Oxygen
    O2rhs = (1/86400)*(f*(CP - RP*g1)*PHYo - n*k520*g2*NH4o - k120*g3*Lo + k2*g4*(Cs - O2o) - BEN/H);
    
    %%%% Advancing in time
    PHYo = PHYo + PHYrhs*dt;
    PO4o = PO4o + PO4rhs*dt;
    PORo = PORo + PORrhs*dt;
    NO3o = NO3o + NO3rhs*dt;
    NORo = NORo + NORrhs*dt;
    NH4o = NH4o + NH4rhs*dt;
    Lo = Lo + Lrhs*dt;
    O2o = O2o + O2rhs*dt;
    
    %Storage of results
    j=j+1;
    if j>jmax-2
        RES1(k+1) = PHYo;
        RES2(k+1) = PO4o;
        RES3(k+1) = PORo;
        RES4(k+1) = NO3o;
        RES5(k+1) = NORo;
        RES6(k+1) = NH4o;
        RES7(k+1) = Lo;
        RES8(k+1) = O2o;
        k=k+1;
    end

end


end