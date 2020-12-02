%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Waqtel modules
%
%%%%%%%%%%%%%%%%
% BIOMASS Module
%%%%%%%%%%%%%%%%

% Variables
%
% PHY --> Phytoplankton Biomass (micgClh"a"/L)
% PO4 --> Dissolved mineral phosphorus assimilable by phytoplankton(mgP/L)
% POR --> Degradable phosphorus not assimilable by phytoplankton(mgP/L)
% NO3 --> Dissolved mineral nitrogen assimilable by phytoplankton(mgN/L)
% NOR --> Degradable nitrogen not assimilable by phytoplankton(mgN/L)
%
% Formulation
%
%   dPHY/dt = (CP - DP)*PHY
%
%   dPO4/dt = fp(dtp*DP - CP)*PHY + k1*g2*POR
%
%   dPOR/dt = pf(1 - dtp)*DP*PHY - k1*g2*POR - FPOR/H
%
%   dNO3/dt = fn(dtn*DP - CP)*PHY + k2*g2*NOR
%
%   dNOR/dt = fn(1 - dtn)*DP*PHY - k2*g2*NOR - FNOR/H
%
% CONSTANST DESCRIPTION
%
%%%%%%%
% CP, algal growth rate (d^-1)
%
%       CP=Cmax*RAY*g1*LNUT*α1
%
% Cmax = algal growth maximum rate at 20°C; typical vale of 2
%
% RAY = effect of sunlight on algal growth; ranges between 0 and 1
% 
%   %Smith formula averaged in the vertical
%
%       RAY = (1/ke*H)*log((Io + sqrt(IK^2 + Io^2))/(Ih + sqrt(IK^2 + Ih^2)))
%
%   With:
%
%   ke = 1.7/Zs;    Extinction coefficient of solar rays in water (Aktins formula), Zs = Secchi depth
% or
%   ke = kpe + beta*PHY;    Moss relation, kpe is the vegetal turbidity without phytoplankton (m^-1)
%                            and beta the Moss Coefficient 0.015 aprox.
%   IK,     Calibration parameter (W/m^2), of an order of magnitude 100
%
%   Io,     The flux density of solar radiation on the surface (W/m^2)
%
%   Ih,     The flux density of solar radiation at the bed bottom (W/m^2)
%       Ih = Io*exp(-ke*H)
%
% g1, effect of temperature on algal growth
%
%   g1 = T/20;   where T is temperature, valid for 5 < T < 25
%
% LNUT
%
%   LNUT = min(PO4/(KP + PO4) , NO3/(KN + NO3))
%
% Whit
%
%   KP,     Phoshate half-saturation constant (mg/L), 0.005 aprox., 0.006 Cienaga Santa Marta paper
%
%   KN,     Nitrate half-saturation constant (mg/L), 0.03 aprox., 0.025 Cienaga Santa Marta paper
%
%
% Note:  Nutrients affect phytoplankton growth PHY only by limiting the factorLNU T.  
% When[PO4] and [NO3] are high enough,LNUT is close to 1 and phytoplankton evolution 
% no longerdepends on nutrients.   In this case,  there is no need to model the cycles 
% of phosphorus andnitrogen for simulating the evolution of phytoplankton.
%
% α1, water toxicity coefficient for algae (1.0 means no toxicity)
%
%%%%%%%
%
%%%%%%%
% DP, algae disappearance rate (d^-1)
%
%       DP = (RP + MP)*g2
%
% RP = algal biomass respiration rate at 20°C (d^-1)
%
% MP = algal biomass disappearance rate at 20°C (d^-1)
%
%   MP = M1 + M2*PHY + α2
%
% With
%       M1 and M2, algal mortality coefficients at 20°C
%       α2 = water toxicity coefficient for algae 
%
% g2 = effect of temperature on algal disappearance
%
%   g2 = (1.05)^(T-20);   where T is temperature, valid for 5 < T < 25
%
%%%%%%%
%
%%%%%%%
% Nitric and phosphoric nutrients
%
% fp, average proportion of phosphorus in the cells of living phytoplankton
% (mgP/micgChl"a")
%
% dtp, proportion of directly assimilable phosphorus in dead phytoplankton
% (%)
%
% k1, transformation rate of POR into PO4 through bacterial mineralization
% (d^-1)
%
% k2, transformation rate of NOR into NO3 through heterotrophic and autotrophic bacterial mineralization
% (d^-1)
%
% fn, average proportion of directly assimilable nitrogen in living phtyoplankton
% (mgN/micgChl"a")
%
% dtn, proportion of directly assimilable nitrogen in dead phytoplankton
% (%)
%
% FPOR, deposition flux of non-algal organic phosphorus
% (g/m^2/s)
%
%       FPOR = WPOR*POR
%
%       WPOR, sedimentation velocity of non-algal organic phosphorus (m/s)
%
% FNOR, deposition flux of non-algal organic nitrogen
% (g/m^2/s)
%
%       FNOR = WNOR*NOR
%
%       WNOR, sedimentation velocity of non-algal organic nitrogen (m/s)
%
%%%%%%%%%%%%%%%%%%%%%%%
% TOY MODEL
%%%%%%%%%%%%%%%%%%%%%%%

% Typical values of velocity(m/s), water depth(m) and temperature(°C)
U = 0.8;
H = 4.0;
T = 25;

%ALGAL GROWTH, CP (d^-1)

% Cmax = algal growth maximum rate at 20°C; typical vale of 2
Cmax = 2;

% ke, extinction coefficient of solar rays in water (Aktins formula)
Zs = 0.15;      %Secchi depth
ke = 1.7/Zs;

% IK,     Calibration parameter (W/m^2)
IK = 100;

%   Io,     The flux density of solar radiation on the surface (W/m^2)
Io = 160;

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
KP = 0.005;

%   KN,     Nitrate half-saturation constant (mg/L), 0.03 aprox., 0.025 Cienaga Santa Marta paper
KN = 0.025;

% α1, water toxicity coefficient for algae (1.0 means no toxicity)
alfa1 = 1.0;


%ALGAL DISAPPEARANCE, DP (d^-1)
% On the order of 0.5 (d^-1)

% RP = algal biomass respiration rate at 20°C (d^-1)
RP = 0.01;

% Algal mortality coefficients M1 and M2
% MP = M1 + M2*PHY + α2
M1 = 0.01;
M2 = 0.001;

% α2, water toxicity coefficient for algae (0.0 means no toxicity)
alfa2 = (1.0 - alfa1);

% g2, effect of temperature on algal disappearance
g2 = (1.05)^(T-20);

%FOR NITRIC AND PHOSPHORIC NUTRIENTS
% fp, average proportion of phosphorus in the cells of living phytoplankton
% (mgP/micgChl"a"); 18 (mgC/mgChl"a") according with Cienaga Grande Paper
fp = (18/1000);

% dtp, proportion of directly assimilable phosphorus in dead phytoplankton
% (%); 0.05 according with Cienaga Grande Paper phytoplankton diet
dtp = 0.5;

% k1, transformation rate of POR into PO4 through bacterial mineralization
% (d^-1)
k1 = 0.01;

% k2, transformation rate of NOR into NO3 through heterotrophic and autotrophic bacterial mineralization
% (d^-1)
k2 = 0.01;

% fn, average proportion of directly assimilable nitrogen in living phtyoplankton
% (mgP/micgChl"a"); 18 (mgC/mgChl"a") according with Cienaga Grande Paper
fn = (18/1000);

% dtn, proportion of directly assimilable nitrogen in dead phytoplankton
% (%); 0.95 according with Cienaga Grande Paper phytoplankton diet
dtn = 0.5;

% WPOR, sedimentation velocity of non-algal organic phosphorus (m/s)
WPOR = (0.1/86400);
% WPOR = (1.0/86400);

% WNOR, sedimentation velocity of non-algal organic nitrogen (m/s)
WNOR = (0.1/86400);
% WNOR = (1.0/86400);

%%%%%%%
%
%%%%%%%
%%%% INITIAL VALUES %%%%
PHYo = 240;     % Phytoplankton Biomass (micgClh"a"/L)
PO4o = 0.06;    % Dissolved mineral phosphorus assimilable by phytoplankton(mgP/L)
PORo = 0.03;    % Degradable phosphorus not assimilable by phytoplankton(mgP/L)
NO3o = 0.03;    % Dissolved mineral nitrogen assimilable by phytoplankton(mgN/L)
NORo = 0.01;    % Degradable nitrogen not assimilable by phytoplankton(mgN/L)

%Tiempo
to = 0;               %segundos
tf = 86400*60;         %segundos
dt = 60;              %segundos
time = to:dt:tf;
tsteps = length(time);

%Storage of results
RES = zeros(tsteps,5);
j=1;
RES(j,1) = PHYo;
RES(j,2) = PO4o;
RES(j,3) = PORo;
RES(j,4) = NO3o;
RES(j,5) = NORo;

%Loop temporal
for t = to+dt:dt:tf
    
    %RIGHT HAND SIDES
    
    % Phytoplankton
    
    % Algae growth
    % ke, extinction coefficient of solar rays in water (Cienaga Grande Paper)
    ke = 2.8 + 0.028*PHYo;
    
    %   Ih,     The flux density of solar radiation at the bed bottom (W/m^2)
    Ih = Io*exp(-ke*H);

    %RAY effect
    RAY = (1/(ke*H))*log((Io + sqrt(IK^2 + Io^2))/(Ih + sqrt(IK^2 + Ih^2)));
    
    %LNUT
    LNUT = min(PO4o/(KP + PO4o) , NO3o/(KN + NO3o));
    
    %CP
    CP = Cmax*RAY*g1*LNUT*alfa1;
    
    % Algae Dissapearance
    % MP = algal biomass disappearance rate at 20°C (d^-1)
    MP = M1 + M2*PHYo + alfa2;
    
    %DP
    DP = (RP + MP)*g2;
    CP - DP
    
    % RHS
    PHYrhs = (1/86400)*(CP - DP)*PHYo;    
    
    % Phosphorus
    % PO4 RHS
    PO4rhs = (1/86400)*(fp*(dtp*DP - CP)*PHYo + k1*g2*PORo);
    
    % POR RHS
    FPOR = WPOR*PORo;
    PORrhs = (1/86400)*(fp*(1-dtp)*DP*PHYo - k1*g2*PORo) - FPOR/H;
    
    %Nitrogen
    %NO3 RHS
    NO3rhs = (1/86400)*(fn*(dtn*DP - CP)*PHYo + k2*g2*NORo);
    
    %NOR RHS
    FNOR = WNOR*NORo;
    NORrhs = (1/86400)*(fn*(1-dtn)*DP*PHYo - k2*g2*NORo) - FNOR/H;
    
    %%%% Advancing in time
    PHYo = PHYo + PHYrhs*dt;
    PO4o = PO4o + PO4rhs*dt;
    PORo = PORo + PORrhs*dt;
    NO3o = NO3o + NO3rhs*dt;
    NORo = NORo + NORrhs*dt;
    
    %Storage of results
    j=j+1;
    RES(j,1) = PHYo;
    RES(j,2) = PO4o;
    RES(j,3) = PORo;
    RES(j,4) = NO3o;
    RES(j,5) = NORo;    
    
end

figure(1)

plot(time,RES(:,1))
xlabel('Time (s)')
ylabel('Cocentration (micg Chl"a"/L)')
legend('Phytoplankton','Location','Northeast')

figure(2)

plot(time,RES(:,2:5))
xlabel('Time (s)')
ylabel('Cocentration (mg/L)')
legend('PO4','POR','NO3','NOR','Location','Northeast')


















