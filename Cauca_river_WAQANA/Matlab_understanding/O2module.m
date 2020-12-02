%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Waqtel modules
%
%%%%%%%%%%%%
% O2 MODULE
%%%%%%%%%%%%

% Variables
%
% O2 --> Disolved oxygen (mgO2/L)
% L  --> Organic Load (mgO2/L)
% NH4--> Amonia Load (mgH4/L)
%
% Formulation
%
% dL/dt = -k1*L
%
% dNH4/dt = -k4*NH4
%
% dO2/dt = k2(Cs - O2) - k1*L -k4*NH4 + P - R - BEN/h
%
% CONSTANTS DESCRIPTION
%
% k2 = Reaeration process constant (d^-1)
%         type of watercourse                     Interval ofk2(d−1) at 20◦C
%         Small ponds and backwaters              0.10-0.23
%         Sluggish streams and large lakes        0.23-0.35
%         Large streams of low flow velocity      0.35-0.46
%         Large streams of normal flow velocity   0.46-0.69
%         Swift streams                           0.69-1.15
%         Rapids and waterfalls                   > 1.15
%
%         k2 = 5.23*U*h^−1.67                         (Tennessee Valley Authority)
%         k2 = 5.33*U^0.67*h^−1.85                    (Owens et al.)
%         k2 = 0.746*U^2.695*h^−3.085*J^−0.823        (Churchill et al.)
%         k2 = 3.9*U^0.5*h^−1.5                       (O′Connor and Dobbins)
%
%       U = velocity magnitude (m/s); h = water depth (m)
%       J = Energy slope (m)
%
%       *O'Connor and Dobbins formula has the best results for shallow
%       rivers
%
%   %Correction by temperature
%
%       k2= [k2]*(1.0241)^(T−20)  ; with T in Celsius
%
% Cs = Oxygen concentration at saturation in water (mgO2/L)
%       Typical value of 9.0 (mgO2/L) at 20°C
%
%       %Elmore and Hayes formula
%       Cs= 14.652 − 0.41022*T + 0.007991*T^2 − 7.7774*10^(−5)*T^3
%
%       %Montgomery formula
%       Cs = 468 / (31.6 * T)
%
%       With T in Celsius
%
% k1 = Kinetic degradation constant of the organic load (d^-1)
% 
% k4 = Kinetic constant of nitrification (d^-1)
%
% P = Photosyntesis production of oxygen (mgO2/d/l)
%
% R = Plant respiration (mgO2/d/l)
%
% BEN = Bhentic demand (gO2/m2/d)
%
%   %Correction by temperature
%
%       BEN = [BEN]*(1.065)^(T−20)
%
% Typical values at 20°C
%
%         Bottom type                             Typical value ofBEN(gO2/m2/d) at 20◦C
%         Filamentous bacteria (10 g/m2)          7
%         Mud from waste water, near to release   4
%         Mud from waste water, far from release  1.5
%         Estuarine silt                          1.5
%         Sand                                    0.5
%         Mineral soil                            0.007
%

%%%%%%%%%%%%%%%%%%%%%%%
% TOY MODEL
%%%%%%%%%%%%%%%%%%%%%%%

% Typical values of velocity(m/s), water depth(m) and temperature(°C)
U = 0.8;
H = 4.0;
T = 25;

% k2 (d^-1) % O'Connor and Dobbins formula
k2 = 3.9*U^0.5*H^(-1.5);
k2= k2*(1.0241)^(T-20);

% Cs, Oxygen concentration at saturation (mgO2/L) % Elmore and Hayes formula
Cs = 14.652 - 0.41022*T + 0.007991*(T^2) - 7.7774e-5*(T^3);

% k1 (d^-1) % Kinetic degradation constant of the organic load L
% Proemdio estación puente autopista / Rio Negro / Paper anexo
% k1 = (0.078 + 0.085 + 0.154 + 0.039)/4;
k1 = 0.01;

% k4 (d^-1) % Kinetic constant of nitrification
k4 = 0.176;       %Paper cienaga grande

% P Photosyntesis production of oxygen (mgO2/d/l) %
P = 0.02;

% R Plant respiration (mgO2/d/l) %
R = 0.01;

% BEN, Bhentic demand (gO2/m2/d) %
BEN = 1.5;
BEN = BEN*(1.065)^(T-20);


%%%% INITIAL VALUES %%%%
O2o = 6.0;      %Oxígeno disuelto (mgO2/L)
Lo = 2.0;        %Carga orgáica (mgO2/L)
NH4o = 0.06;     %Nitrógeno amoniacal (mgH4/L)

%Tiempo
to = 0;               %segundos
tf = 86400*7;         %segundos
dt = 60;              %segundos
time = to:dt:tf;
tsteps = length(time);

%Storage of results
RES = zeros(tsteps,3);
j=1;
RES(j,1) = O2o;
RES(j,2) = Lo;
RES(j,3) = NH4o;

%Loop temporal
for t = to+dt:dt:tf
    
    %RIGHT HAND SIDES
    
    % Organic Load
    Lrhs = (1/86400)*k1*Lo;

    % Amonia Load
    NH4rhs = (1/86400)*k4*NH4o;
    
    % Disolved oxygen
    O2rhs = (1/86400)*(k2*(Cs - O2o) - k1*Lo - k4*NH4o + P - R - (BEN/H));

    %%%% Advancing in time
    Lo = Lo + Lrhs*dt;
    NH4o = NH4o + NH4rhs*dt;
    O2o = O2o + O2rhs*dt;
    
    %Storage of results
    j=j+1;
    RES(j,1) = O2o;
    RES(j,2) = Lo;
    RES(j,3) = NH4o;
    
end
    
plot(time,RES(:,1:3))
xlabel('Time (s)')
ylabel('Cocentration (mg/L)')
legend('Disolved Oxygen','Organic Load','Amonia Load','Location','Northeast')












