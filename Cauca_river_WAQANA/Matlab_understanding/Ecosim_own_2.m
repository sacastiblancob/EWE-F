%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Ecosim
%
% SECOOOOOOND APPROACH, WITH VULNERABILITIES
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECOSIM DYNAMIC SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: This code only computes trophic values for consumers
%
% Equations
%
% For producers
%
%   dBi/dt = f(Bi) - sum_j(Cij) + Ii - (Mi + Fi + ei)Bi
%
% Where
%
%   Ii = Net Inmigration
%
%   Mi = Natural mortatility
%
%   Fi = Fishing Effort over i group
%
%   ei = Emigration rate
%
%   f(Bi) = ri*Bi/(1 + Bi*hi)
%
%   ri: maximum P/B that i can exhibit when Bi is low, and ri/hi is the
%   maximum net primary production rate for pool i when biomass i is not
%   limiting to production (Bi high). ri is the es_rel_PoB_max given
%
%   hi = ((ri/(P/B)i)-1)/Bi ; where P/B comes from Ecopath calculations
%
% For consumers
%
% dBi/dt = gi*sum_j(Cji) - sum_j(Cij) + Ii - (Mi + Fi + ei)Bi
%
%   gi: Net growth efficiency
%   Ii: Biomass immigration rate
%   Mi: non-predation mortality/metabolic rate
%   Fi: Fishing mortality rate
%   ei: Emmigration rate
%
%   Cij: consumption rate of pool i biomass by pool j organisms,i.e., the
%   'flow' from pool i to pool j per unit time. For primary producers first
%   term is replaced by a biomass-dependent production rate.
%
%  Concept of Foragin Arenas inclusion
%
%   Cij = vij*aij*Bi*Bj/(vij + v'ij + aij*Bj)
%
%   aij: rate of effective search for pool type i by predator j
%   vij,v'ij, with default setting Vij=V'ij are prey behavioral exchange
%       rate parameters.
%
%   Walters, Christensen, Pauly, 1997
%   Structuring Dynamic Models of Exploited Ecosystem from Trophic
%   Mass-balance Assessments
%
%   vij = xij*Qij/Bi;   where xij is the ratio of the maximum instantaneous
%                       mortality rate
%   Then,
%
%   aij = 2*Qij*vij/(vij*Bi*Bj-Qij*Bj)
%
%   Setting low ratios of values for the vulnerability ratios vij leads to
%   "bottom-up" control of flow rates from prey to predators, such that
%   increases in prey productivity will lead to prey biomass increases and 
%   then to increased availability to predators. Setting high values for 
%   vij leads to ‘top- down’ control and ‘trophic cascade’ effects 
%   (Carpenter and Kitchell, 1993); increases in top predator abundance 
%   can lead to depressed abundance of smaller forage fishes, and this in 
%   turn to increases in abundance of invertebrates upon which these forage
%   fishes depend.
%
%   OWN BUILD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading Configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nvars
% nvars = csvread('Ecosim_data/Lab_Scenario.csv',0,1,[0,1,0,1]);
% nvars = csvread('Ecosim_data/Tampa_Bay_Scenario.csv',0,1,[0,1,0,1]);
nvars = csvread('Ecosim_data/Cauca_Scenario.csv',0,1,[0,1,0,1]);

% Reading Ecopath Data
% ep_data = h5read('Ecosim_data/Lab.h5','/ep_data');
% ep_data = h5read('Ecosim_data/Tampa_Bay.h5','/ep_data');
ep_data = h5read('Ecosim_data/Cauca.h5','/ep_data');

v = fieldnames(ep_data);
for i = 1:length(v)
    name = v{i};
    myVar = ep_data.(v{i});
    assignin('base',strcat('ep_',name), myVar)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diet, vulnerability and consumption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diet
%ep_diet = h5read('Ecosim_data/Lab.h5','/ep_diet');
% ep_diet = csvread('Ecosim_data/Lab_DC.csv',3,1);
%ep_diet = csvread('Ecosim_data/Tampa_Bay_DC.csv',3,1);
ep_diet = csvread('Ecosim_data/Cauca_DC.csv',3,1);

%vulnerability
% numpred = csvread('Ecosim_data/Lab_vul.csv',1,1,[1,1,1,1]);
% vul = csvread('Ecosim_data/Lab_vul.csv',3,1);
%numpred = csvread('Ecosim_data/Tampa_Bay_vul.csv',1,1,[1,1,1,1]);
%vul = csvread('Ecosim_data/Tampa_Bay_vul.csv',3,1);
numpred = csvread('Ecosim_data/Cauca_vul.csv',1,1,[1,1,1,1]);
vul = csvread('Ecosim_data/Cauca_vul.csv',3,1);

vul(vul==-999) = 0.0;

%consumption
Q = zeros(nvars,numpred);
for i=1:nvars
    Q(i,:) = ep_consumption(1:numpred)'.*ep_diet(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scenario
% es_conf = csvread('Ecosim_data/Lab_Scenario.csv',4,1);
%es_conf = csvread('Ecosim_data/Tampa_Bay_Scenario.csv',4,1);
es_conf = csvread('Ecosim_data/Cauca_Scenario.csv',4,1);

confnames = ["MaxrelPB","Feed_time","Maxrel_feeding_time","Feed_time_adjust_rate",...
    "Fraction_of_other_mortality","Predator_effect_on_feeding_time","Density-dep_catchability_Qmax/Qo",...
    "QBmax/QBo","Switching_power_parameter","Advected"];

es_rel_PoB_max = es_conf(:,1);
es_Ftime = es_conf(:,2);
es_Ftime_max = es_conf(:,3);
es_Ftime_adjust = es_conf(:,4);
es_M0_pred = es_conf(:,5);
es_risk_time = es_conf(:,6);
es_risk_time(es_risk_time==-999) = 0;
es_Q_maxoQ_0 = es_conf(:,7);
es_QB_maxoQB_0 = es_conf(:,8);
es_switch_power = es_conf(:,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting time and initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial values
B = ep_biomass;

%Tiempo
% to = 0;               %segundos
% tf = 86400*90;         %segundos
% dt = 60;              %segundos
to = 0;               %años
tf = 2;         %años
dt = 0.001;              %años
time = to:dt:tf;
tsteps = length(time);

%Para guardar resultados
RES = zeros(length(ep_biomass),tsteps);

%Net growth efficiency gi
%g = 0.5*ones(nvars,1);
% g = P/Q

%initial storing
k=1;
RES(:,k) = ep_biomass;
k=k+1;

h = zeros(nvars,1);
f = zeros(nvars,1);
% Computing f(Bi) for producers
for i=1:nvars
    if (ep_org_type(i) == 1)
        % h_i = ((r_i*(P/B)_i) - 1)/Bi
        h(i) = (es_rel_PoB_max(i)*ep_PoB(i) - 1)/ep_biomass(i);
        % f_i = r_i*B_i/(1 + B_i*h_i)
        f(i) = (es_rel_PoB_max(i)*ep_biomass(i))/(1 + ep_biomass(i)*h(i));
    end
end

a = zeros(nvars,numpred);
v = zeros(nvars,numpred);
% Computing a_ij
for i=1:nvars
    for j=1:numpred
        if (ep_diet(i,j)~= 0)
            v(i,j) = (ep_QoB(j)*ep_biomass(j)*ep_diet(i,j)*vul(i,j))/ep_biomass(i);
            a(i,j) = (2*ep_QoB(j)*ep_biomass(j)*ep_diet(i,j)*v(i,j))/...
              (v(i,j)*ep_biomass(i)*ep_biomass(j) - ...
              ep_QoB(j)*ep_biomass(j)*ep_diet(i,j)*ep_biomass(j));
        else
            a(i,j) = 0;
        end
        
    end
end

C = zeros(nvars,numpred);
for t = to+dt:dt:tf
   % Computing consumption as C_ij = a_ij*B_i*B_j
    for i=1:nvars
        for j=1:numpred
            if(ep_diet(i,j) ~= 0)
                C(i,j) = a(i,j)*v(i,j)*ep_biomass(i)*ep_biomass(j)/...
                 (2*v(i,j) + a(i,j)*ep_biomass(j));
            else
                C(i,j) = 0;
            end
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTING RIGHT HAND SIDES (ONLY TWO FIRST TERMS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RHS = zeros(nvars,1);
    for i=1:nvars
        if (ep_org_type(i) == 2)
            %g_i = P/Q_i
            RHS(i) = ep_PoQ(i)*sum(C(:,i)) - sum(C(i,:));
        elseif (ep_org_type(i) == 1)
            RHS(i) = f(i) - sum(C(i,:));
        else
            RHS(i) = 0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ep_biomass = ep_biomass + RHS*dt;
    %ep_biomass(ep_biomass < 0) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    RES(:,k) = ep_biomass;
    k = k+1;
    
end

plot(time,RES)
ylim([-1,110])
title('Biomass vs Time')
xlabel('Time(years)')
ylabel('Biomass')



















