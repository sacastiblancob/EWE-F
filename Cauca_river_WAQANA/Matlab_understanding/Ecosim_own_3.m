%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Ecosim
%
% SECOOOOOOND APPROACH, WITH VULNERABILITIES AND INCLUDING Q/B USER
% LIMITATION IN CONSUMPTION
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

% Max rel P/B
es_rel_PoB_max = es_conf(:,1);

%Feeding time
es_Ftime = es_conf(:,2);

%Max feeding time
es_Ftime_max = es_conf(:,3);

%Feeding time adjust rate  [0,1]
es_Ftime_adjust = es_conf(:,4);

%Fraction of other mortality sensible to changes in feeding time
es_M0_pred = es_conf(:,5);

%Predator effect on feeding time [0,1]
es_risk_time = es_conf(:,6);
es_risk_time(es_risk_time==-999) = 0;

%Density-Dependent catchability Q_max / Q_0 [>= 1]
es_Q_maxoQ_0 = es_conf(:,7);

%Q/B_max / Q/B_0 for handling time [>= 1]
es_QB_maxoQB_0 = es_conf(:,8);

%Switching power parameter [0,2]
es_switch_power = es_conf(:,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting time and initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tiempo
% to = 0;               %segundos
% tf = 86400*90;         %segundos
% dt = 60;              %segundos
to = 0;               %años
tf = 2;         %años
dt = 0.00001;              %años
time = to:dt:tf;
tsteps = length(time);

%Para guardar resultados
RES = zeros(length(ep_biomass),tsteps);
RESftime = zeros(length(ep_biomass),tsteps);

%initial storing
k=1;
RES(:,k) = ep_biomass;
k=k+1;

h = zeros(nvars,1);
f = zeros(nvars,1);

%%%%%%%% PRED %%%%%%%%
% pred is different if mulstistanza
pred = ep_biomass;

%%% calculate risk time for consumers
CB_base = ep_EatenBy./pred;
CB_last = CB_base;
Q_main = (1 - es_risk_time).*CB_base;
Q_risk = es_risk_time.*CB_base.*(ep_EatenOf./ep_biomass + (1-ep_EE).*ep_PoB + 0.0000000001);

% Computing f(Bi) for producers
for i=1:nvars
    if (ep_org_type(i) == 1)
        % h_i = ((r_i*(P/B)_i) - 1)/Bi
        h(i) = (es_rel_PoB_max(i)*ep_PoB(i) - 1)/ep_biomass(i);
        % f_i = r_i*B_i/(1 + B_i*h_i)
        f(i) = (es_rel_PoB_max(i)*ep_biomass(i))/(1 + ep_biomass(i)*h(i));
    end
end

%
% Calculate Lotka-Volterra base parameters
%
% computing a and v
Dzero = es_QB_maxoQB_0./(es_QB_maxoQB_0 - 1); %% For include Q/B user limitation
a = zeros(nvars,numpred);
v = zeros(nvars,numpred);
% Computing a_ij
for i=1:nvars
    for j=1:numpred
        if (ep_diet(i,j)~= 0)
            v(i,j) = (ep_QoB(j)*ep_biomass(j)*ep_diet(i,j)*vul(i,j))/ep_biomass(i);
            den = (v(i,j)*ep_biomass(i)*pred(j) - ...
              ep_QoB(j)*ep_biomass(j)*ep_diet(i,j)*pred(j));
            if den == 0
                den = 1E-10;
            end
            a(i,j) = (2*Dzero(j)*ep_QoB(j)*ep_biomass(j)*ep_diet(i,j)*v(i,j))/den;
            
        else
            a(i,j) = 0;
        end
        
    end
end

%htime = Pred_i / ((Q/B max / Q/B zero) * B_i * (Q/B)_i)
htime = pred ./ (es_QB_maxoQB_0.*ep_biomass.*ep_QoB);
htime(ep_org_type~=2) = 0;
%
% End Calculate Lotka-Volterra base parameters
%
%
% Initialise relative switching parameters
%
%Pred_den
pred_den = zeros(nvars,1);
for j= 1:numpred
   pred_den(j) = sum(a(:,j).*(ep_biomass.^es_switch_power(j)));
end

%base_time_switch

base_time_switch = zeros(nvars,numpred);
for j=1:numpred
    for i=1:nvars
        base_time_switch(i,j) = a(i,j)*(ep_biomass(i)^es_switch_power(j))/(pred_den(j)+1E-20);
    end
end

%
% End Initialise relative switching parameters
%

%Hden (Q/B max / Q/B zero)/((Q/B max / Q/B zero) + 1)
hden0 = (es_rel_PoB_max./ep_QoB)./((es_rel_PoB_max./ep_QoB) + 1);
hden0(ep_QoB==-999) = 0;
hden = hden0;

%dwe parameter
dwe = 0.5;


C = zeros(nvars,numpred);
C2 = zeros(nvars,numpred);
for t = to+dt:dt:tf
    
    %Update pred, pred is different if multistanza
    pred = ep_biomass;

    %
    % Set relative switching parameters
    %
    %Pred_den
    pred_den = zeros(nvars,1);
    for j= 1:numpred
       pred_den(j) = sum(a(:,j).*(ep_biomass.^es_switch_power(j)));
    end

    %base_time_switch

    rela_switch = zeros(nvars,numpred);
    for j=1:numpred
        for i=1:nvars
            if (base_time_switch(i,j)~=0)
            rela_switch(i,j) = a(i,j)*(ep_biomass(i)^es_switch_power(j))/...
                (pred_den(j)+1E-20)/base_time_switch(i,j);
            end
        end
    end
    
    %
    % End Set relative switching parameters
    %    
    
    %
    % calculate Vulnerable Biomasses
    %
    
    %Effective a
    a_eff = a.*rela_switch.*(es_Ftime(1:numpred)');
        
    %Effective v
    v_eff = v.*(es_Ftime(1:numpred)');
    
    %v_denom parameter
    v_denom = a_eff.*(pred(1:numpred)')./(hden(1:numpred)');
        
    %Vulnerable biomass
    v_biom = (v_eff.*ep_biomass)./(v + v_eff + v_denom);
    v_biom((v + v_eff + v_denom)==0) = 0;
    
    %hdent
    hdent = sum(a_eff.*v_biom,1);
    hdent = [hdent';zeros(nvars-numpred,1)];
    
    %hden
    hden = (1-dwe).*(1 + htime.*hdent) + dwe.*hden;
    
    %v_denom parameter
    v_denom = a_eff.*(pred(1:numpred)')./(hden(1:numpred)');
        
    %Vulnerable biomass
    v_biom = (v_eff.*ep_biomass)./(v + v_eff + v_denom);
    v_biom((v + v_eff + v_denom)==0) = 0;
    
    %
    % end calculate Vulnerable Biomasses
    %
    
    %
    % Computing consumption
    %
    for i=1:nvars
        for j=1:numpred
            if(ep_diet(i,j) ~= 0)
                denom = (v(i,j) + v_eff(i,j) + a_eff(i,j)*pred(j)/hden(j))*hden(j); 
                C(i,j) = a_eff(i,j)*v_eff(i,j)*ep_biomass(i)*pred(j)/denom;
                %C2(i,j) = a_eff(i,j)*v_biom(i,j)*pred(j)/hden(j);
                %C(i,j) = a_eff(i,j)*v_biom(i,j)*pred(j)/hden(j);
            else
                C(i,j) = 0;
            end
        end
    end
    %
    % Computing consumption as C_ij = a_ij*B_i*B_j
    %
    %
    % Updating Switching relative parameters
    %
    % Computing EatenBy (the consumption of group in time)
    ep_EatenBy = [sum(C),zeros(1,nvars-numpred)];
    ep_EatenBy = ep_EatenBy';
    
    % Computing EatenOf (the consumption over group in time)
    ep_EatenOf = sum(C,2);
    
    % Computing risk_rate
    risk_rate = (ep_EatenOf./ep_biomass) + (1 - ep_EE).*ep_PoB + 0.0000000001;
    risk_rate(ep_PoB==-999) = 0.0000000001;
    
    % Computing CB_last y Q_opt
    CB_last = ep_EatenBy./pred;
    Q_opt = Q_main + (Q_risk/risk_rate);
    
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
    
    RHSftime = zeros(nvars,1);
    for i=1:nvars
        if (ep_org_type(i) == 2)
            RHSftime(i) = 0.1*(1-es_Ftime(i)) - ...
                0.9*es_Ftime(i)*es_Ftime_adjust(i)*(1 + Q_opt(i)/CB_last(i));
%             RHSftime(i) = es_Ftime(i)*es_Ftime_adjust(i)*(Q_opt(i)/CB_last(i) - 1);
        else
            RHSftime(i) = 0;
        end
    end
    RHSftime = RHSftime*12; %1/month to 1/year
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ep_biomass = ep_biomass + RHS*dt;
    es_Ftime = es_Ftime + RHSftime*dt;
    for i=1:numpred
        if es_Ftime(i) > es_Ftime_max(i)
            es_Ftime(i) = es_Ftime_max(i);
        end
    end
    
    %ep_biomass(ep_biomass < 0) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    RES(:,k) = ep_biomass;
    RESftime(:,k) = es_Ftime;
    k = k+1;
    
end


plot(time,RES)
ylim([-1,110])
title('Biomass vs Time')
xlabel('Time(years)')
ylabel('Biomass')



















