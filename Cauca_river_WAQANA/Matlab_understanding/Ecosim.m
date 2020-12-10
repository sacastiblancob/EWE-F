%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Ecosim
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECOSIM DYNAMIC SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: This code only computes trophic values for consumers
%
% Equations
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
%   Aij: rate of effective search for pool type i by predator j
%   Vij,V'ij, with default setting Vij=V'ij are prey behavioral exchange
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
%   BASE ON PROFESSOR Ekin Agoklu WORK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Globals

relax = 0;

% Reading Ecopath Data
ep_data = h5read('Ecosim_data/Lab.h5','/ep_data');
%ep_detfate = h5read('Ecosim_data/Lab.h5','/ep_detfate');
ep_detfate = csvread('Ecosim_data/Lab_DetFate.csv',4,1);
%ep_detfate = ep_detfate';
% ep_diet = h5read('Ecosim_data/Lab.h5','/ep_diet');
ep_diet = csvread('Ecosim_data/Lab_DC.csv',3,1);
% ep_diet = ep_diet';
ms_data = h5read('Ecosim_data/Lab.h5','/ms_data');

% % ep_data = h5read('Ecosim_data/Tampa_Bay.h5','/ep_data');
% % ep_detfate = h5read('Ecosim_data/Tampa_Bay.h5','/ep_detfate');
% % ep_detfate = ep_detfate';
% % ep_diet = h5read('Ecosim_data/Tampa_Bay.h5','/ep_diet');
% % ep_diet = ep_diet';
% % ms_data = h5read('Ecosim_data/Tampa_Bay.h5','/ms_data');

v = fieldnames(ep_data);
for i = 1:length(v)
    name = v{i};
    myVar = ep_data.(v{i});
    assignin('base',strcat('ep_',name), myVar)
end

v = fieldnames(ms_data);
for i = 1:length(v)
    name = v{i};
    myVar = ms_data.(v{i});
    assignin('base',strcat('ms_',name), myVar)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Imputs and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading imput data
% % nvars = csvread('Ecosim_data/Tampa_Bay_Scenario.csv',0,1,[0,1,0,1]);
nvars = csvread('Ecosim_data/Lab_Scenario.csv',0,1,[0,1,0,1]);
% % nstanzas = csvread('Ecosim_data/Tampa_Bay_Scenario.csv',2,1,[2,1,2,1]);
nstanzas = csvread('Ecosim_data/Lab_Scenario.csv',2,1,[2,1,2,1]);

%Configutarion

%Base proportion of free nutrients (if==1, means free nutrients is infinity)
Nbasefree = 0.9;

%Maximum P/B rate due to the nutrient concentration
NutPBmax = 1.5;

%Scenario
es_conf = csvread('Ecosim_data/Lab_Scenario.csv',4,1);
% % es_conf = csvread('Ecosim_data/Tampa_Bay_Scenario.csv',4,1);
confnames = ["MaxrelPB","Feed_time","Maxrel_feeding_time","Feed_time_adjust_rate",...
    "Fraction_of_other_mortality","Predator_effect_on_feeding_time","Density-dep_catchability_Qmax/Qo",...
    "QBmax/QBo","Switching_power_parameter","Advected"];
%groups = ["GC_D_0_3","GC_D_3_6","GC_C","GC_B","GC_A","GP_A","GD_B","GD_A"];

es_rel_PoB_max = es_conf(:,1);
es_Ftime = es_conf(:,2);
es_Ftime_max = es_conf(:,3);
es_Ftime_adjust = es_conf(:,4);
es_M0_pred = es_conf(:,5);
es_risk_time = es_conf(:,6);
es_Q_maxoQ_0 = es_conf(:,7);
es_QB_maxoQB_0 = es_conf(:,8);
es_switch_power = es_conf(:,9);

%Vulnerability matrix
numpred = csvread('Ecosim_data/Lab_vul.csv',1,1,[1,1,1,1]);
% % numpred = csvread('Ecosim_data/Tampa_Bay_vul.csv',1,1,[1,1,1,1]);
vul = csvread('Ecosim_data/Lab_vul.csv',3,1);
% % vul = csvread('Ecosim_data/Tampa_Bay_vul.csv',3,1);

%%%%%%% SECTION: MAIN CALCULATIONS BEFORE RUN MODEL !!!!!!!!!!!!!!!!!!

ndetritus = sum(ep_org_type == 0);

j = 0;
detritus_no = zeros(ndetritus);
for i = 1:nvars
  if (ep_org_type(i) == 0)
      j = j + 1;
      detritus_no(j) = i;
  end
end

%Calculate Maximum PoB related values
es_abs_PoB_max = zeros(nvars,1);
for i=1:nvars
    if (ep_org_type(i) == 0) || (ep_org_type(i) == 2)
        es_abs_PoB_max(i) = 0;
    else
        es_abs_PoB_max(i) = es_rel_PoB_max(i) * ep_PoB(i);
    end
end

es_PoB_biomass = zeros(nvars,1);
for i=1:nvars
    if (ep_org_type(i) == 1)
        es_PoB_biomass(i) = (es_abs_PoB_max(i)/ep_PoB(i) - 1)/ep_biomass(i);
    else
        es_PoB_biomass(i) = 0;
    end
end

%Calculate Nutrient Concentrations
%Nutrient biomass
NutBiom = sum(ep_biomass);
NutTot = NutBiom/(1 - Nbasefree);

% Base concentration of free nutrients
NutFreeBase = zeros(nvars,1);
NutFree = NutTot - NutBiom;
for i=1:nvars
    if (ep_org_type(i) ~= 0)
        NutFreeBase(i) = NutFree;
    end
end

NutMin = 0.00101*NutFree;

% Remove Import from Diet
es_QBoutside = zeros(nvars);
for j = 1:nvars
    if (ep_org_type(j) == 2)
        if (ep_diet(nvars+1,j) > 0)
            fractionW0import = (1 - ep_diet(nvars+1,j)/1);
        else
            fractionW0import = 1;
        end

        for i = 1:nvars
            if (fractionW0import == 0)
                ep_diet(i,j) = 0;
            else
                ep_diet(i,j) = ep_diet(i, j);
            end
        end
        es_QBoutside(j) = ep_QoB(j) * (1 - fractionW0import);
    end
end
% end Remove Import from Diet

es_hden = zeros(nvars,1);
for i = 1:nvars
    if (ep_org_type(i) == 2)
        es_Ftime(i) = 1;
        es_hden(i)  = es_QB_maxoQB_0(i) / (es_QB_maxoQB_0(i) + 1);
    end 
end

es_CB_base = zeros(nvars,1);
es_CB_last = zeros(nvars,1);
for i = 1:nvars
    es_CB_base(i) = ep_EatenBy(i) / ep_biomass(i);
    if (es_CB_base == 0)
        es_CB_base(i) = 1;
        es_Ftime_max(i) = 1;
    end
    es_CB_last(i) = es_CB_base(i);
end

% % Initialize stanza parameters if any
B = zeros(nvars,1);

if nstanzas > 0
%
    es_ms_NageS = zeros(max(ms_age_infinity)+1,nstanzas);
    es_ms_WageS = zeros(max(ms_age_infinity)+1,nstanzas);
    es_ms_EggsSplit = zeros(max(ms_age_infinity)+1,nstanzas);
    es_ms_BaseEggsStanza = zeros(nstanzas);
    es_ms_EggsStanza = zeros(nstanzas);
    es_ms_RscaleSplit = zeros(nstanzas);
    es_pred = zeros(nvars,1);
    es_ms_SplitAlpha = zeros(max(ms_age_infinity)+1,nstanzas);
% initialise Split Groups
    for stanza = 1:nstanzas
        Be = 0;
        for age = 1:ms_age_infinity(stanza)+1
            es_ms_NageS(age,stanza) = ms_splitno(age,stanza);
            es_ms_WageS(age,stanza) = ms_Wage(age,stanza);
            
            if (es_ms_WageS(age,stanza) > ms_Wmat_Winf(stanza))
                
                es_ms_EggsSplit(age,stanza) = es_ms_WageS(age,stanza) - ...
                    ms_Wmat_Winf(stanza);
                
                Be = Be + es_ms_NageS(age,stanza)*es_ms_EggsSplit(age,stanza);
                
            end
        end
        es_ms_BaseEggsStanza(stanza) = Be;
        es_ms_EggsStanza(stanza) = Be;
        es_ms_RscaleSplit(stanza) = 1;
    end
    
%   Set split Pred
    for stanza = 1:nstanzas
        for substanza = 1:ms_substanzas(stanza)
            Bt = 1E-30;
            Pt = 1E-30;
            Nt = 1E-30;
            
            if (substanza < ms_substanzas(stanza))
                age_last = ms_age_start(substanza+1,stanza) - 1;
            else
                age_last = ms_age_infinity(stanza);
            end
            
            for age = ms_age_start(substanza,stanza)+1:age_last+1
                Bt = Bt + es_ms_NageS(age,stanza)*es_ms_WageS(age,stanza);
                Pt = Pt + es_ms_NageS(age,stanza)*ms_WWa(age,stanza);
                Nt = Nt + es_ms_NageS(age,stanza);
            end
            
            B(ms_ep_groupno(substanza,stanza)) = Bt;
            es_pred(ms_ep_groupno(substanza,stanza)) = Pt;
        end
    end
%   end Set Split Pred
    
    for stanza = 1:nstanzas
        for substanza = 1:ms_substanzas(stanza)
            if (substanza < ms_substanzas(stanza))
                Agem = ms_age_start(substanza+1,stanza) - 1;
            else
                Agem = ms_age_infinity(stanza);
            end
            
            if substanza == ms_substanzas(stanza)
                Agem = ms_age_infinity(stanza);
            end
            
            for age = ms_age_start(1,stanza) + 1:Agem + 1
                es_ms_SplitAlpha(age,stanza) = (ms_Wage(age+1,stanza) - ...
                    ms_vbM(stanza)*ms_Wage(age,stanza))*...
                    es_pred(ms_ep_groupno(substanza,stanza))/...
                    ep_EatenBy(ms_ep_groupno(substanza,stanza));
            end
        end
        es_ms_SplitAlpha(ms_age_infinity(stanza)+1) = es_ms_SplitAlpha(ms_age_infinity(stanza));
    end
% end initialise Split Groups
end

B = zeros(nvars,1);
%

integrate = zeros(nvars,1);
for i = 1:nvars
    if (ep_isstanza(i) == 1)
        integrate(i) = -i;
    else
        integrate(i) = i;
    end
end

% setpred()
b_pred = ep_biomass;
for i=1:nvars
    if(ep_biomass(i) < 1E-20)
        b_pred(i) = 1E-20;
    end
    if(integrate(i) >= 0)
        es_pred(i) = b_pred(i);
    end
end

for i = 1:nvars
    if (es_risk_time == -999)
        es_risk_time(i) = 0;
    end
end

% calculate risk time for consumers
es_Q_main = zeros(nvars,1);
es_Q_risk = zeros(nvars,1);
for i=1:nvars
    if(ep_org_type(i)==2)
        es_CB_base(i) = ep_EatenBy(i)/es_pred(i);
        es_CB_last(i) = es_CB_base(i);
        es_Q_main(i) = (1 - es_risk_time(i))*es_CB_base(i);
        es_Q_risk(i) = es_risk_time(i)*es_CB_base(i)*(ep_EatenOf(i)/ep_biomass(i) + ...
            ((1 - ep_EE(i))*ep_PoB(i)) + 0.0000000001);
        
    
    end
end

% Calculate Lokta Volterra Effective Search Rates
es_htime=zeros(nvars,1);
for i=1:nvars
    if (ep_org_type(i)==2)
        es_htime(i) = es_pred(i)/(es_QB_maxoQB_0(i)*ep_biomass(i)*ep_QoB(i));
    else
        es_htime(i) = 0;
    end
end

vsize = size(vul);
arena_vulrate = zeros(vsize);
arena_a = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        
        Dzero = es_QB_maxoQB_0(pred)/(es_QB_maxoQB_0(pred) - 1);
        
        if (vul(prey,pred) ~= -999)
            Consumption = ep_biomass(pred)*ep_QoB(pred)*ep_diet(prey,pred);
            arena_vulrate(prey,pred) = vul(prey,pred)*Consumption/ep_biomass(prey);
            
            % Denominator for computing search rate of predator
            Denv = ep_biomass(prey)*es_pred(pred)*arena_vulrate(prey,pred) - ...
                Consumption*es_pred(pred);
            
            if Denv < 1E-20
                Denv = 1E-20;
            end
            
            arena_a(prey,pred) = Dzero*2*Consumption*arena_vulrate(prey,pred)/Denv;
        else
            arena_a(prey,pred) = 0;
        end
    end
end

% Initial switching parameters (InitRelaSwitch)
es_pred_den = zeros(vsize(2),1);
for pred = 1:vsize(2)
    es_pred_den(pred) = 0;
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            es_pred_den(pred) = es_pred_den(pred) + arena_a(prey,pred)*...
                ep_biomass(prey)^es_switch_power(pred);
        end
    end
end

arena_base_time_switch = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            arena_base_time_switch(prey,pred) = arena_a(prey,pred)*ep_biomass(prey)^...
                es_switch_power(pred)/(es_pred_den(pred) + 1E-20);
        end
    end
end

% set arena vulnerability and search rates
arena_Q_arena = zeros(vsize);
arena_Q_link = zeros(vsize);
arena_vul_arena = zeros(vsize);
arena_vul_biom = zeros(vsize);
arena_a_link = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey, pred) ~= -999)

            % total consumptions of predators in arena
            arena_Q_link(prey, pred) = (ep_biomass(pred)* ep_QoB(pred) *...
                ep_diet(prey, pred));

            arena_Q_arena(prey, pred) =  arena_Q_arena(prey, pred) +...
              arena_Q_link(prey, pred);

            % setting of initial vulnerable biomasses
            arena_vul_arena(prey, pred) = (vul(prey, pred) +...
              0.0000000001) * arena_Q_arena(prey, pred) / ep_biomass(prey);

            if (arena_vul_arena(prey, pred) == 0)
                arena_vul_arena(prey, pred) = 1;
            end

            arena_vul_biom(prey, pred) = (vul(prey, pred) +...
              0.0000000001D0 - 1.0D0) * arena_Q_arena(prey, pred)/...
              (2 * arena_vul_arena(prey, pred));

            if (arena_vul_biom(prey, pred) == 0)
                arena_vul_biom(prey, pred) = 1;
            end
      
            % setting of predator search rates
            if (arena_vul_biom(prey, pred) > 0)
                
                % QBmaxQBo / (QbmaxQBo - 1)
                Dzero = es_QB_maxoQB_0(pred)/(es_QB_maxoQB_0(pred) - 1);

                arena_a_link(prey, pred) = Dzero * arena_Q_link(prey, pred) / ... 
                  (arena_vul_biom(prey, pred) * es_pred(pred));
            else
                arena_a_link(prey, pred) = 0;
            end
        end
    end
end

% % end set Arena Vulnerability and Search Rates

% % prepare for run model by calculating initial values
% % for preparing initial values, year is set to null
% % to disregard forcing data

es_fishmort = zeros(nvars,1);
for j = 1:nvars
  es_fishmort(j) = ((ep_landings(j) + ep_discards(j))/ ep_biomass(j));
end

FirstTime = 1;

time = 0.0;
biomass = ep_biomass;
xdot = zeros(nvars,1);
biomeq = zeros(nvars,1);
loss = zeros(nvars,1);

% % hasta aca todo va bien (09-12-2020)

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %DERIV
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocatea for derivs
integsum = zeros(nvars,1);
es_PoB_base = zeros(nvars,1);
es_abs_PoB_max = zeros(nvars,1);
es_PoB_biomass = zeros(nvars,1);
arena_rela_switch = zeros(vsize);
es_pred_den = zeros(nvars,1);
arena_base_time_switch = zeros(vsize);
es_hdent = zeros(nvars,1);
arena_v_denom = zeros(vsize);
arena_a_eff = zeros(vsize);
arena_v_eff = zeros(vsize);
arena_v_biom = zeros(vsize);
es_pp = zeros(nvars,1);
es_qq = zeros(nvars,1);
es_EatenBy = zeros(nvars,1);
es_unassimilated = zeros(nvars,1);
es_M0 = zeros(nvars,1);
es_M2 = zeros(nvars,1);
es_EatenOf = zeros(nvars,1);
flow2detritus = zeros(ndetritus,1);
det_export_rate = zeros(ndetritus,1);
%
% imonth
imonth = 0;
%
derivs
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %END DERIV
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FirstTime = 1;

%%%%% this is the rate of sinks in each state variable in the derivs()
rrate = zeros(nvars,1);
for i = 1:nvars
    rrate(i) = abs(loss(i)) / ep_biomass(i);
end

%%%%% determine whether to integrate or not each state variable
%%%%% depending on the rate of change in one time step
for i = 1:nvars
    if (rrate(i) > 24 && integrate(i) == i)
        integrate(i) = 0;
        % else if (ep_data(i)%org_type == 0) then
        %     integrate(i) = 0
    elseif (ep_isstanza(i) == 1)
        integrate(i) = -i;
    else

    end
end

%%%%% Init is done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% !!!!!!!!!!!!!!!!!!!!!! Init is done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
%   ! Now calculate final rate of change to set integration flags
%   call calculateMaximumPoBRelatedValues ()
%
for var = 1:nvars
    if (ep_org_type(var) == 0 || ep_org_type(var) == 2)
        es_abs_PoB_max(var) = 0;
    else
        es_abs_PoB_max(var) = es_rel_PoB_max(var) * ep_PoB(var);
    end
end

for var = 1:nvars
    if (ep_org_type(var) == 1)
        es_PoB_biomass(var) = (es_abs_PoB_max(var)/ep_PoB(var) - 1) / ep_biomass(var);
    else
        es_PoB_biomass(var) = 0;
    end
end
%
%   end calculateMaximumPoBRelatedValues ()

time = 1.0;
biomass = ep_biomass;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %DERIV
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
derivs
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %END DERIV
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FirstTime = 1.0;

%%%% this is the rate of sinks in each state variable in the derivs()
for i = 1:nvars
  rrate(i) = abs(loss(i))/ep_biomass(i);
end

%%%%% determine whether to integrate or not each state variable
%%%%% depending on the rate of change in one time step
for i = 1:nvars
    if (rrate(i) > 24 && integrate(i) == i)
          integrate(i) = 0;
          % % else if (ep_data(i)%org_type == 0) then
          % %     integrate(i) = 0
    elseif (ep_isstanza(i) == 1)
        integrate(i) = -i;
    else

    end
end

%%%%%%%%%%%% SECTION: RUN MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% hasta aca todo va bien

to = 0;
tf = 2;
% StepsPerMonth = 60*24*30;
StepsPerMonth = 1;
tstep = 1/(12*StepsPerMonth);
noftsteps = tf/tstep;

mat_out = zeros(noftsteps + 1, nvars);
mat_out(1,:) = ep_biomass;

catch_out_monthly = zeros((12 * tf) + 1, nvars);
for j = 1:nvars
    catch_out_monthly(1, j) = ep_landings(j) + ep_discards(j);
end

BB = ep_biomass;

time = 0;
imonth = 0;
step = 0;

%%%%% run model over the specified time frame
for i = 0:(tf - 1)
  for m = 1:12
%i=0;
%m=1;

      % % Clean monthly stanza variables
      BBAvg = zeros(nvars,1);
      LossAvg = zeros(nvars,1);
      EatenByAvg = zeros(nvars,1);
      EatenOfAvg = zeros(nvars,1);
      PredAvg = zeros(nvars,1);

      imonth  = (i * 12) + m;

      % % call calculateFishingMortalities (BB)
      es_Q_mult = zeros(nvars,1);  
      for var = 1:nvars
            if (es_Q_maxoQ_0(var) ~= -999)
                es_Q_mult(var) = es_Q_maxoQ_0(var)/(1 + (es_Q_maxoQ_0(var) - 1) *...
                  BB(var)/ep_biomass(var));
            else
                es_Q_mult(var) = 0;
            end
%             if (force%forcetype(var) == 5 .and. force%fishforce(imonth, var) /= -999) then
%                 es_data(var)%fishmort = ((ep_data(var)%landings + ep_data(var)%discards) &
%                   / ep_data(var)%biomass) * es_data(var)%Q_mult * force%fishforce(imonth, var)
%             elseif (force%forcetype(var) == 4 .and. force%fishforce(imonth, var) /= -999) then
%                 es_data(var)%fishmort = force%fishforce(imonth, var) &
%                   * es_data(var)%Q_mult
%             elseif  (force%fishforce(imonth, var) == -999) then
%                 es_data(var)%fishmort = ((ep_data(var)%landings + ep_data(var)%discards) &
%                   / ep_data(var)%biomass) * es_data(var)%Q_mult
%             end
      end
        % % end call calculateFishingMortalities (BB)
      
      for n = 1:StepsPerMonth
          %n=1;
          time
          if (n == StepsPerMonth)
              UpdateStanzas = 1.0;
          else
              UpdateStanzas = 0.0;
          end

%           ! call the Runge-Kutta 4th order numeric ode solver
%           call rk4 (BB, time, tstep, integrate)
          % % RK4 but actually Euler
          biomass = BB;
          
          % % % call derivs (computing RHS)
          derivs
          % % % end call derivs
          
          BB = biomass;
            
            lossSt = zeros(nvars);
            for k= 1:nvars
                EatenByAvg(k)  = EatenByAvg(k) + es_qq(k);
                EatenOfAvg(k)  = EatenOfAvg(k) + es_M2(k);
                PredAvg(k)     = PredAvg(k) + es_pred(k);
                lossSt(k)      = loss(k);
            end
          
            for k = 1:nvars
                if (integrate(k) == 0)
                    BB(k) = (1 - relax)*biomeq(k) + relax*BB(k);
                    % yt(k) = (1 - relax) * biomeq(i) + relax * B(i);
                elseif (integrate(k) ~= 0)
                    BB(k) = BB(k) + tstep*xdot(k);
                    % yt(i) = B(i) + dh * dydx(i);
                else
                    
                end
%                 BB(k) = BB(k) + tstep*xdot(k);
            end
            
            % averaging here for multistanza calculations
            % and updating foraging times
            for k = 1:nvars

                BBAvg(k)   = BBAvg(k) + BB(k);
                LossAvg(k) = LossAvg(k) + loss(k);

                if (UpdateStanzas == 1.0)
                    BBAvg(k)        = BBAvg(k) / StepsPerMonth;
                    LossAvg(k)      = LossAvg(k) / StepsPerMonth;
                    EatenByAvg(k)   = EatenByAvg(k) / StepsPerMonth;
                    EatenOfAvg(k)   = EatenOfAvg(k) / StepsPerMonth;
                    PredAvg(k)      = PredAvg(k) / StepsPerMonth;
                    lossSt(k)       = lossSt(k) / StepsPerMonth;
                end

            end

            if (UpdateStanzas == 1.0)

                if (nstanzas > 0)
                    % % below was lossSt but corrected upon EwE6 rk4 routine check
                    % % call multistanza (B)
                    for stanza = 1:nstanzas

                        Be = 0;
                        for substanza = 1:ms_substanzas(stanza)

                            Su = exp(-LossAvg(ms_ep_groupno(substanza,stanza)) / 12 /...
                              BBAvg(ms_ep_groupno(substanza,stanza)));

                            Gf = es_EatenBy(ms_ep_groupno(substanza,stanza)) /...
                              es_pred(ms_ep_groupno(substanza,stanza));

                            if (substanza < ms_substanzas(stanza))
                                age_last = ms_age_start(substanza+1,stanza) - 1;
                            else
                                age_last = ms_age_infinity(stanza);
                            end

                            for age = ms_age_start(substanza,stanza)+1:age_last+1
                                es_ms_NageS(age,stanza) = es_ms_NageS(age,stanza)*Su;
                                es_ms_WageS(age,stanza) = ms_vbM(stanza)*es_ms_WageS(age,stanza) +...
                                  Gf*es_ms_SplitAlpha(age,stanza);

                                if (es_ms_WageS(age,stanza) > ms_Wmat_Winf(stanza))
                                    Be = Be + es_ms_NageS(age,stanza)*(es_ms_WageS(age,stanza) -...
                                      ms_Wmat_Winf(stanza));
                                end
                            end 
                        end

                        es_ms_WageS(ms_age_infinity(stanza),stanza)=(Su*ms_Ahat(stanza) +...
                          (1 - Su)*es_ms_WageS(ms_age_infinity(stanza) - 1,stanza)) /...
                          (1 - ms_Rhat(stanza)*Su);

                        es_ms_EggsStanza(stanza) = Be;

                        for substanza = ms_substanzas(stanza):-1:1

                            if (substanza == ms_substanzas(stanza))
                                AgeMax = ms_age_infinity(stanza);
                            else
                                AgeMax = ms_age_start(substanza + 1,stanza) - 1;
                            end

                            if (substanza > 1)
                                AgeMin = ms_age_start(substanza,stanza);
                            else
                                AgeMin = 1;
                            end

                            if (substanza == ms_substanzas(stanza))
                                % Nt = es_ms_NageS(AgeMax,stanza)+es_ms_NageS(AgeMax - 1,stanza);
                                Nt = es_ms_NageS(AgeMax+1,stanza)+es_ms_NageS(AgeMax,stanza);
                                if (Nt == 0)
                                    Nt = 1E-30;
                                end

                                %es_ms_NageS(AgeMax,stanza) = Nt;
                                es_ms_NageS(AgeMax+1,stanza) = Nt;
                                AgeMax = AgeMax - 1;
                            end

                            for age = AgeMax+1:-1:AgeMin+1
                                es_ms_NageS(age,stanza)=es_ms_NageS(age - 1,stanza);
                                es_ms_WageS(age,stanza)=es_ms_WageS(age - 1,stanza);
                            end

                        end

                        if (es_ms_BaseEggsStanza(stanza) > 0)
                            % es_ms_NageS(ms_age_start(1,stanza),stanza)=es_ms_RscaleSplit(stanza)*...
                            %  ms_RzeroS(stanza);
                            es_ms_NageS(ms_age_start(1,stanza)+1,stanza)=es_ms_RscaleSplit(stanza)*...
                              ms_RzeroS(stanza);
                        end

                        % es_ms_NageS(ms_age_start(1,stanza),stanza) = ...
                        %   es_ms_NageS(ms_age_start(1,stanza),stanza)*(es_ms_EggsStanza(stanza)/...
                        %   es_ms_BaseEggsStanza(stanza))^ms_rec_power(stanza);
                        es_ms_NageS(ms_age_start(1,stanza)+1,stanza) = ...
                          es_ms_NageS(ms_age_start(1,stanza)+1,stanza)*(es_ms_EggsStanza(stanza)/...
                          es_ms_BaseEggsStanza(stanza))^ms_rec_power(stanza);

                        % es_ms_WageS(ms_age_start(1,stanza),stanza) = 0;
                        es_ms_WageS(ms_age_start(1,stanza)+1,stanza) = 0;
                    end
                end

                % Update foraging times at the end of each time step
                % % call updateForagingTimes (integrate)
                es_risk_rate = zeros(nvars);
                es_Q_opt = zeros(nvars);
                for var = 1:nvars
                    if (ep_org_type(var) == 2)
                        if (integrate(var) == var || integrate(var) < 0)
                            es_CB_last(var) = EatenByAvg(var) / PredAvg(var);
                        end

                        if (ep_PoB(var) ~= -999)
                            es_risk_rate(var) = EatenOfAvg(var) / BBAvg(var) ...
                                + ((1 - ep_EE(var)) * ep_PoB(var)) + 0.0000000001;
                        else
                %            es_risk_rate(var)% = real(0.0000000001D0, 4)
                            es_risk_rate(var) = 0.0000000001;
                        end

                        es_Q_opt(var) = es_Q_main(var) + es_Q_risk(var)/es_risk_rate(var);

                        if (es_CB_last(var) > 0 && (integrate(var) == var || integrate(var) < 0))

                            es_Ftime(var) = 0.1D0 + 0.9D0 * es_Ftime(var)*...
                              (1 - es_Ftime_adjust(var) +...
                              es_Ftime_adjust(var) * es_Q_opt(var) /...
                              es_CB_last(var));

                            if (es_Ftime(var) > es_Ftime_max(var))
                                es_Ftime(var) = es_Ftime_max(var);
                            end
                        end
                    end
                end

            end
          % % end of RK4 but actually Euler
          mat_out(step + 2, :) = BB;


% %           ! calculate relative change
% %           ! with respect to initial biomasses
% %           for j = 1, nvars
% %               rel_out(step + 2, j) = mat_out(step + 2, j) &
% %                    / ep_data(j)%biomass
% %           end

          step = step + 1;
          time = time + tstep;
      end

% %       mat_out_monthly(imonth + 1, :) = BB;

      for var = 1:nvars
          catch_out_monthly(imonth + 1, var) = BB(var)*es_fishmort(var);
      end

% %       ! calculate relative change with respect to initial biomasses
% %       for j = 1, nvars
% %           rel_out_monthly(imonth + 1, j) &
% %                = mat_out_monthly(imonth + 1, j) / ep_data(j)%biomass
% %       end
  end

end

times = to:tstep:tf;
times = times';
plot(times,mat_out)
legend('1','2','3','4','5','6','7','8');






























