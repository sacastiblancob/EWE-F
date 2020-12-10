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

Globals

% Reading Ecopath Data
ep_data = h5read('Ecosim_data/Lab.h5','/ep_data');
ep_detfate = h5read('Ecosim_data/Lab.h5','/ep_detfate');
ep_detfate = ep_detfate';
ep_diet = h5read('Ecosim_data/Lab.h5','/ep_diet');
ep_diet = ep_diet';
ms_data = h5read('Ecosim_data/Lab.h5','/ms_data');

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
nvars = csvread('Ecosim_data/Lab_Scenario.csv',0,1,[0,1,0,1]);
nstanzas = csvread('Ecosim_data/Lab_Scenario.csv',2,1,[2,1,2,1]);

%Configutarion

%Base proportion of free nutrients (if==1, means free nutrients is infinity)
Nbasefree = 0.99;

%Maximum P/B rate due to the nutrient concentration
NutPBmax = 1.5;

%Scenario
es_conf = csvread('Ecosim_data/Lab_Scenario.csv',4,1);
confnames = ["MaxrelPB","Feed_time","Maxrel_feeding_time","Feed_time_adjust_rate",...
    "Fraction_of_other_mortality","Predator_effect_on_feeding_time","Density-dep_catchability_Qmax/Qo",...
    "QBmax/QBo","Switching_power_parameter","Advected"];
groups = ["GC_D_0_3","GC_D_3_6","GC_C","GC_B","GC_A","GP_A","GD_B","GD_A"];

es_rel_PoB_max = es_conf(:,1);
es_Ftime = es_conf(:,2);
es_Ftime_max = es_conf(:,3);
es_Ftime_adjust = es_conf(:,4);
es_M0_pred = es_conf(:,5);
es_risk_time = es_conf(:,6);
es_Q_maxoQB_0 = es_conf(:,7);
es_QB_maxoQB_0 = es_conf(:,8);
es_switch_power = es_conf(:,9);

%Vulnerability matrix
numpred = csvread('Ecosim_data/Lab_vul.csv',1,1,[1,1,1,1]);
vul = csvread('Ecosim_data/Lab_vul.csv',3,1);

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
    es_ms_NageS = zeros(ms_age_infinity+1,nstanzas);
    es_ms_WageS = zeros(ms_age_infinity+1,nstanzas);
    es_ms_EggsSplit = zeros(ms_age_infinity+1,nstanzas);
    es_ms_BaseEggsStanza = zeros(nstanzas);
    es_ms_EggsStanza = zeros(nstanzas);
    es_ms_RscaleSplit = zeros(nstanzas);
    es_pred = zeros(nvars,1);
    es_ms_SplitAlpha = zeros(ms_age_infinity+1,nstanzas);
% initialise Split Groups
    for stanza = 1:nstanzas
        Be = 0;
        for age = 1:ms_age_infinity+1
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
        es_htime(i) = es_pred(i)/es_QB_maxoQB_0(i)*ep_biomass(i)*ep_QoB(i);
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
derivs



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %DERIV
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % time = 0;
% % biomass = ep_biomass;
% % xdot = zeros(nvars,1);
% % biomeq = zeros(nvars,1);
% % loss = zeros(nvars,1);
% % 
% % % set predator abundance
% % b_derivs = biomass;
% % integsum = zeros(nvars,1);
% % for var = 1:nvars
% %     if (biomass(var) < 1.0e-20)
% %         b_derivs(var) = 1.0e-20;
% %     end
% %     if (integrate(var) >= 0)
% %         es_pred(var) = b_derivs(var);
% %     end
% %     integsum(var) = b_derivs(var);
% % end
% % 
% % % nutrient biomass
% % NutBiom = 0;
% % for var = 1:nvars
% %     NutBiom = NutBiom + biomass(var);
% % end
% % 
% % if (time == 0)
% %     NutFree = NutTot - NutBiom;
% % else
% % % % #ifdef _ForceNutrient_
% % % %     if (time == 1.0) then
% % % %         NutFree = NutTot - NutBiom
% % % %     else
% % % %         NutFree = NutTot * NutrientForce(imonth) - NutBiom
% % % %     end if
% % % % #endif
% %     NutFree = NutTot - NutBiom;
% % end
% % 
% % % amount of free nutrients (NutFree) must not be less than
% % % the background nutrient concentration (NutMin)
% % if (NutFree < NutMin)
% %     NutFree = NutMin;
% % end
% % 
% % j = 0;
% % es_PoB_base = zeros(nvars);
% % for var = 1:nvars
% %     if (ep_org_type(var) ~= 0)
% % % P/B scaled by Max. Relative Production, i.e. es_data%pob_max.
% %         if (time == 0)
% %             Pmult = 1.0;
% %             if (ep_org_type(var) ~= 0)
% %                 es_abs_PoB_max(var) = es_rel_PoB_max(var)* ep_PoB(var);
% %                 es_PoB_biomass(var) = 0;
% %             end
% %         else
% % % % #ifdef _ForcePrimaryProd_
% % % %             if (ep_data(var)%org_type == 1 .and. time /= 1.0) then
% % % %                 j = j + 1
% % % %                 Pmult = PrimaryProdForce(imonth, j)
% % % %             end
% % % % #endif
% % % % #ifndef _ForcePrimaryProd_
% %             Pmult = 1.0;
% % % #endif
% %         end
% %         es_PoB_base(var) = 2*NutFree/(NutFree + NutFreeBase(var))*Pmult*...
% %             es_abs_PoB_max(var)/ (1 + biomass(var) * es_PoB_biomass(var));
% %    end
% % end
% % 
% % % set relative prey switching
% % arena_rela_switch = zeros(vsize);
% % for pred = 1:vsize(2)
% %     es_pred_den(pred) = 0;
% %     for prey = 1:vsize(1)
% %         if (vul(prey,pred) ~= -999)
% %             es_pred_den(pred) = es_pred_den(pred) + arena_a(prey,pred)*...
% %                 biomass(prey)^es_switch_power(pred);
% %             arena_rela_switch(prey,pred)=1;
% %         end
% %     end
% % end
% % 
% % arena_base_time_switch = zeros(vsize);
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         if (vul(prey,pred) ~= -999)
% %             if es_switch_power(pred)>0
% %                 arena_base_time_switch(prey,pred) = arena_a(prey,pred)*ep_biomass(prey)^...
% %                     es_switch_power(pred)/(es_pred_den(pred) + 1E-20)/arena_base_time_switch(prey,pred);
% %             end
% %         end
% %     end
% % end
% % % end set relative prey switching
% % 
% % % Vulnerability Calculations
% % es_hdent = zeros(nvars);
% % arena_v_denom = zeros(vsize);
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         if (vul(prey,pred) ~= -999)
% %             es_hdent(pred) = 0;
% %             arena_v_denom(prey, pred) = 0;
% %         end
% %     end
% % end
% % 
% % dwe = 0.5;
% % arena_a_eff = zeros(vsize);
% % arena_v_eff = zeros(vsize);
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         arena_v_denom(prey, pred) = 0;
% %         if (vul(prey,pred) ~= -999)
% %             arena_a_eff(prey, pred) = arena_a_link(prey, pred) * es_Ftime(pred) *...
% %               arena_rela_switch(prey, pred);
% % 
% %             arena_v_eff(prey, pred) = arena_vul_arena(prey, pred) * es_Ftime(prey);
% % 
% %             arena_v_denom(prey, pred) = arena_v_denom(prey, pred) + arena_a_eff(prey, pred) *...
% %               es_pred(pred) / es_hden(pred);
% %         end
% %     end
% % end
% % 
% % arena_v_biom = zeros(vsize);
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         if (vul(prey,pred) ~= -999)
% %             arena_v_biom(prey, pred) = arena_v_eff(prey, pred) * biomass(prey)/...
% %               (arena_vul_arena(prey, pred) + arena_v_eff(prey, pred) + arena_v_denom(prey, pred));
% %         end
% %     end
% % end
% % 
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         if (vul(prey,pred) ~= -999)
% %             es_hdent(pred) = es_hdent(pred) + arena_a_eff(prey, pred)*arena_v_biom(prey, pred);
% %         end
% %     end
% % end
% % 
% % for pred = 1:vsize(2)
% %     es_hden(pred) = (1 - dwe) * (1 + es_htime(pred)*es_hdent(pred)) + dwe*es_hden(pred);
% % end
% % 
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         arena_v_denom(prey, pred) = 0;
% %         if (vul(prey,pred) ~= -999)
% %             arena_v_denom(prey, pred) = arena_v_denom(prey, pred) + arena_a_eff(prey, pred)*...
% %               es_pred(pred) / es_hden(pred);
% %         end
% %     end
% % end
% % 
% % for pred = 1:vsize(2)
% %     for prey = 1:vsize(1)
% %         if (vul(prey,pred) ~= -999)
% %             arena_v_biom(prey, pred) = arena_v_eff(prey, pred)*biomass(prey)/...
% %               (arena_vul_arena(prey, pred) + arena_v_eff(prey, pred) + arena_v_denom(prey, pred));
% %         end
% %     end
% % end
% % % end Vulnerability Calculations
% % 
% % %%%%% calculate processes; production and consumption
% % % primary production
% % es_pp = zeros(nvars,1);
% % for var = 1:nvars
% %     if (ep_org_type(var) == 1)
% %         es_pp(var) = es_PoB_base(var)*biomass(var);
% %     else
% %         es_pp(var) = 0;
% %     end
% % end
% % 
% % % consumption
% % es_qq = zeros(nvars);
% % es_EatenBy = zeros(nvars);
% % for pred = 1:nvars
% %     if (pred <= vsize(2))
% %         consumption = 0;
% %         for prey = 1:vsize(1)
% %             if (vul(prey, pred) ~= -999)
% %                 consumption = consumption + (arena_a_eff(prey, pred)*...
% %                   arena_v_biom(prey, pred) * es_pred(pred)/es_hden(pred));
% %             end
% %         end
% %         es_qq(pred) = consumption + es_QBoutside(pred)*biomass(pred);
% %         es_EatenBy(pred) = es_qq(pred);
% %     else
% %         es_qq(pred) = 0;
% %         es_EatenBy(pred) = es_qq(pred) / biomass(pred);
% %     end
% % end
% % 
% % % % unassimilated consumption
% % es_unassimilated = zeros(nvars);
% % for var = 1:nvars
% %     if (es_qq(var) ~= 0)
% %         es_unassimilated(var) = es_qq(var) * ep_unass_Q(var);
% %     else
% %         es_unassimilated(var) = 0;
% %     end
% % end
% % 
% % % % non-predation mortalities
% % es_M0 = zeros(nvars);
% % for var = 1:nvars
% %     if (ep_org_type(var) == 1)
% %         es_M0(var) = ((1 - ep_EE(var)) * ep_PoB(var)) * biomass(var);
% %     elseif (ep_org_type(var) == 2)
% %         es_M0(var) = ((1 - ep_EE(var)) * ep_PoB(var))*(1 - es_M0_pred(var) + ...
% %             es_M0_pred(var)* es_Ftime(var)) * biomass(var);
% %     else
% %         es_M0(var) = 0;
% %     end
% % end
% % 
% % % % predation mortalities
% % es_M2 = zeros(nvars);
% % es_EatenOf = zeros(nvars);
% % for prey = 1:nvars
% %     M2_predation = 0;
% %     for pred = 1:vsize(2)
% %         if (vul(prey, pred) ~= -999)
% %             M2_predation = M2_predation + (arena_a_eff(prey, pred)*...
% %               arena_v_biom(prey, pred)*es_pred(pred)/es_hden(pred));
% %         end
% %     end
% %     es_M2(prey) = M2_predation;
% %     es_EatenOf(prey) = M2_predation;
% % end
% % 
% % 
% % % Calculate Detrital Flows
% % flow2detritus = zeros(ndetritus);
% % for j = 1:ndetritus
% %     flow2detritus(j) = 0;
% %     for var = 1:nvars
% %         if (ep_org_type(var) ~= 0) 
% %             % calculate non-predation natural mortality 
% %             % for groups other than detritus
% %             flow2detritus(j) = flow2detritus(j) + ((1 - ep_EE(var))*...
% %               ep_PoB(var))*biomass(var)*ep_detfate(var,j);
% %         end
% % 
% %         if (ep_org_type(var) == 2) 
% %             % if group is a consumer, consider calculating excretion
% %             flow2detritus(j) = flow2detritus(j) + ep_unass_Q(var)*...
% %               es_EatenBy(var)*ep_detfate(var, j);
% %         end
% % 
% %         if (ep_discards(var) ~= 0)
% %             if (FirstTime == 1)
% %                 flow2detritus(j)  = flow2detritus(j) + (ep_landings(var)+...
% %                 ep_discards(var))*(ep_discards(var)/(ep_landings(var) + ep_discards(var)));
% %             else
% %                 % then calculate flows from discarded fish
% %                 flow2detritus(j)  = flow2detritus(j) + biomass(var)*...
% %                   es_fishmort(var)*(ep_discards(var)/(ep_landings(var) + ep_discards(var)));
% %             end
% %         end
% %      end
% % end
% % 
% % %%%% initialize the export rate of detritus out of the system
% % det_export_rate = zeros(ndetritus);
% % if (FirstTime ==1)
% %     for j = 1:ndetritus
% %         det_export_rate(j) = (flow2detritus(j) - es_EatenOf(detritus_no(j))...
% %             + ep_detritus_import(detritus_no(j))) / biomass(detritus_no(j));
% %     end
% % end
% % 
% % FirstTime = 0;
% % % end Calculate Detrital Flows
% % 
% % % sum of sinks for state variables
% % j = 0;
% % for var = 1:nvars
% %     if (ep_org_type(var) ~= 0)
% %         loss(var) = es_M0(var) + es_M2(var) + es_fishmort(var)*biomass(var);
% %     else
% %         j = j + 1;
% %         loss(var) = det_export_rate(j) * biomass(var) + es_EatenOf(var);
% %     end
% % end
% % 
% % % now calculate state equations
% % j = 0;
% % for var = 1:nvars
% %     if (ep_org_type(var) == 2)
% %         xdot(var) = (ep_PoQ(var) * es_qq(var) + biomass(var)*es_PoB_base(var)) - loss(var);
% %         if (loss(var) > 0 && biomass(var) > 0)
% %             biomeq(var) = (ep_PoQ(var)*es_qq(var) + biomass(var)*...
% %               es_PoB_base(var))/(loss(var)/biomass(var));
% %         else
% %             biomeq(var) = 1E-20;
% %         end
% %     elseif (ep_org_type(var) == 1)
% %         xdot(var) = es_pp(var) - loss(var);
% %         if (loss(var) > 0 && biomass(var) > 0)
% %             biomeq(var) = es_pp(var)/(loss(var) / biomass(var));
% %         else
% %             biomeq(var) = 1E-20;
% %         end
% %     else
% %         j = j + 1;
% %         xdot(var) = ep_detritus_import(detritus_no(j)) + flow2detritus(j) - loss(var);
% %         if (loss(var) ~= 0 && biomass(var) > 0 && flow2detritus(j) > 0)
% %             biomeq(var) = (ep_detritus_import(detritus_no(j)) + flow2detritus(j))/...
% %                 (loss(var)/biomass(var));
% %         else
% %             biomeq(var) = 1E-20;
% %         end
% %     end
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %END DERIV
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FirstTime = 1;

%%%%% this is the rate of sinks in each state variable in the derivs()
rrate = zeros(nvars);
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












