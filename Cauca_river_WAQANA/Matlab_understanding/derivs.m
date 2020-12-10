
% Compute derivatives

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DERIV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Globals

% set predator abundance
b_derivs = biomass;
%integsum = zeros(nvars,1);
for var = 1:nvars
    if (biomass(var) < 1.0e-20)
        b_derivs(var) = 1.0e-20;
    end
    if (integrate(var) >= 0)
        es_pred(var) = b_derivs(var);
    end
    integsum(var) = b_derivs(var);
end

% nutrient biomass
NutBiom = 0;
for var = 1:nvars
    NutBiom = NutBiom + biomass(var);
end

if (time == 0 && imonth==0)
    NutFree = NutTot - NutBiom;
else
% % #ifdef _ForceNutrient_
% %     if (time == 1.0) then
% %         NutFree = NutTot - NutBiom
% %     else
% %         NutFree = NutTot * NutrientForce(imonth) - NutBiom
% %     end if
% % #endif
    NutFree = NutTot - NutBiom;
end

% amount of free nutrients (NutFree) must not be less than
% the background nutrient concentration (NutMin)
if (NutFree < NutMin)
    NutFree = NutMin;
end

j = 0;
%es_PoB_base = zeros(nvars,1);
for var = 1:nvars
    if (ep_org_type(var) ~= 0)
% P/B scaled by Max. Relative Production, i.e. es_data%pob_max.
        if (time == 0 && imonth==0)
            Pmult = 1.0;
            if (ep_org_type(var) ~= 0)
                es_abs_PoB_max(var) = es_rel_PoB_max(var)* ep_PoB(var);
                es_PoB_biomass(var) = 0;
            end
        else
% % #ifdef _ForcePrimaryProd_
% %             if (ep_data(var)%org_type == 1 .and. time /= 1.0) then
% %                 j = j + 1
% %                 Pmult = PrimaryProdForce(imonth, j)
% %             end
% % #endif
% % #ifndef _ForcePrimaryProd_
            Pmult = 1.0;
% #endif
        end
        es_PoB_base(var) = 2*NutFree/(NutFree + NutFreeBase(var))*Pmult*...
            es_abs_PoB_max(var)/ (1 + biomass(var) * es_PoB_biomass(var));
   end
end

% set relative prey switching
% arena_rela_switch = zeros(vsize);
% es_pred_den = zeros(nvars,1);
for pred = 1:vsize(2)
    es_pred_den(pred) = 0;
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            es_pred_den(pred) = es_pred_den(pred) + arena_a(prey,pred)*...
                biomass(prey)^es_switch_power(pred);
            arena_rela_switch(prey,pred)=1;
        end
    end
end

% arena_base_time_switch = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            if es_switch_power(pred)>0
                arena_base_time_switch(prey,pred) = arena_a(prey,pred)*ep_biomass(prey)^...
                    es_switch_power(pred)/(es_pred_den(pred) + 1E-20)/arena_base_time_switch(prey,pred);
            end
        end
    end
end
% end set relative prey switching

% Vulnerability Calculations
% es_hdent = zeros(nvars,1);
% arena_v_denom = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            es_hdent(pred) = 0;
            arena_v_denom(prey, pred) = 0;
        end
    end
end

dwe = 0.5;
% arena_a_eff = zeros(vsize);
% arena_v_eff = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        arena_v_denom(prey, pred) = 0;
        if (vul(prey,pred) ~= -999)
            arena_a_eff(prey, pred) = arena_a_link(prey, pred) * es_Ftime(pred) *...
              arena_rela_switch(prey, pred);

            arena_v_eff(prey, pred) = arena_vul_arena(prey, pred) * es_Ftime(prey);

            arena_v_denom(prey, pred) = arena_v_denom(prey, pred) + arena_a_eff(prey, pred) *...
              es_pred(pred) / es_hden(pred);
        end
    end
end

% arena_v_biom = zeros(vsize);
for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            arena_v_biom(prey, pred) = arena_v_eff(prey, pred) * biomass(prey)/...
              (arena_vul_arena(prey, pred) + arena_v_eff(prey, pred) + arena_v_denom(prey, pred));
        end
    end
end

for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            es_hdent(pred) = es_hdent(pred) + arena_a_eff(prey, pred)*arena_v_biom(prey, pred);
        end
    end
end

for pred = 1:vsize(2)
    es_hden(pred) = (1 - dwe) * (1 + es_htime(pred)*es_hdent(pred)) + dwe*es_hden(pred);
end

for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        arena_v_denom(prey, pred) = 0;
        if (vul(prey,pred) ~= -999)
            arena_v_denom(prey, pred) = arena_v_denom(prey, pred) + arena_a_eff(prey, pred)*...
              es_pred(pred) / es_hden(pred);
        end
    end
end

for pred = 1:vsize(2)
    for prey = 1:vsize(1)
        if (vul(prey,pred) ~= -999)
            arena_v_biom(prey, pred) = arena_v_eff(prey, pred)*biomass(prey)/...
              (arena_vul_arena(prey, pred) + arena_v_eff(prey, pred) + arena_v_denom(prey, pred));
        end
    end
end
% end Vulnerability Calculations

%%%%% calculate processes; production and consumption
% primary production
% es_pp = zeros(nvars,1);
for var = 1:nvars
    if (ep_org_type(var) == 1)
        es_pp(var) = es_PoB_base(var)*biomass(var);
    else
        es_pp(var) = 0;
    end
end

% consumption
% es_qq = zeros(nvars,1);
% es_EatenBy = zeros(nvars,1);
for pred = 1:nvars
    if (pred <= vsize(2))
        consumption = 0;
        for prey = 1:vsize(1)
            if (vul(prey, pred) ~= -999)
                consumption = consumption + (arena_a_eff(prey, pred)*...
                  arena_v_biom(prey, pred) * es_pred(pred)/es_hden(pred));
            end
        end
        es_qq(pred) = consumption + es_QBoutside(pred)*biomass(pred);
        es_EatenBy(pred) = es_qq(pred);
    else
        es_qq(pred) = 0;
        es_EatenBy(pred) = es_qq(pred) / biomass(pred);
    end
end

% % unassimilated consumption
% es_unassimilated = zeros(nvars,1);
for var = 1:nvars
    if (es_qq(var) ~= 0)
        es_unassimilated(var) = es_qq(var) * ep_unass_Q(var);
    else
        es_unassimilated(var) = 0;
    end
end

% % non-predation mortalities
% es_M0 = zeros(nvars,1);
for var = 1:nvars
    if (ep_org_type(var) == 1)
        es_M0(var) = ((1 - ep_EE(var)) * ep_PoB(var)) * biomass(var);
    elseif (ep_org_type(var) == 2)
        es_M0(var) = ((1 - ep_EE(var)) * ep_PoB(var))*(1 - es_M0_pred(var) + ...
            es_M0_pred(var)* es_Ftime(var)) * biomass(var);
    else
        es_M0(var) = 0;
    end
end

% % predation mortalities
% es_M2 = zeros(nvars,1);
% es_EatenOf = zeros(nvars,1);
for prey = 1:nvars
    M2_predation = 0;
    for pred = 1:vsize(2)
        if (vul(prey, pred) ~= -999)
            M2_predation = M2_predation + (arena_a_eff(prey, pred)*...
              arena_v_biom(prey, pred)*es_pred(pred)/es_hden(pred));
        end
    end
    es_M2(prey) = M2_predation;
    es_EatenOf(prey) = M2_predation;
end


% call Calculate Detrital Flows
% flow2detritus = zeros(ndetritus,1);
for j = 1:ndetritus
    flow2detritus(j) = 0;
    for var = 1:nvars
        if (ep_org_type(var) ~= 0) 
            % calculate non-predation natural mortality 
            % for groups other than detritus
            flow2detritus(j) = flow2detritus(j) + ((1 - ep_EE(var))*...
              ep_PoB(var))*biomass(var)*ep_detfate(var,j);
        end

        if (ep_org_type(var) == 2) 
            % if group is a consumer, consider calculating excretion
            flow2detritus(j) = flow2detritus(j) + ep_unass_Q(var)*...
              es_EatenBy(var)*ep_detfate(var, j);
        end

        if (ep_discards(var) ~= 0)
            if (FirstTime == 1)
                flow2detritus(j)  = flow2detritus(j) + (ep_landings(var)+...
                ep_discards(var))*(ep_discards(var)/(ep_landings(var) + ep_discards(var)));
            else
                % then calculate flows from discarded fish
                flow2detritus(j)  = flow2detritus(j) + biomass(var)*...
                  es_fishmort(var)*(ep_discards(var)/(ep_landings(var) + ep_discards(var)));
            end
        end
     end
end

for j = 1:ndetritus
    for i = 1:ndetritus
        if (i ~= j)
            % ! calculate flows between detritus groups
            flow2detritus(j) =  flow2detritus(j) + ...
                (ep_DetPassedProp(detritus_no(i)) * ep_biomass(detritus_no(i)) * ep_detfate(i, j));
        end
    end
end

%%%% initialize the export rate of detritus out of the system
% det_export_rate = zeros(ndetritus,1);
if (FirstTime ==1)
    for j = 1:ndetritus
        det_export_rate(j) = (flow2detritus(j) - es_EatenOf(detritus_no(j))...
            + ep_detritus_import(detritus_no(j))) / biomass(detritus_no(j));
    end
end

FirstTime = 0;
% end Calculate Detrital Flows

% sum of sinks for state variables
j = 0;
for var = 1:nvars
    if (ep_org_type(var) ~= 0)
        loss(var) = es_M0(var) + es_M2(var) + es_fishmort(var)*biomass(var);
    else
        j = j + 1;
        loss(var) = det_export_rate(j) * biomass(var) + es_EatenOf(var);
    end
end

% now calculate state equations
j = 0;
for var = 1:nvars
    if (ep_org_type(var) == 2)
        xdot(var) = (ep_PoQ(var) * es_qq(var) + biomass(var)*es_PoB_base(var)) - loss(var);
        if (loss(var) > 0 && biomass(var) > 0)
            biomeq(var) = (ep_PoQ(var)*es_qq(var) + biomass(var)*...
              es_PoB_base(var))/(loss(var)/biomass(var));
        else
            biomeq(var) = 1E-20;
        end
    elseif (ep_org_type(var) == 1)
        xdot(var) = es_pp(var) - loss(var);
        if (loss(var) > 0 && biomass(var) > 0)
            biomeq(var) = es_pp(var)/(loss(var) / biomass(var));
        else
            biomeq(var) = 1E-20;
        end
    else
        j = j + 1;
        xdot(var) = ep_detritus_import(detritus_no(j)) + flow2detritus(j) - loss(var);
        if (loss(var) ~= 0 && biomass(var) > 0 && flow2detritus(j) > 0)
            biomeq(var) = (ep_detritus_import(detritus_no(j)) + flow2detritus(j))/...
                (loss(var)/biomass(var));
        else
            biomeq(var) = 1E-20;
        end
    end
end

%DERIV FUNCIONA BIEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END DERIV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
