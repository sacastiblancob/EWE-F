
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First subroutine after read entries
%For producers, calculateMaximumPoBrelatedValues
aB = 16;
arPoB_max = 2;
aPoB = 2/3;
aaPoB_max = arPoB_max*aPoB;
aPoB_biomass = (aaPoB_max/aPoB - 1)/aB;
aPoB_biomass2 = (arPoB_max - 1)/aB;

% then es_PoB_biomass = (es_rel_PoB_max - 1)/ep_biomass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2, subroutine calculateNutrientConcentration
NutBiom = sum(ep_biomass);  %toda la biomasa sumada

NutBaseFreeProp = 0.99;      % given in the configuration, filenames.mnl

NutTot = NutBiom/(1 - NutBaseFreeProp); %This are the disponible nutrients

NutFree = NutTot  - NutBiom; %disponible nutrients minus all biomass

% and then asign it to NutFreeBase for all detritus groups

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2, subroutine removeImportFromDiet
% this is only for consumers
% es_QB_outside tiene la parte del consumo sobre biomasa que los
% consumidores obtienen de afuera del sistema

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3, en el código hace algunas cosas on the fly en ecosim.F

% % para consumidores
% es_Ftime = 1
% es_hden = es_QB_maxoQB_0/(es_QB_maxoQB_0 + 1);
% es_QB_maxoQB_0 lo introduce uno mismo en la 8va columna de scenario
% descripcion "QBmax/QBo (for handling time) [>1]", está en 1000 en
% Tamba_Bay examplo para todos los grupos, en cuyo caso es_hden da muy
% cercano a 1.0

% % para todos
% es_CB_base = ep_EatenBy / biomass;

% % si es_CB_base = 0, que sería para productores y detritus entonces
% es_CB_base = 1
% es_Ftime_max = 1

% % para todos
% es_CB_last = es_CB_base
% explicacion codigo "es_CB_base is base consumption biomass ratio
% calculated from initial conditions"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocatea B con nvars, y si hay substanzas allocatea es_ms_data con el
% número de grupos con multistanza, y entra a la siguiente rutina
%
% initialiseSplitGroups(B)
%
% entra a un ciclo para cada grupo con multistanza stanza = 1,nstanzas, y luego a otro para
% cada age = 0,ms_age_infinity y hace lo siguiente:
%
% es_ms_NageS = ms_SplitNo para cada stanza y para cada age
%    ms_SplitNo es "number of survivors at age (monthly)"
%
% es_ms_WageS = ms_Wage para cada stanza y para cada age
%    ms_Wage es "relative body weight at age "a" (monthly)"
%
% if es_ms_Wage(age) es mayor que Wmat_Winf para esa stanza entonces
%    es_ms_EggSplit = es_ms_WageS - ms_Wmat_Winf
%
%    Be = Be + es_ms_NageS*es_ms_EggsSplit
% end if
% end ciclo de age
%
% Be termina siendo como el acumulado del número de supervivientes que
% tienen un peso relativo mayor a la relación Wmat_Winf, y eso se lo
% asignan al ms_BaseEggsStanza y al ms_EggsStanza
%
% es_ms_BaseEggsStanza = Be
% es_ms_EggsStanza = Be
% es_ms_RscaleSplit = 1;

% aquí llaman a setSplitPred(B) dentro de initialiseSplitGroups(B)
%
%   WARNING: Esta rutina setSplitPred la llaman cada vez que se completa un mes en la simulación,
%   por eso el parámetro de entrada es B, porque cada mes cambia la
%   proporción de biomasa relativa por substanza
%
%   aqui hacen un ciclo para cada grupo con multistanza y para cada
%     substanza
%
%       inicializan Bt, Pt y Nt en 1.0e-30 (dummy variables)
%       calculan el age_last para cada substanza
%
%       hacen un ciclo de age = age_start,age_last
%
%           Bt = Bt + es_ms_NageS*es_ms_WageS
%           Pt = Pt + es_ms_NageS*ms_WWa
%           Nt = Nt + es_ms_NageS
%       fin del ciclo para age
%       B(grupono_substanza) = Bt
%       es_pred(groupno_substanza) = Pt
%
%       Bt tiene el acumulado de NageS*WageS para cada substanza, que sería
%       "number of survivors at age (monthly)" multiplicado por 
%       "relative body weight at age "a" (monthly)", luego lo guarda en B,
%       esto implicaria que guarda la biomasa relativa de cada substanza
%
%       Pt tiene el acumulado de NageS*WWa para cada substanza, que sería
%       "number of survivors at age (monthly)" multiplicado por
%       "relative consumption (Q) depending on Wage", esto sería que guarda
%       el consumo relativo de cada substanza en es_pred
%    fin de los ciclos de nstanza y substanza
% fin de setSplitPred
%
% inicia un ciclo para cada grupo con substanzas y otro para cada substanza
%
%  Agem es la máxima edad para cada substanza, en meses, aunque para la
%  stanza más vieja le restan un mes, debido a que usan age+1 en algún
%  punto adelante
%
%     inicia un cilo para la edad age=age_star, Agem
%
%       es_ms_SplitAlpha(age) = (ms_Wage(age+1) -
%       vbM*ms_Wage(age))*es_pred/ep_EatenBy      --> para cada substanza
%
%     fin del ciclo para age
%     es_ms_SplitAlpha(age_infinity) = es_ms_SplitAlpha(age_infinity-1)
%
%     es_ms_SplitAlpha tiene ("relative body weight at age "a+1" (monthly)"
%     - "von Bertalanffy metabolic parameter" * "relative body weight at age "a" (monthly)") *
%     ("el consumo relativo de cada substanza (es_pred)"/"the consumption of substanza"
%
% fin de initialiseSplitGroups(B), sale B con la biomasa relativa de cada
% substanza
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% de vuelta en ecosim.F90
%
% acá le asignan integrate=-i si el grupo es multistanza, o integrate=i si
% el grupo no es multistanza
%
% !!!!!! setpred() !!!!!!!!1
%
% allocatean b_pred(nvars), al que le asignan el mismo valor de la biomasa
% de ecopath para cada grupo o 1.0e-20 si la biomasa es menor a 1.0e-20
%
% y le asignan a es_pred la biomasa para los grupos que NO son multistanza,
% como se vió más arriba este valor de es_pred se calcula para los grupos
% que si son multistanza
%
% luego le ponen es_risk_time(i)=0 para aquellos grupos donde el
% es_risk_time es == -999, este es un parametro de entrada, la 6ta columna
% de scenario, uno debe poner -999 para productores y detritus, aquí se
% cambia el -999 por 0.0
%
% !!!!! calculate risk time for consumers !!!!!
%
% ciclo for i=1,nvars
%
% si el grupo es consumidor entonces
%       es_CB_base = ep_EatenBy / es_pred    % consumo / consumo relativo del grupo o algo así 
%       es_CB_last = es_CB_base
%       es_Q_main = (1 - es_risk_time)*es_CB_base  %consumo en tiempo de no riesgo
%       es_Q_risk = es_risk_time * es_CB_base * (ep_EatenOf/ep_biomass + ((1 - ep_EE) * ep_PoB) + 0.00000000001)
%
% end if
% end ciclo para nvars
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% call calculateLotkaVolterraEffectiveSearchRates(vrows,vcols)
% solo se calcula al inicio
%
%  Primero calculan la variable es_htime para cada grupo consumidor, descripción de
%  esta variable "Handling time for predators"
%
%  ciclo para nvars
%    si es consumidor
%      es_htime = es_pred/(es_QB_maxoQB_0 * ep_biomass * ep_QoB)
%      
%      es_htime es consumo relativo (es_pred) / ("maximum relative
%      consumption" * "biomass" * "consumo/biomasa")
%
%    si es productor o detritus
%      es_htime = 0
%    end if
%  end ciclo
%
%  un cilo para cada columna y para cada fila de la matriz de vulnerabilidad
%  do pred=1,vcols
%    do prey=1,vrows
%      
%      Dzero = es_QB_maxoQB_0(pred) / (es_QB_maxoQB_0(pred) - 1)
%      "maximum relative consumption"/("maximum relative consumption" - 1)
%      la parte de consumo que hace falta para llegar al máximo relativo
%
%      si es_vul ~= -999
%        consumption = ep_biomass(pred) * ep_QoB(pred) * ep_diet;
%        consumo de ese depredador de esa presa
%
%        arena_vulrate = (es_vul * consumption)/ep_biomass(prey)
%        relación Q(pred sobre pred)/B(pred) para cada relación predador presa, afectada por la
%        vulnerabilidad
%
%        ! Denominator for calculating search rate of predator
%        Denv = ep_biomass(prey) * es_pred(pred) * arena_vulrate - consumption * es_pred(pred)
%        
%        Despejando Denv
%        "
%        Denv = ep_biomass(prey) * es_pred(pred) * arena_vulrate -
%        ep_biomass(pred) * ep_QoB(pred) * ep_diet * es_pred(pred)
%
%        Denv = ep_biomass(prey) * es_pred(pred) * es_vul * ep_biomass(pred) * ep_QoB(pred) * ep_diet / ep_biomass(prey) -
%        ep_biomass(pred) * ep_QoB(pred) * ep_diet * es_pred(pred)
%
%        Denv = ep_biomass(pred) * ep_QoB(pred) * ep_diet * es_pred(pred) * (es_vul - 1)
%
%        nótese que si es_vul es 1.0, denv da cero
%        "
%
%        si Denv < 1.0e-20
%          Denv = 1.0e-20
%        end if
%
%        arena_a = Dzero * 2 * Consumption * arena_vulrate / Denv
%        para cada relación predador presa
%        "
%        arena_a = (es_QB_maxoQB_0(pred) / (es_QB_maxoQB_0(pred) - 1)) * 2 *
%        ep_biomass(pred) * ep_QoB(pred) * ep_diet * (es_vul *
%        ep_biomass(pred) * ep_QoB(pred) * ep_diet)/ep_biomass(prey) / 
%        ep_biomass(pred) * ep_QoB(pred) * ep_diet * es_pred(pred) * (es_vul - 1)
%
%        arena_a = (es_QB_maxoQB_0(pred) / (es_QB_maxoQB_0(pred) - 1)) * 2 *
%        ep_biomass(pred) * ep_QoB(pred) * ep_diet * es_vul/ep_biomass(prey)* 
%        es_pred(pred) * (es_vul - 1)
%
%      si es_vul = -999
%        arena_a = 0     !para cada relación predador presa
%      end si
%    end do
%  end do
%
% end calculateLotkaVolterraEffectiveSearchRates(vrows,vcols)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% de vuelta en Ecosim.F90
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call initialiseRelativeSwitchingParameters (vrows, vcols)
% 
% ciclo para predadores do pred=1,vcols
%   es_pred_den(pred) = 0
%   ciclo para presas do prey=1,vrows
%     if es_vul ~= -999
%       es_pred_den(pred) = es_pred_den(pred) + arena_a * ep_biomass(prey)^es_switch_power(pred)
%       "
%       es_pred_den = es_pred_den(pred) + "Search rate of predator j for prey i" * "prey
%       biomass" ^ "prey switching power parameter"
%       Descripción: Denominator for calculating prey switching parameters
%       arena_a es la tasa de busqueda del predador sobre la presa para
%       cada relación predador-presa
%       "
%     end
%   end
% end
%
% ciclo para predadores do pred=1,vcols
%   ciclo para presas do prey=1,vrows
%     if es_vul ~= -999
%       arena_base_time_switch = (arena_a * ep_biomass(prey)^es_switch_power(pred))/
%         (es_pred_den(pred) + 1.0e-20)
%     end
%   end
% end
% "Si el predador tiene una sola presa, el arena_base_time_switch da 1 para
% esa relación específica predador presa; de otro modo los valores de esta
% variable deben sumar 1.0 para cada columna (predador)el es_switch_power no parece
% cambiar los resultados"
%
% end initialiseRelativeSwitchingParameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% de vuelta en Ecosim.F90
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% call setArenaVulnerabilityandSearchRates(vrows,vcols)
% solo se ejecuta una vez
%
% inicializa arena_Q_arena en cero
%
%  un cilo para cada columna y para cada fila de la matriz de vulnerabilidad
%  do pred=1,vcols
%    do prey=1,vrows
%      si es_vul ~= -999
%
%        ! total consumptions of predators in arena
%        arena_Q_link = ep_biomass(pred) * ep_QoB(pred) * ep_diet
%
%        arena_Q_arena = arena_Q_arena + arena_Q_link
%
%        ! setting of initial vulnerable biomasses
%        arena_vul_arena = (es_vul + 0.00000000001) * arena_Q_arena / ep_biomass(prey)
%        "consumo del predador relativo a la biomasa total de la presa en
%        esa arena, afectado por la vulnerabilidad"
%
%        si arena_vul_arena es cero, arena_vul_arena = 0.0
%
%        arena_vul_biom = (es_vul + 0.00000000001 - 1) * arena_Q_arena / (2 * arena_vul_arena)
%        "
%        arena_vul_biom = (es_vul - 1) * ep_biomass(pred) * ep_QoB(pred) *
%        ep_diet / (2 * es_vul * ep_biomass(pred) * ep_QoB(pred) * ep_diet / ep_biomass(prey))
%
%        arena_vul_biom = (es_vul - 1) * ep_biomass(prey)/(2 * es_vul)
%        "
%
%        si arena_vul_biom es 0, entonces arena_vul_biom == 1,
%
%        ! setting of predator search rates
%        if arena_vul_biom > 0
%
%          Dzero = es_QB_maxoQB_0(pred) / (es_QB_maxoQB_0(pred) - 1)
%
%          arena_a_link = Dzero * arena_Q_link / (arena_vul_biom * es_pred(pred))
%          "
%          arena_a_link = Dzero * ep_biomass(pred) * ep_QoB(pred) * ep_diet
%          / (((es_vul - 1) * ep_biomass(prey)/(2 * es_vul)) * es_pred(pred))
%          "
%        else
%          arena_a_link = 0
%        end
%      end if
%    end do
%  end do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% de vuelta en Ecosim.F90
%
% allocate(xdot(nvars)); allocate(biomeq(nvars));
% allocate(loss(nvars)); allocate(b_out(nvars));
%
% ciclo para nvars
%   es_fishmort = (ep_landings + ep_discards)/ep_biomass
% end
%
% FirstTime = .TRUE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% call derivs (0., ep_data(:)%biomass, xdot, biomeq, loss, integrate)
% comentarios_derivs.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















