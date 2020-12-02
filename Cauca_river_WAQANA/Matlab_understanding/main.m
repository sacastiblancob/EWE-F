%
% Pontificia Universidad Javeriana
% Maestría en Hidrosistemas
% Taller-Seminario de Investigación 2
% Taller - Modelos Acoplados
% Sergio Castiblanco
%

%%
%Hidrodinámica

%Profunidad en metros, promedio es de 1.6 m, en Munera et. al puede verse
%la variación del nivel en la Boca de la Barra. En la mañana es marea alta,
%en la tarde es marea baja, la variación promedio es de 0.3 metros,
%combinando el promedio y las variaciones se puede construir una variación
%diaria promedio del nivel, asumiendo que el total de la ciénaga cambia de
%la misma forma

%Variación de la marea boca de la barra entre las 1 y las 24
M = [0.7;0.75; 0.78; 0.8; 0.78; 0.78; 0.73; 0.68; 0.6; 0.56; 0.52; 0.5; ...
    0.53; 0.56; 0.56; 0.57; 0.57; 0.56; 0.55; 0.54; 0.56; 0.58;...
    0.61; 0.65; 0.7];
promM = mean(M);

%Nivel de la ciénaga entre las 1 y las 24 (H)
promH = 1.6;
% H = promH + (M-promM);
H = promH*ones(25,1);

%%
% Fitoplancton

%Parte de producción de fitoplancton que se va para abastecimiento de su 
%funcionamiento (respiración)
gammaf = 0.1;       

%Tasa específica de mortandad natural de fitoplancton (h^-1)
muf = 0.4;          

%Velocidad de sedimentación (Álvarez, 1987), promedio 
%Oscillatoria-Lanceaformis son las más similares al tipo de fitoplancton
%encontrado en la CGSM según Manual Fitoplancton CGSM (tipo oscillatoria)
% en (m/h)
%wg = (3600/1E12)*0.082; 
wg = 0.1;

%Tasa máxima posible de crecimiento del fitoplancton (dia^-1)
Vfmax = 10;

%Coeficiente de integral de extinción de luz con la profunidad (m^-1)
% alfa = 2.8 + 0.028*Bf;

%Flujo de la radiación solar fotosintética activa a través de la 
% superficie del mar, 160 (W m^-2) es el máximo anual segun el paper
%En Hernandez-Jimenez se muestra la variación a lo largo de un día de la 
%intensidad lumínica
%Variación de radiación solar entre las 1 y las 24
maxI = 160;
% Io = [0;0; 0; 0; 0; 0; 0.01*maxI; 0.09*maxI; 0.18*maxI; 0.3*maxI; 0.41*maxI;...
%     0.48*maxI; 0.5*maxI; 0.49*maxI; 0.45*maxI; 0.41*maxI; 0.36*maxI; 0.25*maxI;...
%     0.03*maxI;0;0;0;0;0;0];
Io = maxI*ones(25,1);
% Io = 0.5*max(Io) + Io*0.5;

%Luminosidad óptima para fotosíntesis (W m^-2) (Valor según el paper)
Iopt = 95;

%Relación Io/Iopt;
%Ro = Io./Iopt;

%R
%Rh = Ro.*exp(-alfa*H);

%Función f1(I)limitación de crecimiento de fitoplancton por luz
%f1I = (2.718./(H*alfa)).*(exp(-Rh) - exp(-Ro));

%Constantes de semi-saturación de la intensidad del proceso de utilización
%de las formas minerales del nitrógeno y fósforo por fitoplancton
CkPO4 = 0.006;
CkN = 0.025;

%%
%Bacterioplancton

%Tasa máxima posible de crecimiento de bacterias (h^-1)
Vbmax = 1.5;

%Biomasa máxima posible de bacterias (mgC m^-3)
Bbmax = 2500;

%Constante de semisaturación de crecimiento (mgC m^-3)
Bkorg = 3750;

%Constante de semisaturación del proceso de oxidación de la materia
%orgánica muerta y nitrificación por el déficit del contenido de oxígeno en
%el agua de mar
CkO2 = 1.0;

%%
%Zooplancton

%Tasa específica máxima posible de crecimiento de zooplancton (h^-1)

%Tasa específica de mortalidad del zooplancton (h^-1)
%Bajada un orden de magnnitud (orginal 0.7)
%muz = 0.7;
muz = 0.01;


%Tasa específica de excreciones metabólicas de zooplancton (h^-1)
%Bajada un orden de magnitud (original 0.1)
%gammaz = 0.1;
gammaz = 0.01;


%Coeficientes de utilización de los alimentos por zooplancton en su
%crecimiento, coeficientes de asimilación (qué proporción de lo que come se
%vuelve efectivamente biomasa de zooplancton)
om1 = 0.6;
om2 = 0.6;
om3 = 0.6;
om4 = 0.6;
omegas = [om1; om2 ;om3 ; om4];

%Coeficientes de selectividad en el consumo de distintos tipos de
%alimentación por el zooplancton, (dietas, deben sumar 1)
rof = 0.5;     %Cuanto come en proporción al total de Fitoplancton
rob = 0.2;     %Cuanto come en proporción al total de Bacterioplancton
roz = 0.1;     %Cuanto come en proporción al total de si mismo
rod = 0.2;     %Cuanto come en proporción al total de detritos
rhos = [rof ; rob ; roz ; rod];

%Constante de semisaturación por los alimentos en el proceso de crecimiento
%de zooplancton (mgC m^-3)
BK = 4150;

%Tasa máxima de forrageo específico (grazing rate) (dia ^-1)
g = 0.75;

%%
% Materia orgánica muerta

%Coeficiente economico que tiene en cuenta los gastos de energía para el
%crecimiento
Theta = 0.33;

%Aporte de materia orgánica alóctona desde las fuentes externas (mgC m^-3 h^-1)
Qextorg = 0.0;

%%
%Fósforo de fosfatos

%Coeficiente estequeométrico de traspaso de mgC a mgP para la materia
%orgánica muerta (mgP mgC^-1)
betaPC = 0.024;

%Entrada de fosfato de fuentes externas (mgP L^-1 h^-1)
QextPO4 = 0.319;

%Tasa de cambio de concentración de fosfatos por el intercambio de su masa
%entre el agua y los sedimentos del fondo (mgP L^-1 h^-1)
QHPO4 = 0.0;

%Coeficiente de traspaso de m^3 a litros L
betam3L = 0.001;

%%
%Nitrógeno Amoniacal

%Tasa específica de la primera fase de nitrificación (h^-1)
nuN1 = 0.3;

%Coeficiente estequeométrico de traspaso de mgC a mgN (mgN mgC^-1)
betaNC = 0.176;

%Tasa de intercambio del amonio por intercambio con el fondo (mgN L^-1 h^-1)
QextNH4 = 2.67;

%Parte del nitrógeno mineral consumida por el fitoplancton en forma
%amoniacal
% chi = NH4*Phi/(Phi*NH4 + (1-Phi)*CNO3)
%Phi = 0.95;
Phi = 0.95;      %Preferencia del fitoplancton por amonio respecto a nitritos y nitratos

%%
%Nitrógeno de Nitritos

%Tasa específica de la segunda fase de nitrificación (h^-1)
%nuN2 = 3.0;
nuN2 = 0.005;

%Tasa de entrada de nitritos de fuentes externas (mgN L^-1 h^-1)
QextNO2 = 0.0;

%%
%Nitrógeno de Nitratos

%Tasa de entrada de nitratos de fuentes externas (mgN L^-1 h^-1)
QextNO3 = 0.0;

%%
%Oxígeno disuelto

%Coeficiente de transpaso de mgC a mgO2 (mgO2 mgC^-1)
betaO2C = 2.67;

%Coeficiente de transpaso de de oxígeno a primera y segunda forma de
%nitrógeno (mgO2 mgN^-1)
betaO2N1 = 3.4;
betaO2N2 = 1.1;

%Coeficiente de transpaso a oxigeno por producción en el fitoplancton (mgC mgCla^-1)
betaO2cla = 18;

%Consumo de oxígeno desde una columna de agua con la base unitaria como
%resultado de los procesos bioquímicos de oxidación y nitrificación en los
%sedimentos de fonto (mg02 m^-2 h^-1)
% QbotO2 = a*O2^b;
a = 27;
b = 0.66;

%Flujo de oxígeno por intercambio con la atmósfera (mg02 m^-2 h^-1)
% QatmO2 = zetaei*nnu*nt*(O2s - O2)
zetae = 22.0;       %coeficientes de invación evasión (L m^-2 h^-1)
zetai = 11.5;
nnu = 1.0;           %Coeficiente integral del viento
nt = 1.0;           %Coeficiente térmico
O2s = 9.0;          %Concentración de saturación del oxígeno

%%
%%%% CONDICIONES INICIALES %%%%
%Valores típicos (1993-2000)
Bfo = 240;      %Fitoplancton (mgCla m^-3) 
Bbo = 1300;      %Bacterioplancton (mgC m^-3)
Bzo = 3400;      %Zooplancton (mgCla m^-3)
Borgo = 20;    %Materia Orgánica Muerta (mgCla m^-3)
PO4o = 0.09;    %Fósforo de fosfatos (mgP L^-1)
NH4o = 0.06;     %Nitrógeno amoniacal (mgN L^-1)
NO2o = 0.006;   %Nitrógeno de nitritos (mgN L^-1)
NO3o = 0.04;   %Nitrógeno de nitratos (mgN L^-1)
O2o = 6.0;      %Oxígeno disuelto (mg L^-1)

%Tiempo
to = 0;         %horas
tf = 27*100;        %horas
dt = 0.1;
time = to:dt:tf;
tsteps = length(time);

%Storage of results
RES = zeros(tsteps,9);
j=1;
RES(j,1) = Bfo;
RES(j,2) = Bbo;
RES(j,3) = Bzo;
RES(j,4) = Borgo;
RES(j,5) = PO4o;
RES(j,6) = NH4o;
RES(j,7) = NO2o;
RES(j,8) = NO3o;
RES(j,9) = O2o;

%Interpolación para H e Io
Hi = spline(0:24,H,to:dt:tf);
Ioi = spline(0:24,Io,to:dt:tf);
Ioi(abs(Ioi)<0.5) = 0.0;

%Loop temporal
for t = to+dt:dt:tf
    
%     if j==10
%         pause
%     end
    
    ti = mod(t,24);
    %%%% RHS's (Términos de mano derecha) %%%%
    
    %Relaciones Predadores-Presa
    Bs = [Bfo ; Bbo; Bzo ; Borgo];
    sumrhoB = sum(Bs.*rhos);
    
    Ps = (Bs.*rhos)./sumrhoB;
    sumPsB = sum(Bs.*Ps);
    
    Gs = g*(Bs.*Ps)./(BK + sumPsB);
    
    %Function f1(l)
    alfa = 2.8 + 0.028*Bfo;
    %alfa = 20;
    Ioo = spline(0:24,Io,ti);
    if abs(Ioo) < 1.0
        Ioo = 0.0;
    end
    Ro = Ioo/Iopt;
    Rh = Ro*exp(-alfa*spline(0:24,H,ti));
    f1 = (2.718/(alfa*spline(0:24,H,ti)))*(exp(-Rh) - exp(-Ro));
    
    %Function f2(CN,CPO4)
    N = NH4o + NO2o + NO3o;
    f2 = min((N^2/(CkN^2 + N^2)),(PO4o^2/(CkPO4^2 + PO4o^2)));
    
    %Function sigma(Iz,CPO4,CN)
    sigmaf = Vfmax*f1*f2; 
    
    %RHS de Fitoplancton
    Bfrhs = ((1-gammaf)*sigmaf - muf)*Bfo - wg/spline(0:24,H,ti) - Gs(1)*Bzo;
    
    %Mu_{b}
    mub = Vbmax*(Bbo/Bbmax);
    
    %Eps_{ing}
    epsing = O2o/(O2o + CkO2);
    
    %RHS de Bacterioplancton
    Bbrhs = (Vbmax*(Borgo/(Bkorg + Borgo))*epsing - mub)*Bbo - Gs(2)*Bzo;
    
    %RHS de Zooplancton
    Bzrhs = (rhos(1)*Gs(1) + rhos(2)*Gs(2) + (rhos(3)-1)*Gs(3) + ...
        rhos(4)*Gs(4) - gammaz - muz)*Bzo;
    
    %RHS de Materia orgánica muerta
    Borgrhs = muf*Bfo + ((1-rhos(1))*Gs(1) + (1-rhos(2))*Gs(2) + ...
        (1-rhos(3))*Gs(3) - rhos(4)*Gs(4) + muz)*Bzo + ...
        (mub - (Vbmax/Theta)*(Borgo/(Bkorg + Borgo))*epsing)*Bbo + Qextorg;
    
    %RHS de Fósforo de Fosfatos
    PO4rhs = ((1/Theta)-1)*Vbmax*(Borgo/(Bkorg + Borgo))*epsing*Bbo + ...
        gammaz*Bzo - (1-gammaf)*sigmaf*Bfo;
    PO4rhs = PO4rhs*betaPC*betam3L + QextPO4 + QHPO4;
    
    %Parte del nitrógeno mineral consumida por el fitoplancton en forma
    %amoniacal
    chi = NH4o*Phi/(Phi*NH4o + (1-Phi)*NO3o);
    
    %RHS de Nitrógeno Amoniacal
    NH4rhs = ((1/Theta)-1)*Vbmax*(Borgo/(Bkorg + Borgo))*epsing*Bbo + ...
        gammaz*Bzo - (chi-gammaf)*sigmaf*Bfo;
    NH4rhs = NH4rhs*betaNC*betam3L - nuN1*epsing*NH4o + QextNH4;
    
    %lambda_{NO2}
    lamNO2 = NO2o/(NO2o + NO3o);
    
    %RHS de Nitrógeno de Nitritos
    NO2rhs = nuN1*epsing*NH4o - nuN2*epsing*NO2o - ...
        (1-chi)*lamNO2*sigmaf*Bfo*betaNC*betam3L + QextNO2;
    
    %lambda_{NO3}
    lamNO3 = NO3o/(NO2o + NO3o);
    
    %RHS de Nitrógeno de Nitratos
    NO3rhs = nuN2*epsing*NO2o - ...
        (1-chi)*lamNO3*sigmaf*Bfo*betaNC*betam3L + QextNO3;
    nuN2*epsing
    
    %Absorción de oxígeno en el fondo
    QbotO2 = a*(O2o^b);
    
    %Intercambio de oxígeno con la atmósfera
    %QatmO2 = zetae*zetai*nnu*nt*(O2s - O2o);
    QatmO2 = (zetai-zetae)*nnu*nt*(O2s - O2o);
    
    %RHS de Oxígeno Disuelto
    O2rhs = (sigmaf - gammaf)*Bfo - ...
        ((1/Theta)-1)*Vbmax*(Borgo/(Bkorg + Borgo))*epsing*Bbo -...
        gammaz*Bzo;
    O2rhs = O2rhs*betaO2C*betam3L;
    O2rhs = O2rhs - (nuN1*betaO2N1*NH4o + nuN2*betaO2N2*NO2o)*epsing;
    O2rhs = O2rhs - (QbotO2 + QatmO2)/spline(0:24,H,ti);
    
    %%%% Advancing in time
    Bfo = Bfo + Bfrhs*dt;
    Bbo = Bbo + Bbrhs*dt;
    Bzo = Bzo + Bzrhs*dt;
    Borgo = Borgo + Borgrhs*dt;
    PO4o = PO4o + PO4rhs*dt;
    NH4o = NH4o + NH4rhs*dt;
    NO2o = NO2o + NO2rhs*dt;
    NO3o = NO3o + NO3rhs*dt;
    O2o = O2o + O2rhs*dt;
    
    %Storage of results
    j=j+1;
    RES(j,1) = Bfo;
    RES(j,2) = Bbo;
    RES(j,3) = Bzo;
    RES(j,4) = Borgo;
    RES(j,5) = PO4o;
    RES(j,6) = NH4o;
    RES(j,7) = NO2o;
    RES(j,8) = NO3o;
    RES(j,9) = O2o;
 
end

figure(1)

plot(time,RES(:,1:4))
xlabel('Time (h)')
ylabel('Plancton (mgC m^{-3})')
legend('Fitoplancton','Bacterioplancton','Zooplancton','Mat. Orgánica','Location','Northeast')

figure(2)
plot(time,RES(:,5:9))
xlabel('Time (h)')
ylabel('Concentración (mg/L')
legend('PO4','NH4','NO2','NO3','O2','Location','Northeast')













