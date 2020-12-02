%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Waqtel modules
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS FOR SENSITIVITY ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For EUTRO module
%
% Parameters are in the same order as in Waqtel Reference Guide
%
% params = [U,H,T,alfa1,kpe,beta,BEN,M1,M2,k520,k120,KN,KP,n,f,IK,dtn,dtp,fn,fp,k620...
%           k320,RP,WNOR,WLOR,WPOR,Io,Cmax]
params = [0.8,4.0,25,1.0,2.8,0.028,1.5,0.01,0.001,0.1,0.1,0.025,0.005,3.4,0.05,100,...
    0.5,0.5,0.018,0.018,0.01,0.01,0.01,0,0,0,160,2];
%
% Initial conditions
% init = [PHYo,PO4o,PORo,NO3o,NORo,NH4o,Lo,O2o]
init = [240,0.09,0.09,0.04,0.006,0.06,0.02,6];

% Calling eutro
%[RES1,RES2,RES3,RES4,RES5,RES6,RES7,RES8] = eutro(params,init,86400*60,60);

N=36000;

% inits
PHYo = 0 + (1000-0)*rand(N,1);
PO4o = 0 + (5-0)*rand(N,1);
PORo = 0 + (5-0)*rand(N,1);
NO3o = 0 + (5-0)*rand(N,1);
NORo = 0 + (5-0)*rand(N,1);
NH4o = 0 + (5-0)*rand(N,1);
Lo = 0 + (0.1-0)*rand(N,1);
O2o = 0 + (9-0)*rand(N,1);
allinits = [PHYo,PO4o,PORo,NO3o,NORo,NH4o,Lo,O2o];

%params
U = 0 + (2-0)*rand(N,1);
H = 0 + (10-0)*rand(N,1);
T = 18 + (30-18)*rand(N,1);
alfa1 = 0 + (1-0)*rand(N,1);
kpe = 0 + (6-0)*rand(N,1);
beta = 0 + (0.05-0)*rand(N,1);
BEN = 0 + (10-0)*rand(N,1);
M1 = 0 + (0.2-0)*rand(N,1);
M2 = 0 + (0.05-0)*rand(N,1);
k520 = 0 + (1-0)*rand(N,1);
k120 = 0 + (1-0)*rand(N,1);
KN = 0 + (0.01-0)*rand(N,1);
KP = 0 + (0.1-0)*rand(N,1);
n = 0 + (10-0)*rand(N,1);
f = 0 + (0.5-0)*rand(N,1);
IK = 60 + (140-60)*rand(N,1);
dtn = 0 + (1-0)*rand(N,1);
dtp = 0 + (1-0)*rand(N,1);
fn = 0 + (0.05-0)*rand(N,1);
fp = 0 + (0.05-0)*rand(N,1);
k620 = 0 + (1-0)*rand(N,1);
k320 = 0 + (1-0)*rand(N,1);
RP = 0 + (0.1-0)*rand(N,1);
WNOR = 0 + ((0.5/86400)-0)*rand(N,1);
WLOR = 0 + ((0.5/86400)-0)*rand(N,1);
WPOR = 0 + ((0.5/86400)-0)*rand(N,1);
Io = 60 + (180-0)*rand(N,1);
Cmax = 0 + (4-0)*rand(N,1);
allparams = [U,H,T,alfa1,kpe,beta,BEN,M1,M2,k520,k120,KN,KP,n,f,IK,dtn...
    ,dtp,fn,fp,k620,k320,RP,WNOR,WLOR,WPOR,Io,Cmax];

%time
tf = 86400*90;
dt = 60;
TS = 2;

RESPHY = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESPO4 = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESPOR = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESNO3 = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESNOR = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESNH4 = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESL = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));
RESO2 = zeros(N,TS + length(allparams(1,:)) + length(allinits(1,:)));


parfor i=1:N
    i
    %CALL EUTRO
    [RES1,RES2,RES3,RES4,RES5,RES6,RES7,RES8] = eutro(allparams(i,:),allinits(i,:),tf,dt);
    
    %STORING
    RESPHY(i,:) = RES1;
    RESPO4(i,:) = RES2;
    RESPOR(i,:) = RES3;
    RESNO3(i,:) = RES4;
    RESNOR(i,:) = RES5;
    RESNH4(i,:) = RES6;
    RESL(i,:) = RES7;
    RESO2(i,:) = RES8;
end

% save('RES.mat','RESPHY','RESPO4','RESPOR','RESNO3','RESNOR','RESNH4','RESL','RESO2')

% load('RES.mat')
% for i=1:length(params)+length(init)
%     RESPHY(:,i) = normalize(RESPHY(:,i));
%     RESPO4(:,i) = normalize(RESPO4(:,i));
%     RESPOR(:,i) = normalize(RESPOR(:,i));
%     RESNO3(:,i) = normalize(RESNO3(:,i));
%     RESNOR(:,i) = normalize(RESNOR(:,i));
%     RESNH4(:,i) = normalize(RESNH4(:,i));
%     RESL(:,i) = normalize(RESL(:,i));
%     RESO2(:,i) = normalize(RESO2(:,i));
% end
% save('RESNORM.mat','RESPHY','RESPO4','RESPOR','RESNO3','RESNOR','RESNH4','RESL','RESO2')





