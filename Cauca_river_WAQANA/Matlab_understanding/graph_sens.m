%
% PONTIFICIA UNIVERSIDAD JAVERIANA
% EPM-PUJ
% Sergio Castiblanco
% Understanding Waqtel modules
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPH SENSITIVITY ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = "/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES";

S1 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL1.csv",1,1);
S2 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL2.csv",1,1);
S3 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL3.csv",1,1);
S4 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL4.csv",1,1);
S5 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL5.csv",1,1);
S6 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL6.csv",1,1);
S7 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL7.csv",1,1);
S8 = csvread("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/RES/SOBOL8.csv",1,1);

S1n = S1(1,:);
S2n = S2(1,:);
S3n = S3(1,:);
S4n = S4(1,:);
S5n = S5(1,:);
S6n = S6(1,:);
S7n = S7(1,:);
S8n = S8(1,:);

S1no = 100*S1n/sum(S1n);
S2no = 100*S2n/sum(S2n);
S3no = 100*S3n/sum(S3n);
S4no = 100*S4n/sum(S4n);
S5no = 100*S5n/sum(S5n);
S6no = 100*S6n/sum(S6n);
S7no = 100*S7n/sum(S7n);
S8no = 100*S8n/sum(S8n);

pp_names = ["PHY","PO4","POR","NO3","NOR","NH4","L","O2","U","H","T",...
"alfa1","kpe","beta","BEN","M1","M2","k520","k120","KN","KP","n","f",...
"IK","dtn","dtp","fn","fp","k620","k320","RP","WNOR","WLOR","WPOR",...
"Io","Cmax"];
par_cat = categorical(pp_names);

% y = [S1no;S2no;S3no;S4no;S5no;S6no;S7no;S8no];
% 
% bar(y,'stacked')
% legend(pp_names)

close all

figure(1)
bar(par_cat,S1no)
title('Sensitivity Analysis, PHY (SOBOL)')

figure(2)
bar(par_cat,S2no)
title('Sensitivity Analysis, PO4 (SOBOL)')

figure(3)
bar(par_cat,S3no)
title('Sensitivity Analysis, POR (SOBOL)')

figure(4)
bar(par_cat,S4no)
title('Sensitivity Analysis, NO3 (SOBOL)')

figure(5)
bar(par_cat,S5no)
title('Sensitivity Analysis, NOR (SOBOL)')

figure(6)
bar(par_cat,S6no)
title('Sensitivity Analysis, NH4 (SOBOL)')

figure(7)
bar(par_cat,S7no)
title('Sensitivity Analysis, OL (SOBOL)')

figure(8)
bar(par_cat,S8no)
title('Sensitivity Analysis, O2 (SOBOL)')











