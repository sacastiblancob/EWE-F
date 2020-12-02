########################################################################
## SENSITIVITY ANALYSIS FOR WAQTEL OUTPUTS USING SOBOL, AMAE, ETC.
## Using GSA.UN (PhD Arenas-Bautista M.C.)
########################################################################
## 
##  PONTIFICIA UNIVERSIDAD JAVERIANA
##  EPM-PUJ
##  Sergio Castiblanco
##  Understanding Waqtel modules
##  Sensitivity analysis
##
########################################################################
##
## Libraries
##
library(GSA.UN)
library(stats)
library(e1071)
library(R.matlab)
homedir <- "/media/aldair/FILES/Telemac/waqtel/Cauca_river/Sensitivity/Scripts/"
##
## Loading results from matlab (RES.mat)
##
#RES <- readMat("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Matlab_understanding/RES.mat")
RES <- readMat("/media/aldair/FILES/Telemac/waqtel/Cauca_river/Matlab_understanding/RESNORM.mat")
#PHY <- RES$RESPHY
#PO4 <- RES$RESPO4
#POR <- RES$RESPOR
#NO3 <- RES$RESNO3
#NOR <- RES$RESNOR
#NH4 <- RES$RESNH4
#L <- RES$RESL
#O2 <- RES$RESO2
#
# Number of variables
#
N <- 8
#
# Number of parameters
#
P <- 28
#
# STEPS
#
steps <- 2
#
# Parameters names
#
pp_names = c("PHY","PO4","POR","NO3","NOR","NH4","L","O2","U","H","T",
"alfa1","kpe","beta","BEN","M1","M2","k520","k120","KN","KP","n","f",
"IK","dtn","dtp","fn","fp","k620","k320","RP","WNOR","WLOR","WPOR",
"Io","Cmax")
#
# Computing Loop
#
dirout <- paste(homedir,'RES/',sep='')
for(i in 1:N){
    #Result
    RTRAC <- RES[[i]]
    RTRAC[is.na(RTRAC) == TRUE] <- 0.0
    #Parameters
    parameters_set <- RTRAC[,1:(N+P)]
    #Out_set
    out_set <- RTRAC[,(N+P+1):(N+P+steps)]
    #Statistical moments
    GSAtool(
        parameters_set,
        out_set,
        pp_names,
        steps = steps,
        save = TRUE,
        dir = dirout
        )
    #outdir
    file.rename(paste(dirout,'AMAE.csv',sep=''),paste('AMAE',as.character(i),'.csv',sep=''))
    file.rename(paste(dirout,'AMAK.csv',sep=''),paste('AMAK',as.character(i),'.csv',sep=''))
    file.rename(paste(dirout,'AMAR.csv',sep=''),paste('AMAR',as.character(i),'.csv',sep=''))
    file.rename(paste(dirout,'AMAV.csv',sep=''),paste('AMAV',as.character(i),'.csv',sep=''))
    file.rename(paste(dirout,'SOBOL.csv',sep=''),paste('SOBOL',as.character(i),'.csv',sep=''))
}









