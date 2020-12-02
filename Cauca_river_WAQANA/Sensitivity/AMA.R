#' @name
#' AMA
#'
#' @title
#' AMA indices
#'
#' @description
#' This function calculates the AMA indices: AMAE, AMAV, AMAV
#' and AMAK.
#'
#' @param data_Bstat a data frame of dimensions t x 6, here t is the number of
#' temporary steps and each column corresponds to a statistical measure: mean,
#' variance, skewness, kurtosis and excess kurtosis.
#'
#' @param CM A list of arrays, each array corresponds to the conditional
#' moments calculated with the mean, variance, skewness, kurtosis. Each
#' array has dimensions of steps, t, p.
#'
#' @param pp_names vector that contains the names of the parameters (pp)
#'
#' @param steps number of divisions of the parametric range
#'
#' @return
#' A list of four matrices, which corresponds to AMAE, AMAV, AMAR and AMAK indices.
#' Each matrix has dimensions of t x pp.
#'
#' @export
#'
#' @author
#' Camila Garcia-Echeverri <cagarciae@unal.edu.co> \cr
#' Maria Cristina Areas-Bautista <mcarenasb@unal.edu.co> \cr
#'
#'
#' Hydrodynamics of the natural media research group - HYDS
#' National University of Colombia -  Bogota
#'
#' @references
#' Dell’Oca, A., Riva, M., & Guadagnini, A. (2017). Moment-based metrics for global sensitivity
#' analysis of hydrological systems. Hydrology and Earth System Sciences, 21(12), 6219–6234.
#'  https://doi.org/10.5194/hess-21-6219-2017
#'
#' @examples
#' data("data_Bstat", "CM", "pp_names")
#' AMA_indices <- AMA(data_Bstat, CM, pp_names, steps= 15)
#'
#'

AMA <- function(data_Bstat, CM, pp_names, steps = 100){

  t <- dim(data_Bstat)[1]
  pp <- length(pp_names)

  me <- array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
  var <- array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
  skw <- array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
  kurt <- array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))

  amae <- matrix(nrow = t, ncol=pp, dimnames = list(NULL,pp_names))
  amav <- matrix(nrow = t, ncol=pp, dimnames = list(NULL,pp_names))
  amar <- matrix(nrow = t, ncol=pp, dimnames = list(NULL,pp_names))
  amak <- matrix(nrow = t, ncol=pp, dimnames = list(NULL,pp_names))


  for(f in 1:t){
    for (s in 1:pp){

      me[,f,s] <- (abs(data_Bstat[2][f,1]-CM$CM_mean[,f,s]))
      var[,f,s] <- (abs(data_Bstat[3][f,1]-CM$CM_var[,f,s]))
      skw[,f,s] <- (abs(data_Bstat[4][f,1]-CM$CM_skw[,f,s]))
      kurt[,f,s] <- (abs(data_Bstat[5][f,1]-CM$CM_kurt[,f,s]))

      amae[f,s] <- (1/data_Bstat[f,2]*mean(me[,f,s],na.rm=TRUE))
      amav[f,s] <- (1/data_Bstat[f,3]*mean(var[,f,s],na.rm=TRUE))
      amar[f,s] <- (1/data_Bstat[f,4]*mean(skw[,f,s],na.rm=TRUE))
      amak[f,s] <- (1/data_Bstat[f,5]*mean(kurt[,f,s],na.rm=TRUE))
    }
  }

  AMA_indices <- list("AMAE"=amae, "AMAV"=amav, "AMAR"=amar, "AMAK"=amak)

  return(AMA_indices)
}
