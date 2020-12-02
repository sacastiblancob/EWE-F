#' @name
#' GSAtool
#'
#' @title
#' Global Sensitivity Analysis tool
#'
#' @description
#' This function performs the global sensitivity analysis starting from the gross results of the model.
#'
#'
#' @param parameters_set matrix of dimensions n x pp, where n is the
#' number of runs and pp is the number of parameters.
#'
#' @param out_set matrix of dimensions n x t, where n is the number of
#' runs and t is the number of temporary steps.
#'
#' @param pp_names a strings vector with the names of the parameters of the model
#'
#' @param steps number of divisions of the parametric range.
#'
#' @param save T to save the results in .csv files, by default save=F.
#'
#' @param dir a directory  to save the results
#'
#' @return
#' a list containing two outputs: SOBOL and AMA indices.
#'
#' @import stats e1071 utils
#'
#' @export
#'
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
#' Sobol, I. M. (2001). Global sensitivity indices for nonlinear mathematical models
#' and their Monte Carlo estimates. Mathematics and Computers in Simulation, 55(1–3),
#' 271–280. https://doi.org/10.1016/S0378-4754(00)00270-6
#'
#' @examples
#' data("parameters_set", "out_set", "pp_names")
#' \donttest{
#'
#' GSA_results <- GSAtool(parameters_set, out_set, pp_names, steps = 15, save=FALSE)
#' }

GSAtool <- function(parameters_set, out_set, pp_names, steps = 100, save=FALSE, dir=NULL){


  data_Bstat <- Bstat(out_set)

  CM <- Cond_Moments(parameters_set, out_set , pp_names, steps = steps)

  SOBOL_indices <- SOBOL(data_var = data_Bstat[,3], CM_mean = CM$CM_mean, CM_var = CM$CM_var, pp_names = pp_names)

  AMA_indices <- AMA(data_Bstat , CM, pp_names, steps = steps)

  if (save==TRUE){
    (save_results(SOBOL = SOBOL_indices[[1]], amae = AMA_indices$AMAE, amav = AMA_indices$AMAV,
                  amar = AMA_indices$AMAR, amak = AMA_indices$AMAK, dir=dir))
  }

  GSA <- list(SOBOL_indices, AMA_indices)

  return(GSA)
}
