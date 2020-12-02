#' @name
#' save_results
#'
#' @title
#' Save GSA results
#'
#' @description
#' This function helps to save the results in .csv format
#'
#' @param SOBOL SOBOL index
#'
#' @param SOBOL_total SOBOL_total
#'
#' @param amae AMAE index
#'
#' @param amav AMAV index
#'
#' @param amar AMAR index
#'
#' @param amak AMAK index
#'
#' @param dir a directory  to save the results
#'
#' @author
#' Camila Garcia-Echeverri <cagarciae@unal.edu.co> \cr

#'
#' Hydrodynamics of the natural media research group - HYDS
#' National University of Colombia -  Bogota
#'

save_results <- function(SOBOL=NULL, SOBOL_total=NULL, amae=NULL, amav=NULL, amar=NULL, amak=NULL, dir){

  if(missing(SOBOL)){
  }else{
  a <- paste(dir, 'SOBOL.csv', sep = .Platform$file.sep)
  utils::write.csv(SOBOL, file=a)}

  if(missing(SOBOL_total)){
  }else{
  a <- paste(dir,'SOBOL_total.csv', sep = .Platform$file.sep)
  utils::write.csv(SOBOL_total, file=a)}

  if(missing(amae)){
  }else{
  a <- paste(dir,'AMAE.csv', sep = .Platform$file.sep)
  utils::write.csv(amae, file=a)}

  if(missing(amav)){
  }else{
  a <- paste(dir,'AMAV.csv', sep = .Platform$file.sep)
  utils::write.csv(amav, file=a)}

  if(missing(amar)){
  }else{
  a <- paste(dir,'AMAR.csv', sep = .Platform$file.sep)
  utils::write.csv(amar, file=a)}

  if(missing(amak)){
  }else{ a <- "missing"
  a <- paste(dir, 'AMAK.csv', sep = .Platform$file.sep)
  utils::write.csv(amak, file=a)}


  # all_indices=c("SOBOL","AMAE","AMAV", "AMAR",
  #     "AMAK","Mean","Variance","Skewness","Kurtosis")
  # merge_index<-matrix(nrow = pp, ncol=9, dimnames=list(pp_names,all_indices))

  return(invisible())
}
