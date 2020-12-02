#' @name
#' Bstat
#'
#' @title
#' Basic statistical measures of a mathematical model results
#'
#' @description
#' This function calculates the mean, variance, skewness, kurtosis and
#' excess kurtosis of a model output, this output can be given for
#' different temporal periods (days, months or years).
#'
#'
#' @param out_set matrix of dimensions n x t, where n equals the number of
#' runs and t is equal to the number of temporary steps.
#'
#'
#' @return
#' a data frame of dimensions t x 6, here t is the number of temporary steps
#' and each column corresponds to a statistical measure: mean, variance,
#' skewness, kurtosis and excess kurtosis.
#'
#' @export
#'
#' @author
#' Camila Garcia-Echeverri <cagarciae@unal.edu.co> \cr
#'
#' Hydrodynamics of the natural media research group - HYDS
#' National University of Colombia -  Bogota
#'
#' @examples
#' data("out_set")
#' data_Bstat <- Bstat(out_set)
#'

Bstat <- function(out_set){

  if( !is.matrix(out_set)){
    warning("out_set must be a matrix")
  }else{


    sufix <- c(Mean="mean",Variance="var",Skewness="skw",Kurtosis="kurt", Kurtosis_ex="kurt_ex")

    for(i in 1:5){
      a <- paste(sufix[i],"_t",sep="")
      aa <- vector("numeric",0)
      assign(a,aa)
    }

    t <- dim(out_set)[2]

    #calculates the first four statistical moments
    for(f in 1:t){
      mean_t[f] <- mean(out_set[,f])
      var_t[f] <- stats::var(out_set[,f])
      skw_t[f] <- e1071::skewness(out_set[,f])
      kurt_t[f] <- mean((out_set[,f]-mean_t[f])^4)/mean((out_set[,f]-mean_t[f])^2)^2
      kurt_ex_t[f] <- e1071::kurtosis(out_set[,f])
    }

    data_Bstat <- data.frame(t<-c(1:t), Media_t<-mean_t,
                            Var_t<-var_t, Skw_t<-skw_t, Kurt_t<-kurt_t,
                            Kurt_ex_t<-kurt_ex_t)

    return(data_Bstat)
  }
}
