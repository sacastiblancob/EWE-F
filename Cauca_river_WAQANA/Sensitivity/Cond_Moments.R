#' @name
#' Cond_Moments
#'
#' @title
#' Conditional statistical moments of a model output
#'
#' @description
#' This function evaluates the first four statistical moments
#' after grouping the model output by different parametric
#' ranges.
#'
#'
#' @param parameters_set matrix of dimensions n x pp, where n is the
#' number of runs and pp is the number of parameters.
#'
#' @param out_set matrix of dimensions n x t, where n is the number of
#' runs and t is the number of temporary steps.
#'
#' @param pp_names vector that contains the names of the parameters.
#'
#' @param steps number of divisions of the parametric range.
#'
#'
#' @return
#' A list of arrays, each array has dimensions of steps, t, pp.
#'
#' @export
#'
#' @author
#' Camila Garcia-Echeverri <cagarciae@unal.edu.co> \cr
#' Maria Cristina Areas-Bautista <mcarenasb@unal.edu.co> \cr

#'
#' Hydrodynamics of the natural media research group - HYDS
#' National University of Colombia -  Bogota
#'
#' @examples
#' data("parameters_set", "out_set", "pp_names")
#' \donttest{
#' CM <- Cond_Moments(parameters_set, out_set, pp_names, steps=15)
#' }

Cond_Moments <- function(parameters_set, out_set , pp_names, steps = 100){

  if(!is.matrix(parameters_set) | !is.matrix(out_set)){
    warning("parameters_set and out_set must be matrices")
  }else if (!is.vector(pp_names)) {
    warning("pp_names must be a character vector")
  }else if (dim(parameters_set)[1] != dim(out_set)[1]){
    warning("There is a mismatch in the number of runs between parameters_set and out_set")
  }else{

    pp <- dim(parameters_set)[2]
    t <- dim(out_set)[2]

    CM_count<-array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
    CM_mean<-array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
    CM_var<-array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
    CM_skw<-array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))
    CM_kurt<-array(dim=c(steps,t,pp),dimnames = list(NULL,NULL,pp_names))

    #generates the number of divisions to analyze the model output
    div <- dim(parameters_set)[1]/steps
    sorted <- apply(parameters_set,2,sort)

    ranges_pp <- matrix(nrow=3*steps,ncol=pp)
    colnames(ranges_pp) <- pp_names

    #generates the parameter ranges for each division
    for (s in 1:pp){
      x <- 0
      a <- 1
      for(i in 1:steps){
        x <- (3*i)
        ranges_pp[(x-2),s] <- min(sorted[a:(div*i),s])
        ranges_pp[(x-1),s] <- max(sorted[a:(div*i),s])
        ranges_pp[(x),s] <- ranges_pp[(x-1),s]-ranges_pp[(x-2),s]

        a <- a+div
      }
    }

    ranges_pp<-as.data.frame(ranges_pp)



    for(f in 1:t){

      #calculates conditional moments for each range
      for (s in 1:pp){

        count__<-vector("numeric",0); mean__<-vector("numeric",0); var__<-vector("numeric",0);
        skw__<-vector("numeric",0); kurt__<-vector("numeric",0)

        for(i in 1:steps){

          a<-parameters_set[,s]

          count__[i]<-(sum(a>(ranges_pp[((3*i)-2),s])&a<=(ranges_pp[((3*i)-1),s])))/length(a)

          mean__[i]<-mean(out_set[a>(ranges_pp[((3*i)-2),s])&a<=(ranges_pp[((3*i)-1),s]),f], na.rm = T)

          var__[i]<-stats::var(out_set[a>(ranges_pp[((3*i)-2),s])&a<=(ranges_pp[((3*i)-1),s]),f], na.rm = T)

          skw__[i]<-e1071::skewness(out_set[a>(ranges_pp[((3*i)-2),s])&a<=(ranges_pp[((3*i)-1),s]),f], na.rm = T)

          kurt__[i]<-(e1071::kurtosis(out_set[a>(ranges_pp[((3*i)-2),s])&a<=(ranges_pp[((3*i)-1),s]),f], na.rm = T))+3

        }
        CM_count[,f,s]<-count__
        CM_mean[,f,s]<-mean__
        CM_var[,f,s]<-var__
        CM_skw[,f,s]<-skw__
        CM_kurt[,f,s]<-kurt__
      }
    }
    CM <- list("CM_count"=CM_count, "CM_mean"=CM_mean, "CM_var"=CM_var, "CM_skw"=CM_skw, "CM_kurt"=CM_kurt)


    return(CM)
  }
}


