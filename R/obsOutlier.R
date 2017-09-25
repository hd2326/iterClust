#' Outlier Observation Evaluation
#' 
#' A sample outlier observation evaluation framework (described in "Examples"
#' section, used as default in iterClust framework). Customized frameworks can
#' be defined following rules specified in "Usage", "Arguments" and "Value"
#' sections.
#' 
#' @param obsEval, return value of obsEval
#' @param iteration (positive integer) specifies current iteration
#' 
#' @return a boolean vector, specifies whether an observation is outlier
#' 
#' @keywords obsOutlier
#' @examples
#' obsOutlier <- function(obsEval, iteration) return(obsEval < 0*iteration-1)
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @export

obsOutlier <- function(obsEval, iteration) return(obsEval < 0*iteration-1)

