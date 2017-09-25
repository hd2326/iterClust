#' Cluster Heterogeneity Evaluation
#' 
#' A sample cluster heterogeneity evaluation framework (described in "Examples"
#' section, used as default in iterClust framework). Customized frameworks can
#' be defined following rules specified in "Usage", "Arguments" and "Value"
#' sections.
#' 
#' @param clustEval, return value of clustEval
#' @param iteration (positive integer) specifies current iteration
#' 
#' @return a boolean vector, specifies whether clusters are heterogenous
#' 
#' @keywords clustHetero
#' @examples
#' clustHetero <- function(clustEval, iteration){
#'     return(clustEval > 0*iteration+0.15)}
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @export

clustHetero <- function(clustEval, iteration){
    return(clustEval > 0*iteration+0.15)}

