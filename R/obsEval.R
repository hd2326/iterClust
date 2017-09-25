#' Observation-wise Clustering Robustness Evaluation
#' 
#' A sample observation-wise clustering robustness evaluation framework
#' (described in "Examples" section, used as default in iterClust framework).
#' Customized frameworks can be defined following rules specified in "Usage",
#' "Arguments" and "Value" sections.
#' 
#' @param dset (numeric matrix) features in rows and observations in columns
#' @param clust optimal return value of coreClust
#' @param iteration (positive integer) specifies current iteration
#' 
#' @return a numeric vector, specifies the clustering robustness (higher value
#' means more robust) of each observation under the optimal clustering scheme
#' 
#' @keywords obsEval
#' @examples
#' obsEval <- function(dset, clust, iteration){
#'     dist <- as.dist(1 - cor(dset))
#'     obsEval <- vector("numeric", length(clust))
#'     return(silhouette(clust, dist)[, "sil_width"])}
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom cluster silhouette
#' 
#' @export

obsEval <- function(dset, clust, iteration){
    dist <- as.dist(1 - cor(dset))
    obsEval <- vector("numeric", length(clust))
    return(silhouette(clust, dist)[, "sil_width"])}
