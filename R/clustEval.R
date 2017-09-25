#' Cluster-wise Clustering Robustness Evaluation
#' 
#' A sample cluster-wise clustering robustness evaluation framework (described
#' in "Examples" section, used as default in iterClust framework). Customized
#' frameworks can be defined following rules specified in "Usage", "Arguments"
#' and "Value" sections.
#' 
#' @param dset (numeric matrix) features in rows and observations in columns
#' @param iteration (positive integer) specifies current iteration
#' @param clust return value of coreClust
#' 
#' @return a numeric vector, specifies the clustering robustness (higher value
#' means more robust) of each clustering scheme
#' 
#' @keywords clustEval
#' @examples
#' clustEval <- function(dset, iteration, clust){
#'     dist <- as.dist(1 - cor(dset))
#'     clustEval <- vector("numeric", length(clust))
#'     for (i in 1:length(clust)){
#'         clustEval[i] <- mean(silhouette(clust[[i]], dist)[, "sil_width"])}
#'     return(clustEval)}
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @export

clustEval <- function(dset, iteration, clust){
    dist <- as.dist(1 - cor(dset))
    clustEval <- vector("numeric", length(clust))
    for (i in 1:length(clust)){
        clustEval[i] <- mean(silhouette(clust[[i]], dist)[, "sil_width"])}
    return(clustEval)}
