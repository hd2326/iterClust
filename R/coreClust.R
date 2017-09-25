#' Clustering
#' 
#' A sample clustering framework (described in "Examples" section, used as
#' default in iterClust framework). Customized frameworks can be defined
#' following rules specified in "Usage", "Arguments" and "Value" sections.
#' 
#' @param dset (numeric matrix) features in rows and observations in columns
#' @param iteration (positive integer) specifies current iteration
#' 
#' @return a list, each element contains clustering vectors (named numeric
#' vector with observation names as name and corresponding cluster number as
#' element) under a specific clustering parameter
#' 
#' @keywords coreClust
#' @examples
#' coreClust <- function(dset, iteration){
#'     dist <- as.dist(1 - cor(dset))
#'     range=seq(2, (ncol(dset)-1), by = 1)
#'     clust <- vector("list", length(range))
#'     for (i in 1:length(range)) clust[[i]] <- pam(dist, range[i])$clustering
#'     return(clust)}
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @export

coreClust <- function(dset, iteration){
    dist <- as.dist(1 - cor(dset))
    range=seq(2, (ncol(dset)-1), by = 1)
    clust <- vector("list", length(range))
    for (i in 1:length(range)) clust[[i]] <- pam(dist, range[i])$clustering
    return(clust)}