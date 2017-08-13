#' Iterative Clustering
#' 
#' A framework for performing clustering analysis iteratively
#' 
#' ####################
#' General Idea
#' ####################
#' 
#' In a scenario where populations A, B1, B2 exist, pronounce differences between A and B may mask subtle differences between B1 and B2. To solve this problem, so that heterogeneity can be better detected, clustering analysis needs to be performed iteratively, so that, for example, in iteration 1, A and B are seperated and in iteration 2, B1 and B2 are seperated.
#' 
#' ####################
#' General Work Flow
#' ####################
#' 
#' Iteration Start ==>>
#' 
#' featureSelectfun (feature selection) ==>>
#' 
#' minFeatureSize (confirm enough features are selected) ==>>
#' 
#' clustHeterofun (confirm heterogeneity) ==>>
#' 
#' clustfun (generate several clustering schemes, only for heterogenous clusters) ==>>
#' 
#' clustEvalfun (pick the optimal clustering scheme) ==>>
#' 
#' minClustSize (remove clusters with few observations) ==>>
#' 
#' obsEvalfun (evaluate how each bservations are clustered) ==>>
#' 
#' obsOutlierfun (remove poorly clustered observations) ==>>
#' 
#' results in Internal Variables ==>>
#' 
#' ith Iteration End
#' 
#' ####################
#' Internal Variables (IV)
#' ####################
#' 
#' iterClust has the following internal variables (IV) which can be used in externally defined functions:
#' 
#' cluster: (list) the return value, described in "Value" section
#' 
#' depth: (numeric) current round of iteration
#' 
#' @param dset (numeric matrix) features in rows and observations in columns
#' @param maxIter (positive integer) specifies maximum iterations to be performed
#' @param minFeatureSize (positive integer) specifies minimum features needed
#' @param featureSelectfun (function) takes a dataset, depth(IV) and cluster$feature(IV), returns a character array, containing features used for clustering analysis
#' @param clustfun (function) takes a dataset and depth(IV), returns a list, containing clustering vectors under different clustering parameters
#' @param minClustSize (positive integer) specifies minimum cluster size
#' @param clustEvalfun (function) takes a dataset, depth(IV) and clustfun result, returns a numeric vector, evaluating the robustness (higher value means more robust) of each clustering scheme
#' @param clustHeterofun (function) takes depth(IV) and clustEvalfun result, returns a boolean vector, deciding whether a cluster is considered as heterogenous
#' @param obsEvalfun (function) takes a dataset and optimal clustfun result determined by clustEvalfun, returns a numeric vector, evaluating the clustering robustness of each observation
#' @param obsOutlierfun (function) takes depth(IV) and obsEvalfun result, returns a boolean vector, deciding whether an observation is outlier
#' 
#' @return a list with the following structure containing iterClust result
#' 
#' --> $cluster (list) --> Iter[i] (list) --> Cluster[j], (character array) names of observations belong to each cluster
#' 
#' |
#'
#' --> $feature (list) --> Iter[i] (list) --> Cluster[j]inIter[i-1], (character array) features used to split each cluster in the previous iteration thereby produce the current clusters
#' 
#' 
#' @keywords iterClust
#' @examples 
#' library(tsne)
#' library(cluster)
#' library(ConsensusClusterPlus)
#' library(bcellViper)
#' 
#' data(bcellViper)
#' exp <- exprs(dset)
#' pheno <- as.character(dset@phenoData@data$description)
#' exp <- exp[, pheno %in% names(table(pheno))[table(pheno) > 5]]
#' pheno <- pheno[pheno %in% names(table(pheno))[table(pheno) > 5]]
#' #load bcellViper expression and phenotype annotation
#' 
#' featureSelectfun <- function(dset, iteration, feature) return(rownames(dset))
#' 
#' clustfun <- function(dset, iteration){
#'   dist <- as.dist(1 - cor(dset))
#'   range=seq(2, (ncol(dset)-1), by = 1)
#'   clust <- vector("list", length(range))
#'   for (i in 1:length(range)) clust[[i]] <- pam(dist, range[i])$clustering
#'   return(clust)}
#'   
#' clustEvalfun <- function(dset, iteration, clust){
#'   dist <- as.dist(1 - cor(dset))
#'   clustEval <- vector("numeric", length(clust))
#'   for (i in 1:length(clust)){
#'       clustEval[i] <- mean(silhouette(clust[[i]], dist)[, "sil_width"])}
#'   return(clustEval)}
#'   
#' clustHeterofun <- function(clustEval, iteration){
#'     return(clustEval > 0*iteration+0.15)}
#' 
#' obsEvalfun <- function(dset, clust, iteration){
#'   dist <- as.dist(1 - cor(dset))
#'   obsEval <- vector("numeric", length(clust))
#'   return(silhouette(clust, dist)[, "sil_width"])}
#' 
#' obsOutlierfun <- function(obsEval, iteration) return(obsEval < 0*iteration-1)
#' #sample external functions for iterClust
#' 
#' c <- iterClust(exp, maxIter=3,
#'                minFeatureSize=100, featureSelectfun=featureSelectfun,
#'                minClustSize=5, clustfun=clustfun,
#'                clustEvalfun=clustEvalfun, clustHeterofun=clustHeterofun,
#'                obsEvalfun=obsEvalfun, obsOutlierfun=obsOutlierfun)
#' #iterClust
#' 
#' set.seed(1)
#' consensusClust = ConsensusClusterPlus(exp, maxK = 10,
#'                                       reps = 100, clusterAlg = "pam",
#'                                       distance = "pearson", plot = FALSE)
#' ICL <- calcICL(consensusClust, plot = FALSE)
#' ICL <- sapply(2:10, function(k, ICL){
#'   s <- ICL$clusterConsensus[grep(k, ICL$clusterConsensus[, "k"]),
#'                             "clusterConsensus"]
#'   mean(s[is.finite(s)])}, ICL=ICL)
#' #consensus clustering
#' 
#' dist <- as.dist(1 - cor(exp))
#' set.seed(1)
#' tsne <- tsne(dist, perplexity = 20, max_iter = 500)
#' #project data on 2D-tSNE space
#' 
#' pdf("bcellViper.pdf", width = 12, height = 6)
#' par(mfrow = c(2, 4), mar = c(4.5, 4.5, 2, 1))
#' for (j in 1:length(c$cluster)){
#'     COL <- structure(rep(1, ncol(exp)), names = colnames(exp))
#'     for (i in 1:length(c$cluster[[j]])) COL[c$cluster[[j]][[i]]] <- i+1
#'     plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
#'          xlab = "Dim1", ylab = "Dim2",
#'          main = paste("iterClust, iter=", j, sep = ""))
#'     text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5, col = COL)
#'     legend("topleft", legend = "Outliers", fill = 1, bty = "n")
#'     plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
#'          xlab = "Dim1", ylab = "Dim2",
#'          main = paste("PAM, k=", length(c$cluster[[j]]), sep = ""))
#'     text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5,
#'          col = pam(dist, k = length(c$cluster[[j]]))$clustering)}
#' #compare iterClust with PAM, which is core clustering function in iterClust
#' 
#' plot(c(2:10), ICL, xlab = "#Clusters", ylab = "Cluster Consensus Score", 
#'      col = c(1, 2, rep(1, 7)), ylim = c(0.8, 1),
#'      cex.lab = 1.5, pch = 16, main = "")
#' plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
#'      xlab = "Dim1", ylab = "Dim2",main = "Consensus Clustering+PAM, k=3")
#' text(tsne[, 1], tsne[, 2], labels = pheno,
#'      cex = 0.5, col = consensusClust[[3]]$consensusClass)
#' plot(c(2:10), ICL, xlab = "#Clusters", ylab = "Cluster Consensus Score",
#'      col = c(rep(1, 5), 2, 1, 1), ylim = c(0.8, 1),
#'      cex.lab = 1.5, pch = 16, main = "")
#' plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
#'      xlab = "Dim1", ylab = "Dim2",main = "Consensus Clustering+PAM, k=7")
#' text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5,
#'      col = consensusClust[[7]]$consensusClass)
#' #compare iterClust with other framework, e.g. consensus clustering
#' 
#' dev.off()
#' #visualize results
#' 
#' @author DING, HONGXU (poulainding@gmail.com)
#' 
#' @export

iterClust <- function(dset, maxIter=10,
                      minFeatureSize=100, featureSelectfun, 
                      minClustSize=10, clustfun, clustEvalfun, clustHeterofun,
                      obsEvalfun, obsOutlierfun){
  depth <- 1
  cluster <- list()
  feature <- list()
  clustHetero <- list()
  #internal variables
  
  feat <- featureSelectfun(dset, depth, feature)#feature selection
  clust <- clustfun(dset[feat, ], depth)#clustering
  clustEval <- clustEvalfun(dset[feat, ], depth, clust)#evaluate each clustering scheme
  clust <- list(clust=structure(clust[[which.max(clustEval)]], names = colnames(dset)),
                clustEval=clustEval[[which.max(clustEval)]],
                featureSelect=feat)#optimal clustering scheme
  clustHetero[[depth]] <- clustHeterofun(clust$clustEval, depth)#whether the dataset is subsetable
  if (clustHetero[[depth]] & length(feat) >= minFeatureSize){
    message(paste("iteration:", depth, sep = ""))
    obsEval <- obsEvalfun(dset[feat, ], clust$clust, depth)#evaluate each observation
    c <- clust$clust[!obsOutlierfun(obsEval, depth)]#remove outliers
    c <- lapply(names(table(c)), function(x, c) names(c)[c == x], c=c)
    cluster[[depth]] <- c[sapply(c, function(x) length(x)) >= minClustSize]
    names(cluster[[depth]]) <- paste("Cluster", 1:length(cluster[[depth]]), sep = "")#organize clusters
    feature[[depth]] <- list(clust$featureSelect)
    names(feature[[depth]]) <- "OriginalDataset"#organize features
  #1st iteration
    
    depth <- depth + 1
    clust <- lapply(cluster[[depth-1]], function(x, dset, depth, feature, featureSelectfun, clustfun, clustEvalfun){
      feat <- featureSelectfun(dset[, x], depth, feature)#feature selection
      clust <- clustfun(dset[feat, x], depth)#clustering
      clustEval <- clustEvalfun(dset[feat, x], depth, clust)#evaluate each clustering scheme
      list(clust=structure(clust[[which.max(clustEval)]], names = x),
           clustEval=clustEval[[which.max(clustEval)]],
           featureSelect=feat)#optimal clustering scheme
    }, dset=dset, depth=depth, feature=feature, featureSelectfun=featureSelectfun, clustfun=clustfun, clustEvalfun=clustEvalfun)
    clustHetero[[depth]] <- clustHeterofun(sapply(clust, function(x) x$clustEval), depth)#whether the dataset is subsetable
    while(sum(clustHetero[[depth]]) > 0 & depth <= maxIter & sum(sapply(clust, function(x, t) length(x$featureSelect) >= minFeatureSize, t=minFeatureSize)) > 0){
      message(paste("iteration:", depth, sep = ""))
      c <- f <- list()
      for (i in 1:length(clust)){
        if (!clustHetero[[depth]][i] | length(clust[[i]]$featureSelect) < minFeatureSize) {c[[i]] <- list(names(clust[[i]]$clust)); f[[i]] <- "NA, Stopped"}
        else{
          obsEval <- obsEvalfun(dset[clust[[i]]$featureSelect, names(clust[[i]]$clust)], clust[[i]]$clust, depth)#evaluate each observation
          c[[i]] <- clust[[i]]$clust[!obsOutlierfun(obsEval, depth)]#remove outliers
          c[[i]] <- lapply(names(table(clust[[i]]$clust)), function(x, c) names(c)[c == x], c=clust[[i]]$clust)
          f[[i]] <- clust[[i]]$featureSelect}
      }
      c <- Reduce(c, c)
      cluster[[depth]] <- c[sapply(c, function(x) length(x)) >= minClustSize]
      names(cluster[[depth]]) <- paste("Cluster", 1:length(cluster[[depth]]), sep = "")#organize clusters
      feature[[depth]] <- f
      names(feature[[depth]]) <- paste("Cluster", 1:length(feature[[depth]]), "inIter", depth-1, sep = "")#organize features
  #ith iteration
      
      depth <- depth + 1
      clust <- lapply(cluster[[depth-1]], function(x, dset, depth, feature, featureSelectfun, clustfun, clustEvalfun){
        feat <- featureSelectfun(dset[, x], depth, feature)#feature selection
        clust <- clustfun(dset[feat, x], depth)#clustering
        clustEval <- clustEvalfun(dset[feat, x], depth, clust)#evaluate each clustering scheme
        list(clust=structure(clust[[which.max(clustEval)]], names = x),
             clustEval=clustEval[[which.max(clustEval)]],
             featureSelect=feat)#optimal clustering scheme
      }, dset=dset, depth=depth, feature=feature, featureSelectfun=featureSelectfun, clustfun=clustfun, clustEvalfun=clustEvalfun)
      clustHetero[[depth]] <- clustHeterofun(sapply(clust, function(x) x$clustEval), depth)#whether the dataset is subsetable
    }
    names(cluster) <- paste("Iter", 1:length(cluster), sep = "")
    names(feature) <- paste("Iter", 1:length(feature), sep = "")
    return(list(cluster=cluster, feature=feature))
  }else{
    message("Homogenous")
    return("Homogenous")
  }
}
