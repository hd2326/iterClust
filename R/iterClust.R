#' Iterative Clustering
#' 
#' A framework for performing clustering analysis iteratively
#' 
#' ####################
#' General Idea
#' ####################
#' 
#' In a scenario where populations A, B1, B2 exist, pronounce differences
#' between A and B may mask subtle differences between B1 and B2. To solve this
#' problem, so that heterogeneity can be better detected, clustering analysis
#' needs to be performed iteratively, so that, for example, in iteration 1, A
#' and B are seperated and in iteration 2, B1 and B2 are seperated.
#' 
#' ####################
#' General Work Flow
#' ####################
#' 
#' ith Iteration Start ==>>
#' 
#' featureSelect (feature selection) ==>>
#' 
#' minFeatureSize (confirm enough features are selected) ==>>
#' 
#' clustHetero (confirm heterogeneity) ==>>
#' 
#' coreClust (generate several clustering schemes to be evaluated) ==>>
#' 
#' clustEval (pick optimal clustering scheme generated in previous step) ==>>
#' 
#' minClustSize (remove clusters with few observations) ==>>
#' 
#' obsEval (evaluate how each observations are clustered) ==>>
#' 
#' obsOutlier (remove poorly clustered observations) ==>>
#' 
#' results in Internal Variables (IV) ==>>
#' 
#' ith Iteration End
#' 
#' ####################
#' Internal Variables (IV)
#' ####################
#' 
#' The following IVs are used in user-defined functions in each iteration:
#' 
#' cluster: (list) the return value, described in "Value" section
#' 
#' depth: (numeric) current round of iteration
#' 
#' @param dset (numeric matrix or data.frame) features in rows and observations
#' in columns, or SummarizedExperiment0 and ExpressionSet object
#' @param maxIter (positive integer) specifies maximum number iterations to be
#' performed
#' @param minFeatureSize (positive integer) specifies minimum number of features
#' needed
#' @param featureSelect (function) takes a dataset, depth(IV) and
#' cluster$feature(IV), returns a character array, containing features used for
#' clustering analysis
#' @param coreClust (function) takes a dataset and depth(IV), returns a list,
#' containing clustering vectors under different clustering parameters
#' @param minClustSize (positive integer) specifies minimum cluster size
#' @param clustEval (function) takes a dataset, depth(IV) and coreClust result,
#' returns a numeric vector, evaluating the robustness (higher value means more
#' robust) of each clustering scheme
#' @param clustHetero (function) takes depth(IV) and clustEval result, returns a
#' boolean vector, deciding whether a cluster is considered as heterogenous
#' @param obsEval (function) takes a dataset and optimal coreClust result
#' determined by clustEval, returns a numeric vector, evaluating the clustering
#' robustness of each observation
#' @param obsOutlier (function) takes depth(IV) and obsEval result, returns
#' a boolean vector, deciding whether an observation is outlier
#' 
#' @return a list with the following structure containing iterClust result
#' 
#' --> $cluster (list) $Iter[i] (list) $Cluster[j],
#' (character array) names of observations belong to each cluster
#'
#' --> $feature (list) $Iter[i] (list) $Cluster[j]inIter[i-1],
#' (character array) features used to split each cluster in the previous
#' iteration thereby produce the current clusters
#' 
#' 
#' @keywords iterClust
#' @examples 
#' library(tsne)
#' library(cluster)
#' library(bcellViper)
#' 
#' data(bcellViper)
#' exp <- exprs(dset)
#' pheno <- as.character(dset@phenoData@data$description)
#' exp <- exp[, pheno %in% names(table(pheno))[table(pheno) > 5]]
#' pheno <- pheno[pheno %in% names(table(pheno))[table(pheno) > 5]]
#' #load bcellViper expression and phenotype annotation
#' 
#' c <- iterClust(exp, maxIter=3, minClustSize=5)
#' #iterClust
#' 
#' dist <- as.dist(1 - cor(exp))
#' set.seed(1)
#' tsne <- tsne(dist, perplexity = 20, max_iter = 500)#' 
#' for (j in 1:length(c$cluster)){
#'     COL <- structure(rep(1, ncol(exp)), names = colnames(exp))
#'     for (i in 1:length(c$cluster[[j]])) COL[c$cluster[[j]][[i]]] <- i+1
#'     plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
#'          xlab = "Dim1", ylab = "Dim2",
#'          main = paste("iterClust, iter=", j, sep = ""))
#'     text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5, col = COL)
#'     legend("topleft", legend = "Outliers", fill = 1, bty = "n")}
#' #visualize results
#' 
#' @author DING, HONGXU (hd2326@columbia.edu)
#' 
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
#' @importFrom Biobase exprs
#' 
#' @export

iterClust <- function(dset, maxIter=10,
                      minFeatureSize=100,
                      featureSelect=iterClust::featureSelect, 
                      minClustSize=10,
                      coreClust=iterClust::coreClust,
                      clustEval=iterClust::clustEval,
                      clustHetero=iterClust::clustHetero,
                      obsEval=iterClust::obsEval,
                      obsOutlier=iterClust::obsOutlier){
    
    if(is.matrix(dset)) dset <- dset
    else if (is.data.frame(dset)) dset <- as.matrix(dset)
    else if(methods::is(dset, "SummarizedExperiment")) dset <- assay(dset)
    else if(methods::is(dset, "ExpressionSet")) dset <- exprs(dset)
    else message("dset should be matrix, SummarizedExperiment or ExpressionSet")
    #class regularization
    
    depth <- 1
    cluster <- list()
    feature <- list()
    clustHetero <- list()
    #internal variables
    
    feat <- featureSelect(dset, depth, feature)
    #feature selection
    clust <- coreClust(dset[feat, ], depth)
    #clustering
    clustEval <- clustEval(dset[feat, ], depth, clust)
    #evaluate each clustering scheme
    clust <- list(clust=structure(clust[[which.max(clustEval)]],
                                  names = colnames(dset)),
                  clustEval=clustEval[[which.max(clustEval)]],
                  featureSelect=feat)
    #optimal clustering scheme
    clustHetero[[depth]] <- clustHetero(clust$clustEval, depth)
    #whether the dataset is subsetable
    if (clustHetero[[depth]] & length(feat) >= minFeatureSize){
        message(paste("iteration:", depth, sep = ""))
        obsEval <- obsEval(dset[feat, ], clust$clust, depth)
        #evaluate each observation
        cc <- clust$clust[!obsOutlier(obsEval, depth)]
        #remove outliers
        cc <- lapply(names(table(cc)), function(x, cc){
            names(cc)[cc == x]}, cc=cc)
        cluster[[depth]] <- cc[unlist(lapply(cc, function(x){
            length(x) >= minClustSize}))]
        names(cluster[[depth]]) <- paste("Cluster",
                                         seq_len(length(cluster[[depth]])),
                                         sep = "")
        #organize clusters
        feature[[depth]] <- list(clust$featureSelect)
        names(feature[[depth]]) <- "OriginalDataset"
        #organize features
        #1st iteration
        
        depth <- depth + 1
        clust <- lapply(cluster[[depth-1]], function(x, dset, depth, feature,
                                                     featureSelect, coreClust,
                                                     clustEval){
            feat <- featureSelect(dset[, x], depth, feature)
            #feature selection
            clust <- coreClust(dset[feat, x], depth)
            #clustering
            clustEval <- clustEval(dset[feat, x], depth, clust)
            #evaluate each clustering scheme
            list(clust=structure(clust[[which.max(clustEval)]], names = x),
                 clustEval=clustEval[[which.max(clustEval)]],
                 featureSelect=feat)
            #optimal clustering scheme
        }, dset=dset, depth=depth, feature=feature,
        featureSelect=featureSelect, coreClust=coreClust,
        clustEval=clustEval)
        clustHetero[[depth]] <- clustHetero(unlist(lapply(clust, function(x){
            x$clustEval})),
            depth)
        #whether the dataset is subsetable
        while(sum(clustHetero[[depth]]) > 0 & 
              depth <= maxIter & 
              sum(unlist(lapply(clust, function(x, t){
                  length(x$featureSelect) >= minFeatureSize
              }, t=minFeatureSize))) > 0){
            message(paste("iteration:", depth, sep = ""))
            cc <- f <- list()
            for (i in seq_len(length(clust))){
                if (!clustHetero[[depth]][i] | 
                    length(clust[[i]]$featureSelect) < minFeatureSize){
                    cc[[i]] <- list(names(clust[[i]]$clust))
                    f[[i]] <- "NA, Stopped"}
                else{
                    obsEval <- obsEval(
                        dset[clust[[i]]$featureSelect, names(clust[[i]]$clust)],
                        clust[[i]]$clust,
                        depth)
                    #evaluate each observation
                    cc[[i]] <- clust[[i]]$clust[!obsOutlier(obsEval, depth)]
                    #remove outliers
                    cc[[i]] <- lapply(names(table(clust[[i]]$clust)),
                                      function(x, cc) names(cc)[cc == x],
                                      cc=clust[[i]]$clust)
                    f[[i]] <- clust[[i]]$featureSelect}
            }
            cc <- Reduce(c, cc)
            cluster[[depth]] <- cc[unlist(lapply(cc, function(x){
                length(x) >= minClustSize}))]
            names(cluster[[depth]]) <- paste("Cluster",
                                             seq_len(length(cluster[[depth]])),
                                             sep = "")
            #organize clusters
            feature[[depth]] <- f
            names(feature[[depth]]) <- paste("Cluster",
                                             seq_len(length(feature[[depth]])),
                                             "inIter",
                                             depth-1,
                                             sep = "")
            #organize features
            #ith iteration
            
            depth <- depth + 1
            clust <- lapply(cluster[[depth-1]], function(x, dset, depth,
                                                         feature, featureSelect,
                                                         coreClust, clustEval){
                feat <- featureSelect(dset[, x], depth, feature)
                #feature selection
                clust <- coreClust(dset[feat, x], depth)
                #clustering
                clustEval <- clustEval(dset[feat, x], depth, clust)
                #evaluate each clustering scheme
                list(clust=structure(clust[[which.max(clustEval)]], names = x),
                     clustEval=clustEval[[which.max(clustEval)]],
                     featureSelect=feat)
                #optimal clustering scheme
            }, dset=dset, depth=depth,
            feature=feature, featureSelect=featureSelect,
            coreClust=coreClust, clustEval=clustEval)
            clustHetero[[depth]] <- clustHetero(
                unlist(lapply(clust, function(x) x$clustEval)),
                depth)
            #whether the dataset is subsetable
        }
        names(cluster) <- paste("Iter", seq_len(length(cluster)), sep = "")
        names(feature) <- paste("Iter", seq_len(length(feature)), sep = "")
        return(list(cluster=cluster, feature=feature))
    }else{
        message("Homogenous")
        return("Homogenous")
    }
}
