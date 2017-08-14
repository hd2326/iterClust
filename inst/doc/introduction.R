### R code from vignette source 'introduction.rnw'

###################################################
### code chunk number 1: dataprep
###################################################
library(iterClust)
library(bcellViper)
data(bcellViper)
exp <- exprs(dset)
pheno <- as.character(dset@phenoData@data$description)
exp <- exp[, pheno %in% names(table(pheno))[table(pheno) > 5]]
pheno <- pheno[pheno %in% names(table(pheno))[table(pheno) > 5]]

dim(exp)

table(pheno)


###################################################
### code chunk number 2: cluster
###################################################
library(cluster)


###################################################
### code chunk number 3: featureSelectfun
###################################################
featureSelectfun <- function(dset, iteration, feature) return(rownames(dset))


###################################################
### code chunk number 4: clustfun
###################################################
clustfun <- function(dset, iteration){
    dist <- as.dist(1 - cor(dset))
    range=seq(2, (ncol(dset)-1), by = 1)
    clust <- vector("list", length(range))
    for (i in 1:length(range)) clust[[i]] <- pam(dist, range[i])$clustering
    return(clust)}


###################################################
### code chunk number 5: clustEvalfun
###################################################
clustEvalfun <- function(dset, iteration, clust){
    dist <- as.dist(1 - cor(dset))
    clustEval <- vector("numeric", length(clust))
    for (i in 1:length(clust)){
        clustEval[i] <- mean(silhouette(clust[[i]], dist)[, "sil_width"])}
    return(clustEval)}


###################################################
### code chunk number 6: clustHeterofun
###################################################
clustHeterofun <- function(clustEval, iteration){
    return(clustEval > 0*iteration+0.15)}


###################################################
### code chunk number 7: obsEvalfun
###################################################
obsEvalfun <- function(dset, clust, iteration){
    dist <- as.dist(1 - cor(dset))
    obsEval <- vector("numeric", length(clust))
    return(silhouette(clust, dist)[, "sil_width"])}


###################################################
### code chunk number 8: obsOutlierfun
###################################################
obsOutlierfun <- function(obsEval, iteration) return(obsEval < 0*iteration-1)


###################################################
### code chunk number 9: iterClust
###################################################
c <- iterClust(exp, maxIter=3,
               minFeatureSize=100, featureSelectfun=featureSelectfun,
               minClustSize=5, clustfun=clustfun,
               clustEvalfun=clustEvalfun, clustHeterofun=clustHeterofun,
               obsEvalfun=obsEvalfun, obsOutlierfun=obsOutlierfun)

names(c)

names(c$cluster)

names(c$cluster$Iter1)

c$cluster$Iter1$Cluster1

names(c$feature)

names(c$feature$Iter1)

names(c$feature$Iter2)

c$feature$Iter2$Cluster1inIter1[1:10]


###################################################
### code chunk number 10: consensusclust
###################################################
library(ConsensusClusterPlus)
set.seed(1)
consensusClust = ConsensusClusterPlus(exp, maxK = 10,
                                      reps = 100, clusterAlg = "pam",
                                      distance = "pearson", plot = FALSE)
ICL <- calcICL(consensusClust, plot = FALSE)
ICL <- sapply(2:10, function(k, ICL){
    s <- ICL$clusterConsensus[grep(k, ICL$clusterConsensus[, "k"]),
                              "clusterConsensus"]
    mean(s[is.finite(s)])}, ICL=ICL)


###################################################
### code chunk number 11: tsne
###################################################
library(tsne)
dist <- as.dist(1 - cor(exp))
set.seed(1)
tsne <- tsne(dist, perplexity = 20, max_iter = 500)


###################################################
### code chunk number 12: fig1
###################################################
par(mfrow = c(1, 2))
for (j in 1:length(c$cluster)){
    COL <- structure(rep(1, ncol(exp)), names = colnames(exp))
    for (i in 1:length(c$cluster[[j]])) COL[c$cluster[[j]][[i]]] <- i+1
    plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
         xlab = "Dim1", ylab = "Dim2",
         main = paste("iterClust, iter=", j, sep = ""))
    text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5, col = COL)
    legend("topleft", legend = "Outliers", fill = 1, bty = "n")}


###################################################
### code chunk number 13: fig2
###################################################
par(mfrow = c(1, 2))
for (j in 1:length(c$cluster)){
    plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
         xlab = "Dim1", ylab = "Dim2",
         main = paste("PAM, k=", length(c$cluster[[j]]), sep = ""))
    text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5,
         col = pam(dist, k = length(c$cluster[[j]]))$clustering)}


###################################################
### code chunk number 14: fig3
###################################################
par(mfrow = c(2, 2))
plot(c(2:10), ICL, xlab = "#Clusters", ylab = "Cluster Consensus Score", 
     col = c(1, 2, rep(1, 7)), ylim = c(0.8, 1),
     cex.lab = 1.5, pch = 16, main = "")
plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
     xlab = "Dim1", ylab = "Dim2",main = "Consensus Clustering+PAM, k=3")
text(tsne[, 1], tsne[, 2], labels = pheno,
     cex = 0.5, col = consensusClust[[3]]$consensusClass)
plot(c(2:10), ICL, xlab = "#Clusters", ylab = "Cluster Consensus Score",
     col = c(rep(1, 5), 2, 1, 1), ylim = c(0.8, 1),
     cex.lab = 1.5, pch = 16, main = "")
plot(tsne[, 1], tsne[, 2], cex = 0, cex.lab = 1.5,
     xlab = "Dim1", ylab = "Dim2",main = "Consensus Clustering+PAM, k=7")
text(tsne[, 1], tsne[, 2], labels = pheno, cex = 0.5,
     col = consensusClust[[7]]$consensusClass)

