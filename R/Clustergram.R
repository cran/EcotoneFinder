######################################## Clustergram functions ############################################

######################################## Basic clustergram ################################################

#' Clustergram base function
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations
#' @param k.range A vector with the number of clusters to plot the clustergram
#'   for.
#' @param clustering.function Which clustering method to be used. Default is
#'   k-means. Can be FCM is set to clustergram.vegclust. See details
#' @param clustergram.plot Type of plot for the output. See details.
#' @param line.width Graphical parameter. Width of the lines.
#' @param add.center.points Logical. Should the cluster means be plotted (as
#'   points).
#' @param ... Additional arguments to be passed to the clustering function.
#'
#' @details This is the clustergram function created by Matthias Schonlau. See:
#'   Schonlau M. The clustergram: A graph for visualizing hierarchical and
#'   nonhierarchical cluster analyses. The Stata Journal. 2002;2:391–402.
#'
#'   It is reproduced in this package for convenience. This package also provide
#'   extensions of the clustergram method for fuzzy-c-means clustering and for
#'   the evolution of the main fuzzy indices. These extensions take the form of
#'   additional options to be passed in the clustering.function argument and the
#'   clustergram.plot argument.
#'
#'   It is also recommended to run the clustergram analysis several times and
#'   compare the obtained outputs, as they may vary significantly.
#'
#' @return A clustergram plot of the inputed data
#'
#' @export
#'
#' @examples
#'    ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## 6 clustergram plots
#'    for (i in 1:6) clustergram(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                               k.range = 2:10, line.width = .2)
#'
clustergram <- function(Data, k.range = 2:10 ,
                        clustering.function = clustergram.kmeans,
                        clustergram.plot = clustergram.plot.matlines,
                        line.width = .004, add.center.points = TRUE, ...)
{
  n <- dim(Data)[1]
  PCA.1 <- Data %*% stats::princomp(Data)$loadings[,1]
  COL <- colorspace::heat_hcl(n)[order(PCA.1)]
  line.width <- rep(line.width, n)
  Y <- NULL
  X <- NULL
  centers.points <- list()
  for(k in k.range)
  {
    k.clusters <- clustering.function(Data, k, ...)
    clusters.vec <- k.clusters$cluster
    the.centers <- k.clusters$centers
    noise <- unlist(tapply(line.width, clusters.vec, cumsum))[order(seq_along(clusters.vec)[order(clusters.vec)])]
    y <- the.centers[clusters.vec] + noise
    Y <- cbind(Y, y)
    x <- rep(k, length(y))
    X <- cbind(X, x)
    centers.points[[k]] <- data.frame(y = the.centers , x = rep(k , k))
  }
  x.range <- range(k.range)
  y.range <- range(PCA.1)
  clustergram.plot(X,Y, k.range,
                   x.range, y.range , COL,
                   add.center.points , centers.points)
}

################################# Clustergram.kmeans base function: ############################

#' Type function that clustergram takes for clustering.
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations
#' @param k Number of desired groups for the k-means clustering.
#' @param ... Additional parameters to be passed in the kmeans function (from
#'   the stats package).
#'
#' @return A list containing the cluster vector and the centers matrix (see
#'   kmeans function).
#'
#' @details This is the type of function that the clustergram function uses for
#'   clustering. The return list is internally used by the clustergram to build
#'   the clustergram plot.
#'
#' @export
#'
#' @examples
#'   ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## 6 clustergram plots
#'    for (i in 1:6) clustergram(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                               clustering.function = clustergram.kmeans,
#'                               k.range = 2:10, line.width = .2)
#'
clustergram.kmeans <- function(Data, k, ...)
{
  cl <- stats::kmeans(Data, k,...)
  cluster <- cl$cluster
  centers <- cl$centers %*% stats::princomp(Data)$loadings[,1]
  return(list(
    cluster = cluster,
    centers = centers
  ))
}

#################################### New clustering fuction for FCM ########################################

#' Vegclust function for clustergram
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations.
#' @param k Number of desired groups for the FCM clustering.
#' @param method Clustering method for the vegclust function.
#' @param ... Additional parameters to be passed to the vegclust function.
#'
#' @return A list containing the cluster vector and the centers matrix (see
#'   vegclust function).
#'
#' @details This is an implementation of Fuzzy c-means clustering (with the
#'   vegclust function of the vegclust package) for the clustergram function.
#'   The return list is internally used by the clustergram to build the
#'   clustergram plot.
#'
#' @export
#'
#' @examples
#' \donttest{
#'   ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## clustergram plot
#'    clustergram(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                               clustering.function = clustergram.vegclust,
#'                               k.range = 2:10, line.width = .2)
#'  }
#'
clustergram.vegclust <- function(Data, k, method = method, ...)
{
  cl <- vegclust::vegclust(Data, mobileCenters = k, method = "FCM", ...)
  colnames(cl$memb) <- c(1:length(cl$memb))
  cl.crisp <- as.integer(colnames(cl$memb)[max.col(cl$memb, ties.method = "first")])
  names(cl.crisp) <- rownames(cl$memb)
  cluster <- cl.crisp
  centers <- as.matrix(cl$mobileCenters) %*% stats::princomp(Data)$loadings[,1]
  return(list(
    cluster = cluster,
    centers = centers
  ))
}

############################# New clustering function for FCM and fuzzy indices ##########################

#' Vegclust clustering with fuzzy indices computation for clustergram
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations.
#' @param k Number of desired groups for the FCM clustering.
#' @param method Clustering method for the vegclust function.
#' @param ... Additional parameters to be passed to the vegclust function.
#'
#' @return A list containing the cluster vector, the centers matrix and a vector
#'   of four fuzzy indices (partition coefficient (PC), normalized partition
#'   coefficient (PCN), partition entropy (PE) and normalized partition entropy
#'   (PEN)). See vegclust and veclustIndex functions.
#'
#' @details Additionally to the FCM clustering, the function compute the main
#'   fuzzy indices to help with the decision on the optimal number of cluster in
#'   the data. Maximum values of PCN or minimum values of PEN can be used as
#'   criteria to choose the number of clusters.
#'
#' @export
#'
#' @examples
#' \donttest{
#'   ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## clustergram plots with fuzzy indices plots:
#'    clustergramInd(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                                  clustering.function = clustergram.vegclust.Ind,
#'                                  clustergram.plot = clustergram.plot.matlines,
#'                                  FuzzyIndice.plot = FuzzyIndice.plot.matlines,
#'                                  k.range = 2:10, line.width = .2)
#'  }
#'
clustergram.vegclust.Ind <- function(Data, k, method = "FCM", ...)
{
  cl <- vegclust::vegclust(Data, mobileCenters = k, method = method)
  colnames(cl$memb) <- c(1:length(cl$memb))
  cl.crisp <- as.integer(colnames(cl$memb)[max.col(cl$memb, ties.method = "first")])
  names(cl.crisp) <- rownames(cl$memb)
  ind <- vegclust::vegclustIndex(cl)
  Ind <- data.frame(Indice = names(ind), Value = ind, row.names = NULL)
  cluster <- cl.crisp
  centers <- as.matrix(cl$mobileCenters) %*% stats::princomp(Data)$loadings[,1]
  return(list(
    cluster = cluster,
    centers = centers,
    Indices = Ind
  ))
}

################################ New clustering fuction for c-means #########################

#' cmeans function for clustergram
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations.
#' @param k Number of desired groups for the c-means clustering.
#' @param method Clustering method for the cmeans function.
#' @param ... Additional parameters to be passed to the cmeans function.
#'
#' @details This is an implementation of Fuzzy c-means clustering (with the
#'   cmeans function of the e1071 package) for the clustergram function. The
#'   return list is internally used by the clustergram to build the clustergram
#'   plot.
#'
#' @return A list containing the cluster vector and the centers matrix (see
#'   cmeans function).
#'
#' @export
#'
#' @examples
#' \donttest{
#' ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## 6 clustergram plots
#'    for (i in 1:6) clustergram(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                               clustering.function = clustergram.cmeans,
#'                               k.range = 2:10, line.width = .2)
#'  }
#'
clustergram.cmeans <- function(Data, k, method = "cmeans", ...)
{
  cl <- e1071::cmeans(Data, centers = k, method = method, ...)
  cl.crisp <- as.integer(colnames(cl$membership)[max.col(cl$membership, ties.method = "first")])
  names(cl.crisp) <- rownames(cl$membership)
  cluster <- cl.crisp
  centers <- as.matrix(cl$centers) %*% stats::princomp(Data)$loadings[,1]
  return(list(
    cluster = cluster,
    centers = centers
  ))
}

######################## New clustering fuction for c-means and fuzzy indices #########################

#' cmeans clustering with fuzzy indices computation for clustergram
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations.
#' @param k Number of desired groups for the FCM clustering.
#' @param method Clustering method for the cmeans function.
#' @param ... Additional parameters to be passed to the cmeans function.
#'
#' @details Additionally to the FCM clustering with the cmeans function (e1071
#'   package), the function compute the main fuzzy indices to help with the
#'   decision on the optimal number of cluster in the data. The indices are
#'   computed with the vegclustIndex function of the vegclust package. Maximum
#'   values of PCN or minimum values of PEN can be used as criteria to choose
#'   the number of clusters.
#'
#' @return A list containing the cluster vector, the centers matrix and a vector
#'   of four fuzzy indices (partition coefficient (PC), normalized partition
#'   coefficient (PCN), partition entropy (PE) and normalized partition entropy
#'   (PEN)). See vegclust and veclustIndex functions.
#'
#' @export
#'
#' @examples
#' \donttest{
#'   ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## clustergram plots with fuzzy indices plots:
#'    clustergramInd(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                                  clustering.function = clustergram.cmeans.Ind,
#'                                  clustergram.plot = clustergram.plot.matlines,
#'                                  FuzzyIndice.plot = FuzzyIndice.plot.matlines,
#'                                  k.range = 2:10, line.width = .2)
#'  }
#'
clustergram.cmeans.Ind <- function(Data, k, method = "cmeans", ...)
{
  cl <- e1071::cmeans(Data, centers = k, method = method, ...)
  cl.crisp <- as.integer(colnames(cl$membership)[max.col(cl$membership, ties.method = "first")])
  names(cl.crisp) <- rownames(cl$membership)
  ind <- vegclust::vegclustIndex(cl$membership)
  Ind <- data.frame(Indice = names(ind), Value = ind, row.names = NULL)
  cluster <- cl.crisp
  centers <- as.matrix(cl$centers) %*% stats::princomp(Data)$loadings[,1]
  return(list(
    cluster = cluster,
    centers = centers,
    Indices = Ind
  ))
}

################################## Basic clustergram plot function ###########################

#' Plot function for clustergram
#'
#' @param X vector of the k number of cluster for the x-axis
#' @param Y corrdinates of the cluster centers on the y-axis
#' @param k.range x axis breaks.
#' @param x.range x axis range.
#' @param y.range y axis range (PCA scores).
#' @param COL colour palette.
#' @param add.center.points Should the centers be plotted as points. (see
#'   clustergram)
#' @param centers.points matrix of centers position.
#'
#' @details Internal clustergram plot function. The input arguments are computed
#'   by the clustergram function directly.
#'
#' @return A clustergram plot.
#'
#' @export
#'
#' @examples
#' ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## 6 clustergram plots
#'    for (i in 1:6) clustergram(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                               clustering.function = clustergram.kmeans,
#'                               clustergram.plot = clustergram.plot.matlines,
#'                               k.range = 2:10, line.width = .2)
#'
clustergram.plot.matlines <- function(X,Y, k.range,
                                      x.range, y.range , COL,
                                      add.center.points , centers.points)
{
  graphics::plot(0,0, col = "white", xlim = x.range, ylim = y.range,
       axes = FALSE,
       xlab = "Number of clusters (k)", ylab = "PCA weighted Mean of the clusters", main = c("Clustergram of the PCA-weighted Mean of" ,"the clusters k-mean clusters vs number of clusters (k)"))
  graphics::axis(side = 1, at = k.range)
  graphics::axis(side = 2)
  graphics::abline(v = k.range, col = "grey")
  graphics::matlines(t(X), t(Y), pch = 19, col = COL, lty = 1, lwd = 1.5)
  if(add.center.points)
  {
    xx <- plyr::ldply(centers.points, rbind)
    graphics::points(xx$y~xx$x, pch = 19, col = "red", cex = 1.3)
  }
}


################################## clustergram plot adaptation to fuzzy indices ############################

############################## Internal plot function :

#' Plot function for fuzzy indices with clustergram.
#'
#' @param Z Fuzzy indices matrix.
#' @param k.range x axis breaks.
#' @param x.range x axis range.
#' @param z.range y axis range for the fuzzy indices.
#'
#' @return A plot with the evolution of the fuzzy indices given the number of
#'   fuzzy clusters that were applied to the data.
#'
#' @details This function provide the tools to add a fuzzy indices evolution
#'   plot together with the normal clustegram plot with the evolution of the
#'   relative positions of the cluster centers.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## clustergram plots with fuzzy indices plots:
#'    clustergramInd(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                                  clustering.function = clustergram.vegclust.Ind,
#'                                  clustergram.plot = clustergram.plot.matlines,
#'                                  FuzzyIndice.plot = FuzzyIndice.plot.matlines,
#'                                  k.range = 2:10, line.width = .2)
#'  }
#'
FuzzyIndice.plot.matlines = function(Z, k.range,
                                     x.range, z.range)
{
  graphics::plot(Z[,1], Z[,2], col = "white", xlim = x.range, ylim = z.range,
       xlab = "Number of clusters (k)", ylab = "Fuzzy Partition Entropy", main = c("Partition Entropy decrease" ,"with increasing number of clusters (k)"),
       pch = NA)
  graphics::abline(v = k.range, col = "grey")
  graphics::axis(side = 1, at = k.range)
  graphics::points(Z[,1], Z[,2], pch = 16, col = "black")
  graphics::lines(Z[,1], Z[,2], pch = 16, col = "red", lty = 1, lwd = 1.5)
  graphics::points(Z[,1], Z[,3], pch = 16, col = "black")
  graphics::lines(Z[,1], Z[,3], pch = 16, col = "blue", lty = 1, lwd = 1.5)
}

############################# Clustergram plot parameters :

#' Clustergram with fuzzy indices plot
#'
#' @param Data Should be a scales matrix. Where each column belongs to a
#'   different dimension of the observations.
#' @param k.range A vector with the number of clusters to plot the clustergram
#'   for.
#' @param clustering.function Which clustering method to be used. Default is
#'   k-means. Can be FCM is set to clustergram.vegclust. See details
#' @param clustergram.plot Type of plot for the clustergram output. See details.
#' @param FuzzyIndice.plot Type of plot for the fuzzy indices output. See
#'   details.
#' @param line.width Graphical parameter. Width of the lines.
#' @param add.center.points Logical. Should the cluster means be plotted (as
#'   points).
#' @param ... Additional arguments to be passed to the clustering function.
#'
#' @details This clustergram fuction produces an additional plot with the
#'   evolution of the main fuzzy indices (normalized partition coefficient (PCN)
#'   and normalized partition entropy (PEN)). Maximum values of PCN or minimum
#'   values of PEN can be used as criteria to choose the number of clusters.
#'
#' @return A clustergram plot and a fuzzy indices evolution plot of the inputed
#'   data
#'
#' @export
#'
#' @examples
#' \donttest{
#' ####### Example data:
#'    SyntheticTrial <- SyntheticData(SpeciesNum = 100,
#'                                    CommunityNum = 3, SpCo = NULL,
#'                                    Length = 500,
#'                                    Parameters = list(a=c(40, 80, 50),
#'                                                      b=c(100,250,400),
#'                                                      c=rep(0.03,3)),
#'                                    dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'    ######## clustergram plots with fuzzy indices plots:
#'    clustergramInd(as.matrix(SyntheticTrial[,2:ncol(SyntheticTrial)]),
#'                                  clustering.function = clustergram.vegclust.Ind,
#'                                  clustergram.plot = clustergram.plot.matlines,
#'                                  FuzzyIndice.plot = FuzzyIndice.plot.matlines,
#'                                  k.range = 2:10, line.width = .2)
#'  }
#'
clustergramInd = function(Data, k.range = 2:10 ,
                           clustering.function = clustergram.kmeans,
                           clustergram.plot = clustergram.plot.matlines,
                           FuzzyIndice.plot = FuzzyIndice.plot.matlines,
                           line.width = .004, add.center.points = TRUE, ...)
{
  n <- dim(Data)[1]
  PCA.1 <- Data %*% stats::princomp(Data)$loadings[,1]
  COL <- colorspace::heat_hcl(n)[order(PCA.1)]
  line.width <- rep(line.width, n)
  Y <- NULL
  X <- NULL
  Z <- data.frame(matrix(NA, nrow = 9, ncol = 3,
                        dimnames = list(c(1:9), c("k","PCN","PEN"))))
  centers.points <- list()
  for(k in k.range)
  {
    k.clusters <- clustering.function(Data, k, ...)
    clusters.vec <- k.clusters$cluster
    the.centers <- k.clusters$centers
    Fuzzy.Indices <- k.clusters$Indices
    noise <- unlist(tapply(line.width, clusters.vec, cumsum))[order(seq_along(clusters.vec)[order(clusters.vec)])]
    y <- the.centers[clusters.vec] + noise
    Y <- cbind(Y, y)
    x <- rep(k, length(y))
    X <- cbind(X, x)
    Z[k,1] <- k
    Z[k,2] <- Fuzzy.Indices[2,2]
    Z[k,3] <- Fuzzy.Indices[4,2]
    centers.points[[k]] <- data.frame(y = the.centers , x = rep(k , k))
  }
  x.range <- range(k.range)
  y.range <- range(PCA.1)
  z.range <- range(Z[,2:3], na.rm = TRUE)
  clustergram.plot(X,Y, k.range,
                   x.range, y.range , COL,
                   add.center.points , centers.points)
  FuzzyIndice.plot(Z, k.range,
                   x.range, z.range)
}

