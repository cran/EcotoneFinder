########################################## Networks and Distance analyses #######################################

###################################### For unique set of data ##########################################

#' Networks for ecotones and communities
#'
#' @param ecotonefinder A list containing elements named in the same way than
#'   EcotoneFinder function outcomes. Must contain cmeans results or vegclust
#'   results.
#' @param threshold If count = T, the membership grade threshold used to sort
#'   the species in the different clusters.
#' @param plot.type Which graphical representation to be plotted. Among
#'   "percentage", "corrplot", "heatmap","network".
#' @param method The membership computation method to be used. One of "cmeans"
#'   or "vegclust". Must be present in the ecotonefinder list.
#' @param dist.method Distance method for the computation of a distance matrix,
#'   when dist = "raw" and dist = "relative".
#' @param plot If plot = "species", the distances are computed between the
#'   species in the data. If plot = "community", the distances are computed
#'   between the cluster centroids.
#' @param order.sp Vector providing the order in which to arrange the species.
#'   If NULL, the column order will be kept.
#' @param dist The type of data on which distance calculations are made from. If
#'   dist = "raw", the distance matrix is computed from the membership matrix
#'   directly. if dist = "relative", the distance matrix is computed from the
#'   relative memberships grades of each species in the clusters (between 0 and
#'   1). If dist = "count", the species are assigned to clusters according to
#'   the threshold and the distance matrix is computed from the number of common
#'   species between the different clusters. See details.
#' @param no.plot Logical. Should the plot be displayed?. Set to TRUE to gain
#'   computation time with large community matrix.
#' @param network.group Grouping parameter for the network. Can be user defined
#'   (see qgraph documentation for details) but must be a factor of the same
#'   lenght as the nodes of the graph.
#' @param ... Additional arguments to be passed to the plotting functions, see
#'   details.
#'
#' @details The NetworkEco function provides wqys to explore the relations
#'   between fuzzy clusters. Several options are implemented. If dist = "raw",
#'   it computes a distance matrix from the membership grade matrix directly. If
#'   dist = "relative", the membership grades are standardized so that the sum
#'   of the membership grades of a given species equals to 1 for every points
#'   along the gradient (which corresponds to a percentage ot membership in each
#'   cluster). If dist = "count", the standardized membership grades of the are
#'   used to assign species in the community to a unique cluster and the number
#'   of common species between pairs of clusters is counted. The assignement of
#'   species to clusters is done by listing all the species that score a
#'   membership grade higher than the specified threshold in a cluster. The
#'   resulting list of species are then compared to one another.
#'
#'   The function also allows the computation of distances between species
#'   rather than between clusters, when plot = "species". This can only be done
#'   from the memberships grades (raw or relative) and this argument will be
#'   disregarded if dist = "count".
#'
#'   Several methods of visualisation are implemented: "percentage", "corrplot",
#'   "heatmap" and "network". If "percentage", a barplot (using ggplot2) of the
#'   standardized memberships grades per fuzzy cluster is plotted. It always
#'   plot the standardized membership grades regardless of the chosen dist
#'   option, but if dist = "count" or dist = "raw" are chosen, the function
#'   still compute the corresponding distance matrices and return them ti the
#'   output list. For time efficiency, it is not recommended to plot it when the
#'   number of species in the community is large (>100). "corrplot" and
#'   "heatmap" produce correlation matrix and heat map. The "network" is based
#'   on the qgraph function of the qgraph package. The ... argument may be used
#'   to pass additional arguments to the plotting functions (for graphical
#'   purposes).
#'
#' @return A list containing the percentage matrix, the distance matrix and the
#'   network object (depending of the arguments passed to the function)
#'
#' @export
#'
#' @examples
#'  #### Artificial dataset:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.015,3)),
#'                                  pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  ## Analyses:
#'  SyntheticEcoFinder <- EcotoneFinder(data = SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all",
#'                                      groups = 3, standardize = "hellinger",
#'                                      diversity = "all")
#'
#'  ## Percentage plot:
#'  SyntheticNetwork <- NetworkEco(SyntheticEcoFinder, threshold = .3, method = "cmeans",
#'                                 plot.type = "percentage", dist = "count")
#'
#'  ## Heatmap plot:
#'  SyntheticNetwork <- NetworkEco(SyntheticEcoFinder, plot.type = "heatmap",
#'                                 method = "cmeans", dist = "raw", plot = "species")
#'
#'  ## Network:
#'  # From raw membership grades:
#'  SyntheticNetwork <- NetworkEco(SyntheticEcoFinder, plot.type = "network",
#'                                 method = "cmeans", dist = "raw", plot = "species")
#'
#'  # From number of species per clusters:
#'  SyntheticNetwork <- NetworkEco(SyntheticEcoFinder, plot.type = "network", threshold = .3,
#'                                 method = "cmeans", dist = "count", plot = "community",
#'                                 layout = "spring")
#'
NetworkEco <- function (ecotonefinder, threshold = 0.80, plot.type = c("percentage", "corrplot", "heatmap","network"),
                        method = c("cmeans", "vegclust"), dist.method = "inner_product", plot = c("species","community"),
                        order.sp = NULL, dist = c("count", "relative", "raw"), no.plot = FALSE, network.group = NULL, ...)
{
  returnList <- list()
  returnList[["threshold"]] <- threshold
  if(is.null(ecotonefinder$cmeans) && is.null(ecotonefinder$vegclust)) {
    stop("NetworkEcoSeries require cmeans or vegclust analyses in the input list")
  }
  if (length(method) > 1 || (method != "cmeans" && method != "vegclust")) {
    stop("method should be either cmeans or vegclust")
  }
  if (method == "cmeans") {
    x <- "cmeans"
    y <- "centers"
  }
  if (method == "vegclust") {
    x <- "vegclust"
    y <- "mobileCenters"
  }
  if (any(dist == c("count","relative")) || any(plot.type == "percentage")) {
    PercentMemb <- list()
    for (i in 1:ncol(ecotonefinder[[x]][[y]])) {
      PercentMemb[[i]] <- t(ecotonefinder[[x]][[y]])[i,]/sum(t(ecotonefinder[[x]][[y]])[i,])
    }
    PercentMemb <- purrr::reduce(PercentMemb, rbind)
    rownames(PercentMemb) <- colnames(ecotonefinder$cmeans$centers)
    returnList[["percentage"]] <- PercentMemb
  }
  if (any(plot.type == "percentage") && no.plot == FALSE) {
    MPercentMemb <- reshape::melt(PercentMemb)
    if (is.null(order.sp) == TRUE) {
      MPercentMemb[,1] <- ordered(MPercentMemb[,1], levels = colnames(ecotonefinder$cmeans$centers))
    } else { MPercentMemb[,1] <- ordered(MPercentMemb[,1], levels = order.sp) }
    Plot <-  ggplot2::ggplot(data=MPercentMemb) +
      ggplot2::geom_bar(ggplot2::aes(x = MPercentMemb[,1], y = MPercentMemb[["value"]], group = MPercentMemb[,1], fill = MPercentMemb[,1]), stat = "identity") +
      ggplot2::facet_grid(~MPercentMemb[,2]) +
      ggplot2::ylab("Normalised membership grades")
    print(Plot)
  }
  if (dist == "relative") {
    if (plot == "community") {
      PercentDist <- philentropy::distance(t(PercentMemb), method = dist.method)
      colnames(PercentDist) <- colnames(PercentMemb)
      rownames(PercentDist) <- colnames(PercentMemb)
    }
    if (plot == "species") {
      PercentDist <- philentropy::distance(PercentMemb, method = dist.method)
      colnames(PercentDist) <- rownames(PercentMemb)
      rownames(PercentDist) <- rownames(PercentMemb)
    }
    returnList[["PercentDist"]] <- PercentDist
  }
  if (dist == "count") {
    MList <- list()
    for (i in 1:ncol(PercentMemb)) {
      MList[[i]] <- as.character(rownames(PercentMemb[PercentMemb[,i] >= threshold,]))
    }
    CountDist <- data.frame(matrix(NA, nrow = ncol(PercentMemb), ncol = ncol(PercentMemb),
                                dimnames = list(colnames(PercentMemb), colnames(PercentMemb))))
    for (i in 1:ncol(PercentMemb)) {
      for (j in 1:ncol(PercentMemb)) {
        CountDist[i,j] <- length(intersect(grep(paste(MList[[i]], collapse = "|"), MList[[j]], value = TRUE),
                                         grep(paste(MList[[j]], collapse = "|"), MList[[i]], value = TRUE)))
      }
    }
    returnList[["Count"]] <- CountDist
  }
  if (dist == "raw") {
    if (plot == "community") {
      Dist <- philentropy::distance(ecotonefinder[[x]][[y]], method = dist.method)
      colnames(Dist) <- rownames(ecotonefinder[[x]][[y]])
      rownames(Dist) <- rownames(ecotonefinder[[x]][[y]])
    }
    if (plot == "species") {
      Dist <- philentropy::distance(t(ecotonefinder[[x]][[y]]), method = dist.method)
      colnames(Dist) <- colnames(ecotonefinder[[x]][[y]])
      rownames(Dist) <- colnames(ecotonefinder[[x]][[y]])
    }
    returnList[["Distance"]] <- Dist
  }
  if (any(plot.type == "heatmap")) {
    if (dist == "raw") {
      stats::heatmap(Dist, ...)
    }
    if (dist == "count") {
      stats::heatmap(as.matrix(CountDist), ...)
    }
    if (dist == "relative") {
      stats::heatmap(as.matrix(PercentDist), ...)
    }
  }
  if (any(plot.type == "corrplot")) {
    if (dist == "raw") {
      corrplot::corrplot(as.matrix(Dist), ...)
    }
    if (dist == "count") {
      corrplot::corrplot(as.matrix(CountDist), ...)
    }
    if (dist == "relative") {
      corrplot::corrplot(as.matrix(PercentDist), ...)
    }
  }
  if (any(plot.type == "network")) {
    if (dist == "count") {
      Qgraph <- qgraph::qgraph(stats::as.dist(CountDist), groups=network.group, DoNotPlot=no.plot, ...)
    }
    if (dist == "raw") {
      Qgraph <- qgraph::qgraph(Dist, groups=network.group, DoNotPlot=no.plot, ...)
    }
    if (dist == "relative") {
      Qgraph <- qgraph::qgraph(PercentDist, groups=network.group, DoNotPlot=no.plot, ...)
    }
    returnList[["Network"]] <- Qgraph
  }
  return(returnList)
}

################################### For data series ####################################################

#' Networkeco for data series
#'
#' @param ecotonefinder A list containing elements named in the same way than
#'   EcotoneFinderSeries function outcomes. Must contain cmeans results or
#'   vegclust results.
#' @param threshold If dist = "count", the membership grade threshold used to
#'   sort the species in the different clusters.
#' @param method The membership computation method to be used. Must be present
#'   in the ecotonefinder list.
#' @param plot.type Which graphical representation to be plotted. Among
#'   "percentage", "corrplot", "heatmap","network"
#' @param plot If plot = "species", the distances are computed between the
#'   species in the data. If plot = "community", the distances are computed
#'   between the cluster centroids.
#' @param no.plot Logical. Should the plot be displayed?. Set to TRUE to gain
#'   computation time with large community matrix.
#' @param order.sp Vector providing the order in which to arrange the species.
#'   If NULL, the column order will be kept.
#' @param dist.method Distance method for the computation of a distance matrix,
#'   when dist = "raw" or dist = "percent".
#' @param dist The type of data on which distance calculations are made from. If
#'   dist = "raw", the distance matrix is computed from the membership matrix
#'   directly. if dist = "relative", the distance matrix is computed from the
#'   relative memberships grades of each species in the clusters (between 0 and
#'   1). If dist = "count", the species are assigned to clusters according to
#'   the threshold and the distance matrix is computed from the number of common
#'   species between the different clusters. See details.
#' @param network.group If network.group = "site" the nodes of the networks will
#'   be colored according to the different times or sites of the series. If
#'   network.group = "cluster" the nodes of the network will be colored
#'   according to the different fuzzy clusters. Can be user defined (see qgraph
#'   documentation for details) but must be a factor of the same lenght as the
#'   nodes of the graph.
#' @param method.corr If plot.type = "corrplot", the method to be used for the
#'   corrplot. Must be one of "circle", "square", "ellipse", "number", "shade",
#'   "color", "pie". Default to "number".
#' @param ... Additional arguments to be passed to the plotting functions, see
#'   details.
#'
#' @details NetworkEcoSeries is a generalisation of the NetworkEco function to
#'   analyses space/time series. The ... argument may be used to pass additional
#'   arguments to the plotting functions (for graphical purposes).
#'
#' @return A list containing the percentage matrix, the distance matrix and the
#'   network object (depending of the arguments passed to the function)
#'
#' @export
#'
#' @examples
#'   \donttest{
#'   SyntheticTrialSeries <- SyntheticDataSeries(CommunityPool = 40,
#'                                               CommunityNum = 4, SpCo = NULL,
#'                                               Length = 500, SeriesNum = 5,
#'                                               Parameters = list(a=rep(60, 4),
#'                                                               b=c(0,200,350,500),
#'                                                               c=rep(0.03,4)),
#'                                               pal = c("#008585", "#B8CDAE", "#E6C186", "#C7522B"),
#'                                               replacement = TRUE,
#'                                               Parameters.repl = TRUE)
#'
#'   EcoTimeSeriesTrial <- EcotoneFinderSeries(data = SyntheticTrialSeries,
#'                                             dist = "Distance",
#'                                             method = c("cmeans","vegclust"),
#'                                             series = "Time", groups = 4,
#'                                             standardize = "hellinger", na.rm=TRUE)
#'
#' #### Network from the common number of species above membership threshold between clusters:
#'   SyntheticNetworkSeries <- NetworkEcoSeries(EcoTimeSeriesTrial, threshold = .2,
#'                                              method = "cmeans", plot.type = "network",
#'                                              plot = "community", dist = "count",
#'                                              network.group = "cluster",
#'                                              dist.method = "inner_product",
#'                                              no.plot = FALSE, layout = "spring",
#'                                              shape = "ellipse",
#'                                              palette = "colorblind")
#'
#' #### Network of relations between species from their raw membership values in each cluster:
#'   SyntheticNetworkSeries <- NetworkEcoSeries(EcoTimeSeriesTrial, threshold = .2,
#'                                              method = "cmeans", plot.type = "network",
#'                                              plot = "species", dist = "raw",
#'                                              dist.method = "inner_product",
#'                                              no.plot = FALSE, layout = "spring",
#'                                              shape = "ellipse",
#'                                              palette = "colorblind")
#'  }
#'
NetworkEcoSeries <- function (ecotonefinder, threshold = 0.80, method = c("cmeans","vegclust"), plot.type = c("percentage", "heatmap", "corrplot", "network"),
                              plot = c("species", "community"), no.plot = FALSE, order.sp = NULL, dist.method = "inner_product",
                              dist = c("count", "relative", "raw"), network.group = c("site","cluster"), method.corr = "number", ...)
{
  returnList <- list()
  returnList[["threshold"]] <- threshold
  CommunityNum <- as.numeric(lapply(ecotonefinder, function (n) n$groups))
  Name <- NULL
  for (k in 1:length(ecotonefinder)) {
    for (j in 1:CommunityNum[k]) {
      Name[paste(k,j)] <- paste(names(ecotonefinder[k]),"Cluster",j)
    }
  }
  if(suppressWarnings(any(lapply(ecotonefinder, function (n) is.null(n$cmeans))) == TRUE &&
                      any(lapply(ecotonefinder, function (n) is.null(n$vegclust))) == TRUE)) {
    stop("NetworkEcoSeries require cmeans or vegclust analyses in the input list")
  }
  if (length(method) > 1 || (method != "cmeans" && method != "vegclust")) {
    stop("method should be either cmeans or vegclust")
  }
  if (method == "cmeans") {
    x="cmeans"
    y="centers"
  }
  if (method == "vegclust") {
    x="vegclust"
    y="mobileCenters"
  }
  if (length(unique(as.numeric(lapply(lapply(ecotonefinder, function (n) n[[x]][[y]]), length)))) != 1) {
    Memb <- lapply(ecotonefinder, function (n) n[[x]][[y]])
    SpeList <- unique(as.character(unlist(lapply(Memb, colnames))))
    for (k in 1:length(Memb)) {
      for (i in 1:length(SpeList)) {
        if (is.null(Memb[[k]][[SpeList[i]]])) {
          Memb[[k]][[SpeList[i]]] <- 0
        }
      }
    }
  } else {
    Memb <- lapply(ecotonefinder, function (n) n[[x]][[y]])
  }
  if (any(dist == c("count","relative")) || any(plot.type == "percentage")) {
    PercentMemb <- list()
    for (k in 1:length(ecotonefinder)) {
      for (i in 1:ncol(Memb[[k]])) {
        PercentMemb[[names(ecotonefinder)[k]]][[i]] <- t(Memb[[k]])[i,]/sum(t(Memb[[k]])[i,])
      }
      PercentMemb[[names(ecotonefinder)[k]]] <- purrr::reduce(PercentMemb[[names(ecotonefinder)[k]]], rbind)
      rownames(PercentMemb[[names(ecotonefinder)[k]]]) <- colnames(Memb[[k]])
      PercentMemb[[names(ecotonefinder)[k]]] <- replace(PercentMemb[[names(ecotonefinder)[k]]], is.nan(PercentMemb[[names(ecotonefinder)[k]]]), 0)
    }
    returnList[["percentage"]] <- PercentMemb
  }
  if (any(plot.type == "percentage") && no.plot == FALSE) {
    for (k in 1:length(ecotonefinder)) {
      MPercentMemb <- reshape::melt(PercentMemb[[names(ecotonefinder)[k]]])
      if (is.null(order.sp) == TRUE) {
        MPercentMemb[,1] <- ordered(MPercentMemb[,1], levels = colnames(ecotonefinder[[names(ecotonefinder)[k]]][[x]][[y]]))
      } else { MPercentMemb[,1] <- ordered(MPercentMemb[,1], levels = order.sp) }
      Plot = ggplot2::ggplot(data = MPercentMemb) +
        ggplot2::geom_bar(ggplot2::aes(x=MPercentMemb[,1], y = MPercentMemb[["value"]], group = MPercentMemb[,1], fill = MPercentMemb[,1]), stat = "identity") +
        ggplot2::facet_grid(~MPercentMemb[,2]) +
        ggplot2::ggtitle(paste(names(ecotonefinder)[k])) +
        ggplot2::ylab("Normalised membership grades")
      print(Plot)
    }
  }
  if (dist == "count") {
    MList <- list()
    for (k in 1:length(ecotonefinder)) {
      for (i in 1:ncol(PercentMemb[[names(ecotonefinder)[k]]])) {
        MList[[names(ecotonefinder)[k]]][[i]] <- as.character(rownames(PercentMemb[[names(ecotonefinder)[k]]][PercentMemb[[names(ecotonefinder)[k]]][,i] >= threshold,]))
      }
    }
    CallList <- do.call(EcotoneFinder::rbindna, MList)
    CountDist <- data.frame(matrix(NA, nrow = length(CallList), ncol = length(CallList)))
    for (i in 1:length(CallList)) {
      if (is.na(CallList[i]) == FALSE) {
        for (j in 1:length(CallList)) {
          if (is.na(CallList[j]) == FALSE) {
            CountDist[i,j] <- length(intersect(grep(paste(CallList[[i]], collapse = "|"), CallList[[j]], value = TRUE),
                                              grep(paste(CallList[[j]], collapse = "|"), CallList[[i]], value = TRUE)))
          }
        }
      }
    }
    CountDist <- CountDist[!apply(is.na(CountDist), 1, all),!apply(is.na(CountDist), 1, all)]
    rownames(CountDist) <- Name
    colnames(CountDist) <- Name
    returnList[["Count"]] <- CountDist
  }
  if (dist == "relative") {
    if (plot == "community") {
      PercentDist <- do.call(rbind, as.data.frame(PercentMemb))
      PercentDist <- philentropy::distance(PercentDist, method=dist.method)
      rownames(PercentDist) <- Name
      colnames(PercentDist) <- Name
    }
    if (plot == "species") {
      PercentDist <- do.call(cbind, as.data.frame(PercentMemb))
      PercentDist <- philentropy::distance(PercentDist, method = dist.method)
      rownames(PercentDist) <- rownames(as.data.frame(PercentMemb))
      colnames(PercentDist) <- rownames(as.data.frame(PercentMemb))
    }
    returnList[["PercentDist"]] <- PercentDist
  }
  if (dist == "raw")  {
    if (plot == "community") {
      Dist <- do.call(rbind, rbind(Memb))
      Dist <- philentropy::distance(Dist, method = dist.method)
      rownames(Dist) <- Name
      colnames(Dist) <- Name
    }
    if (plot == "species") {
      Dist <- t(do.call(rbind, rbind(Memb)))
      Name <- rownames(Dist)
      Dist <- philentropy::distance(Dist, method=dist.method)
      rownames(Dist) <- Name
      colnames(Dist) <- Name
    }
    returnList[["Distance"]] <- Dist
  }
  if (plot.type == "heatmap" && no.plot == FALSE) {
    if (dist == "raw") {
      stats::heatmap(as.matrix(Dist), ...)
    }
    if (dist == "count") {
      stats::heatmap(as.matrix(CountDist), ...)
    }
    if (dist == "relative") {
      stats::heatmap(as.matrix(PercentDist), ...)
    }
  }
  if (plot.type == "corrplot" && no.plot == FALSE) {
    if (dist == "raw") {
      corrplot::corrplot(as.matrix(Dist), method = method.corr, ...)
    }
    if (dist == "count") {
      corrplot::corrplot(as.matrix(CountDist), method = method.corr, ...)
    }
    if (dist == "relative") {
      corrplot::corrplot(as.matrix(PercentDist), method = method.corr, ...)
    }
  }
  if (plot.type == "network") {
    if (dist == "count") {
      if (any(network.group == "cluster")) {
        for (k in 1:length(ecotonefinder)) {
          for (i in 1:CommunityNum[k]) {
            network.group[grep(paste("Cluster", i), rownames(CountDist), value = FALSE)] <- paste("Cluster", i)
          }
        }
      }
      if (any(network.group == "site")) {
        for (i in 1:length(ecotonefinder)) {
          network.group[grep(names(ecotonefinder)[i], rownames(CountDist), value = FALSE)] <- names(ecotonefinder)[i]
        }
      }
      Qgraph <- qgraph::qgraph(stats::as.dist(CountDist), groups = network.group, DoNotPlot = no.plot, ...)
    }
    if (dist == "raw") {
        if (any(network.group=="cluster")) {
          for (k in 1:length(ecotonefinder)) {
            for (i in 1:CommunityNum[k]) {
              network.group[grep(paste("Cluster", i), rownames(Dist), value = FALSE)] <- paste("Cluster", i)
            }
          }
        }
        if (any(network.group == "site")) {
          for (i in 1:length(ecotonefinder)) {
            network.group[grep(names(ecotonefinder)[i], rownames(Dist), value = FALSE)] <- names(ecotonefinder)[i]
          }
        }
      Qgraph <- qgraph::qgraph(Dist, groups = network.group, DoNotPlot = no.plot, ...)
    }
    if (dist == "relative") {
        if (any(network.group == "cluster")) {
          for (k in 1:length(ecotonefinder)) {
            for (i in 1:CommunityNum[k]) {
              network.group[grep(paste("Cluster", i), rownames(PercentDist), value = FALSE)] <- paste("Cluster", i)
            }
          }
        }
        if (any(network.group == "site")) {
          for (i in 1:length(ecotonefinder)) {
            network.group[grep(names(ecotonefinder)[i], rownames(PercentDist), value = FALSE)] <- names(ecotonefinder)[i]
          }
        }
      Qgraph <- qgraph::qgraph(PercentDist, groups = network.group, DoNotPlot = no.plot, ...)
    }
    returnList[["Network"]] <- Qgraph
  }
  return(returnList)
}
