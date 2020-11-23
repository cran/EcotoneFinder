################################ Extracting centroids ##############################################

#' Visualisation of fuzzy centroids:
#'
#' @param ecotonefinder A list containing elements named in the same way
#'  than EcotoneFinder function outcomes. Must contain “cmeans”, “fanny”
#'  or “vegclust” results.
#' @param method The fuzzy clustering results from which the centroids
#'   will be extracted.
#' @param normalized Method to normalise the centroid values, either
#'   by “species” or “cluster”. If “none”, the centroids are plotted
#'   without transformation. See details.
#' @param position Set the positions of the bars for the barchart.
#'   This is passed down to the geom_bar function of ggplot. Default
#'   is set to “dodge”.
#' @param threshold Threshold for centroid contribution value under
#'   which the species will not be plotted. Can be used to simplify
#'   plots containing many species. See Details.
#' @param plot Logical. Should the plot be displayed. If FALSE,
#'   the centroids matrix is returned without plotting.
#' @param col Colour vector for the plot. Should be of the same
#'   length that the number of fuzzy clusters.
#' @param labels Character vectors of labels for the legend. Must
#'   be of the same length that the number of fuzzy clusters.
#' @param return.plot Logical. Should the GGplot object be
#'   stored internally (e.g. for multi-ploting). Default is
#'   TRUE.
#' @param main Main title for the plot. See plot.
#' @param xlab A title for the x-axis. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param cex.x cex for the x-axis labels.
#'
#' @details This function extracts and plots the fuzzy centroids species
#'  contributions, according to user-defined normalisation steps and
#'  threshold value. The contributions of the different species in the
#'  fuzzy centroids may be used as a proxy for community compositions.
#'  The cmeans function (cmeans package) and vegclust function (vegclust
#'  package) internally compute the centroid compositions and their
#'  outputs are directly used by the ExtractCentroid function.
#'  The fanny function (cluster package), however, does not provide
#'  internal centroids calculation. They are computed here as:
#'  \deqn{Centroid[cluster j] = \sum[ij] (Membership[ij] x Observation[ij]) / \sum[j] Membership[j]
#'  }
#'  Where the centroid of a cluster is the mean of all observations,
#'  weighted by their degree of belonging to the cluster.
#'  The obtained species contributions to the centroids of the fuzzy
#'  clusters can then be plotted as they are, if normalised = “none”.
#'  To obtain more intuitive units for the interpretation of the
#'  species contributions, two normalisation methods are proposed.
#'  If normalised = “cluster”, the species contributions are given
#'  in percent per clusters (i.e. the sum of all species contributions
#'  in each cluster centroid equals 100). If normalised = “species”,
#'  each species has its contributions summed to 100 (i.e. each species
#'  is in percent per cluster).
#'  For normalised = “none” and normalised = “cluster”, a threshold
#'  value can be specified. Species that do not score above this
#'  threshold will not be displayed on the resulting plot. This can
#'  be used to simplify the outputs, for dataset containing large
#'  number of species.
#'
#' @return A matrix containing the cluster centroids.
#'
#' @export
#'
#' @examples
#' ##### Artificial dataset & analyses:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 20, CommunityNum = 4,
#'                                  SpCo = NULL ,Length = 500,
#'                                  Parameters=list(a = rep(60, 4),
#'                                                  b = c(0,150,350,500),
#'                                                  c = rep(0.015,4)),
#'                                  dev.c = 0.007,
#'                                  pal = c("#008585", "#B8CDAE", "#E6C186", "#C7522B"))
#'
#'  EcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
#'                             dist = SyntheticTrial$Distance,
#'                             method = "all", groups=4,
#'                             standardize = "hellinger",
#'                             diversity="all")
#'
#'  ##### Centroid plot without normalisation:
#'  Centroid <- ExtractCentroid(EcoFinder, method = "fanny",
#'                              normalized = "none",  threshold = 0,
#'                              plot = TRUE, position = "dodge",
#'                              col = colorspace::heat_hcl(4))
#'
#'  ##### Centroid plot normalised by clusters:
#'  Centroid <- ExtractCentroid(EcoFinder, method = "fanny",
#'                              normalized = "cluster",  threshold = 0,
#'                              plot = TRUE, position = "dodge",
#'                              col = colorspace::heat_hcl(4))
#'
ExtractCentroid <- function (ecotonefinder, method = c("fanny", "cmeans", "vegclust"), normalized = c("species", "cluster", "none"),
                             position = "dodge", threshold = 0, plot = TRUE, col = NULL, return.plot = TRUE, labels = ggplot2::waiver(),
                             main = "Community composition", xlab = "species", ylab = "Centroid contribution",
                             cex.x = 12)
{
  if (any(is.na(names(ecotonefinder[method])))) {
    stop("ecotonefinder must comntains the elements required in 'method'.")
  }
  returnList <- list()
  if (method == "cmeans") {
    if (normalized == "species") {
      GGData <- as.data.frame(t(vegan::decostand(t(ecotonefinder$cmeans$centers), method = "total")*100))
    }
    if (normalized == "cluster") {
      GGData <- as.data.frame(vegan::decostand(ecotonefinder$cmeans$centers, method = "total")*100)
      GGData <- GGData[,as.numeric(which(base::colSums(GGData) > threshold))]
    }
    if (normalized == "none") {
      GGData <- as.data.frame(ecotonefinder$cmeans$centers)
      GGData <- GGData[,as.numeric(which(base::colSums(GGData) > threshold))]
    }
    GGData$Cluster <- base::rownames(GGData)
    GGData <- reshape::melt(GGData, id.vars = "Cluster")
    colnames(GGData) <- c("Cluster", "Species", "Contribution")
    GGData$Cluster <- as.factor(GGData$Cluster)
    returnList[["centroids"]] <- GGData
  }
  if (method == "vegclust") {
    if (normalized == "species") {
      GGData <- as.data.frame(t(vegan::decostand(t(ecotonefinder$vegclust$mobileCenters), method = "total")*100))
    }
    if (normalized == "cluster") {
      GGData <- as.data.frame(vegan::decostand(ecotonefinder$vegclust$mobileCenters, method = "total")*100)
      GGData <- GGData[,as.numeric(which(base::colSums(GGData) > threshold))]
    }
    if (normalized == "none") {
      GGData <- as.data.frame(ecotonefinder$vegclust$mobileCenters)
      GGData <- GGData[,as.numeric(which(base::colSums(GGData) > threshold))]
    }
    GGData$Cluster <- base::rownames(GGData)
    GGData <- reshape::melt(GGData, id.vars = "Cluster")
    colnames(GGData) <- c("Cluster", "Species", "Contribution")
    GGData$Cluster <- as.factor(GGData$Cluster)
    returnList[["centroids"]] <- GGData
  }
  if (method == "fanny") {
    FannyCenters <- list()
    FannyCenters <- matrix(NA, ecotonefinder$groups,
                           ncol(ecotonefinder$data))
    for (j in 1:ecotonefinder$groups){
      FannyCenters[j,] <- colSums(ecotonefinder$fanny$membership[,j] * ecotonefinder$data) / sum(ecotonefinder$fanny$membership[,j])
    }
    colnames(FannyCenters) <- colnames(ecotonefinder$data)
    if (normalized == "species") {
      GGData <- as.data.frame(t(vegan::decostand(t(FannyCenters), method = "total")*100))
    }
    if (normalized == "cluster") {
      GGData <- as.data.frame(vegan::decostand(FannyCenters, method = "total")*100)
      GGData <- GGData[,as.numeric(which(base::colSums(GGData) > threshold))]
    }
    if (normalized == "none") {
      GGData <- as.data.frame(FannyCenters)
      GGData <- GGData[,as.numeric(which(base::colSums(GGData) > threshold))]
    }
    GGData$Cluster <- base::rownames(GGData)
    GGData <- reshape::melt(GGData, id.vars = "Cluster")
    colnames(GGData) <- c("Cluster", "Species", "Contribution")
    GGData$Cluster <- as.factor(GGData$Cluster)
    returnList[["centroids"]] <- GGData
  }
  if (is.null(col)) {
    col <- colorspace::heat_hcl(length(levels(GGData$Cluster)))
  }
  GG <- ggplot2::ggplot(GGData) +
    ggplot2::geom_bar(ggplot2::aes(x = GGData[["Species"]], y = GGData[["Contribution"]], group = GGData[["Cluster"]], fill = GGData[["Cluster"]]),
                      position = position, stat = "summary", fun = mean) +
    ggplot2::scale_fill_manual(values=col, labels = labels) +
    ggplot2::ggtitle(main) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = cex.x),
          plot.title = ggplot2::element_text(hjust = 0.5))
  if (plot == TRUE) {
    print(GG)
  }
  if (return.plot == TRUE) {
    returnList[["GGplot"]] <- GG
  }
  return(returnList)
}
