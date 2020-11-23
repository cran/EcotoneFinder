############################# Internal data structure exploration ##################################

#' Tools for internal data structure exploration
#'
#' @param data A community or environmental matrix, containing
#'   species or variables as columns and sites as rows.
#' @param distance.method The distance method to be used for the
#'   calculation of the distance matrix. Must be one of
#'   philentropy::distance
#' @param transpose Logical. If TRUE, the distance matrix is
#'   calculated between species. If FALSE, it is calculated between
#'   sites.
#' @param symm Logical indicating if x should be treated symmetrically;
#'   can only be true when x is a square matrix. See stats::heatmap.
#' @param plot The kind of plot produced by the function.
#'   Can be “heatmap” or “network”.
#' @param palette The colour palette for the network, if spinglass = TRUE.
#'   Must be one of the palettes supported by qgraph.
#'   Default to “colorblind”.
#' @param spinglass Logical. Whether or not to run a spinglass
#'   algorithm to produce statistical groups for the network.
#'   The spinglass algorithm is performed with the CommunityNetwork
#'   function.
#' @param run Number of runs for the spinglass algorithm. Higher numbers
#'   produce more trustable results but rapidly increase computation time.
#'   Default to 10.
#' @param spinglass.groups If spinglass = TRUE, the type of
#'   grouping to use from the results of the spinglass algorithm.
#'   See Details.
#' @param manual.groups If spinglass = FALSE, an object that indicates
#'   which nodes belong together. Can be a list in which each element
#'   is a vector of integers identifying the numbers of the nodes that
#'   belong together, or a factor.
#' @param return.network Logical. If TRUE, the qgraph object is returned
#'   as output of the function.
#' @param ... Additionnal parameters for heatmap or qgraph.
#'
#' @return A plot corresponding to the plot argument.
#'
#' @export
#'
#' @examples
#'  ### Artificial data:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                  b=c(0,250,500),
#'                                                  c=rep(0.01,3)),
#'                                  pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  ## Network of species, with raw spinglass groups:
#'  DistEco(SyntheticTrial[,2:ncol(SyntheticTrial)], transpose = TRUE,
#'          plot = c("network"), spinglass = TRUE, run = 10,
#'          spinglass.groups = c("raw"))
#'
#'  ## Heatmap of species:
#'  DistEco(SyntheticTrial[,2:ncol(SyntheticTrial)], transpose = TRUE,
#'          symm = FALSE, plot = c("heatmap"))
#'
#'
DistEco <- function(data, distance.method = "inner_product", transpose = TRUE, symm = FALSE,
                    plot = c("heatmap", "network"), palette = "colorblind",
                    spinglass = TRUE, run = 10, spinglass.groups = c("rounded", "raw"),
                    manual.groups = NULL, return.network = TRUE, ...)
{
  if (length(plot) > 1 || !(plot %in% c("heatmap", "network"))) {
    stop("plot should be either `heatmap` or `network`")
  }
  if (transpose == TRUE) {
    Dist <- philentropy::distance(t(data), method = distance.method)
    colnames(Dist) <- colnames(data)
    rownames(Dist) <- colnames(data)
  } else {
    Dist <- philentropy::distance(data, method = distance.method)
    colnames(Dist) <- rownames(data)
    rownames(Dist) <- rownames(data)
  }
  if (plot == "network") {
    if (spinglass == TRUE) {
      Qgraph=qgraph::qgraph(stats::as.dist(Dist), layout = "spring",
                            groups = NULL, DoNotPlot = TRUE, ...)
      SpinglassCommunity <- NetworkCommunity(Qgraph, run = run)
      if (spinglass.groups == "raw") {
        network.group <- as.factor(SpinglassCommunity$Memberships$Mean)
      }
      if (spinglass.groups == "rounded") {
        network.group <- as.factor(SpinglassCommunity$Memberships$RoundedMean)
      }
      Qgraph <- qgraph::qgraph(stats::as.dist(Dist), layout = "spring",
                            groups = network.group,
                            palette = palette, ...)

    } else {
      Qgraph <- qgraph::qgraph(stats::as.dist(Dist), layout = "spring",
                            groups = manual.groups, ...)
    }
    if (return.network == TRUE) {
      return(Qgraph)
    }
  }
  if (plot == "heatmap") {
    stats::heatmap(Dist, symm = symm, ...)
  }
}
