################################## Spinglass analysis and network communities ##############################

#' Perform Spinglass algorythm and find networks communities
#'
#' @param networkeco A network object (either qgraph or igraph) or a list
#'   created by the NetworkEco or NetworkEcoSeries functions.
#' @param run Number of runs for the spinglass algorithm. Computation may be
#'   heavy for high numbers.
#'
#' @details The function perform spinglass algorithm on the provided network.
#'   (see spinglass.community() function of the igraph package for more details)
#'   The provided graph is internally transformed into a igraph object if
#'   needed. The function returns a number of summary statistics from the n runs
#'   of the spinglass algorithm. Each run of the spinglass algorithm is done
#'   with a different seed, to ensure different outputs. The seeds are internally
#'   recycled by the \emph{with_seed} fuction of the \emph{withr} package, so
#'   that the global environment is not modified.
#'   The frequencies at which a number of communities are recognised in the
#'   network and the average assignements (rounded or not) of the nodes into
#'   these communities are returned by the function. The latter can help
#'   to statistically define groups for network graphical representations.
#'
#' @return A list containing the number of runs, the number of possible
#'   communities defined by the spinglass algorithm (with frequencies) and the
#'   mean and rounded mean of the assignement of the nodes of the network to
#'   these communities.
#'
#' @export
#'
#' @examples
#'  #### Artificial data:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.01,3)),
#'                                  pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  # Building first network:
#'  Network <- DistEco(SyntheticTrial[,2:ncol(SyntheticTrial)],
#'                     transpose = TRUE, plot = c("network"), spinglass = FALSE,
#'                     return.network = TRUE)
#'
#'  ### Spinglass algorithm (increase number of run for better accuracy):
#'  SpinglassTrial <- NetworkCommunity(Network, run = 5)
#'
#'  ### Network with spinglass groups:
#'  DistEco(SyntheticTrial[,2:ncol(SyntheticTrial)], transpose = TRUE,
#'          plot = c("network"), spinglass = FALSE, return.network = FALSE,
#'          manual.groups = as.factor(SpinglassTrial$Memberships$RoundedMean))
#'
#'
NetworkCommunity <- function (networkeco, run = 100)
{
  returnlist <- list()
  returnlist[["RunNumber"]] <- run
  if (class(networkeco) == "qgraph" || class(networkeco) == "igraph") {
    Network <- networkeco
  } else { if (any(lapply(networkeco, function (x) class(x)) == "qgraph") ||
               any(lapply(networkeco, function (x) class(x)) == "igraph")) {
    if (any(lapply(networkeco, function (x) class(x)) == "qgraph") == TRUE) {
      Network <- networkeco[[as.numeric(which(lapply(networkeco, function (x) class(x)) == "qgraph"))]]
    } else {
      Network <- networkeco[[as.numeric(which(lapply(networkeco, function (x) class(x)) == "igraph"))]]
    }
  }  else {
    stop("networkeco must be (or contain) a qgraph or igraph object")
  }
  }
  if (class(Network) != "igraph") {
    Network <- igraph::as.igraph(Network, attributes = TRUE)
  }
  if (is.null(networkeco$Count) && is.null(networkeco$Distance) && is.null(networkeco$PercentDist)) {
    Node <- as.data.frame(Network[1])
  } else { if (is.null(networkeco$Count) && is.null(networkeco$Distance)) {
    Node <- networkeco$PercentDist
  } else { if (is.null(networkeco$Count)) {
    Node <- networkeco$Distance
  } else { Node <- networkeco$Count }
  }
  }
  SpinglassCommunities <- data.frame(matrix(NA, nrow = run, ncol = nrow(Node)+1,
                                         dimnames = list(1:run, c(rownames(Node), "MaxNumber"))))
  for (i in 1:run) {
    withr::with_seed(i,
    spinglassTest <- igraph::spinglass.community(Network)
    )
    SpinglassCommunities[i,1:(ncol(SpinglassCommunities)-1)] <- spinglassTest$membership
    SpinglassCommunities[i,ncol(SpinglassCommunities)] <- max(spinglassTest$membership)
  }
  for (i in 1:nrow(SpinglassCommunities)) {
    comvec <- as.numeric(SpinglassCommunities[i,-ncol(SpinglassCommunities)])
    comvec <- as.factor(comvec)
    ord <-  as.list(unique(comvec))
    names(ord) <- as.character(1:length(ord))
    levels(comvec) <- ord
    SpinglassCommunities[i,-ncol(SpinglassCommunities)] <- as.numeric(comvec)
  }
  returnlist[["Spinglass"]] <-  SpinglassCommunities
  FreqCommunity <- list()
  FreqCommunity[["ID"]] <- c(min(SpinglassCommunities$MaxNumber):max(SpinglassCommunities$MaxNumber))
  for (i in min(SpinglassCommunities$MaxNumber):max(SpinglassCommunities$MaxNumber)) {
    FreqCommunity[["Count"]][i] <- nrow(SpinglassCommunities[SpinglassCommunities$MaxNumber == i,])
    FreqCommunity[["Proportions"]][i] <- nrow(SpinglassCommunities[SpinglassCommunities$MaxNumber == i,])/nrow(SpinglassCommunities)
  }
  FreqCommunity[["Count"]] <- stats::na.omit(FreqCommunity[["Count"]])
  attributes(FreqCommunity[["Count"]])$na.action <- NULL
  FreqCommunity[["Proportions"]] <- stats::na.omit(FreqCommunity[["Proportions"]])
  attributes(FreqCommunity[["Proportions"]])$na.action <- NULL
  returnlist[["MaxCommunityNumber"]] <-  FreqCommunity
  Mean <- colMeans(as.matrix(SpinglassCommunities[,-ncol(SpinglassCommunities)]))
  RoundedMean <- round(Mean)
  returnlist[["Memberships"]][["Mean"]] <- Mean
  returnlist[["Memberships"]][["RoundedMean"]] <- RoundedMean
  return(returnlist)
}
