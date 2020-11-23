##################################### EcotoneFinder Analyses ##############################################

##################################### Unique set of data ################################################

#' Wraper function to perform ecological gradient analysis
#'
#' @param data A dataframe containing species as columns and sites as rows.
#'   May also contain environmental parameters, as long as the parameters are
#'   as columns and the site as rows.
#' @param dist A vector or column containing the gradient along which the
#'   analysis will be done. Must be of the same lenght as data.
#' @param method One of c("dca", "fanny", "vegclust", "diversity", "cmeans", "all").
#'   Tell the function which analysis to perform. See details.
#' @param groups Interger. The desired number of clusters if any of the
#'   clustering method is selected.
#' @param m.exp Integer. The membership exponent for any of the clustering
#'   method.
#' @param standardize Standardize method to apply to the data before further
#'   analysis (for fanny and cmeans). Must be one of decostand methods (see
#'   decostand).
#' @param seed Integer or NULL. Set a seed for initial membership matrix for
#'   cmeans and vegclust algorithms. Recomended for time series, to keep
#'   a more consistent labelling of fuzzy clusters along the gradient.
#'   See Details.
#' @param diversity diversity indice to be calculated. See details.
#' @param na.rm Logical. Should NAs be removed.
#'
#' @details EcotoneFinder is a wraper function to perform multiple ecological
#'   gradient analysis at once. The implemented methods are Detrended
#'   Correspondance Analysis (DCA) - see Brownstein et al. 2013 - Fuzzy C-Means
#'   (FCM) - see DeCaceres et al., 2010 - and the calculation of diversity
#'   indices - see Jost, 2007. The DCA is intenally performed by the decorana
#'   function of the vegan package. The FCM analyses can be performed by the
#'   fanny function (cluster package), the vegclust function (vegclust package)
#'   and the cmeans function (e1071 package) - for comparison purposes as the
#'   outcome of the analyses might differ.
#'   The vegclust and cmeans algorithms use random number generators to create
#'   a matrix of initial centers. Setting a seed with the \emph{seed} argument
#'   guaranty the reproducibility of the outputs. This argument can also be set
#'   to NULL, to preserve the randomness of initial centres.
#'
#'   It is recommended to stadardize data before applying fanny or cmeans
#'   analysis. See decostand documentation (vegan package) for information on
#'   standardisation methods. Must be one of:
#'   "total","max","freq","normalize","range","pa","chi.square","hellinger","log".
#'   If "log" is chosen, the user will be asked to provide the base to be used
#'   upon launching the function.
#'
#'   Several diversity incides have been implemented, as they are supposed to
#'   react to ecological gradients. It includes the Shannon index, the Pielou
#'   eveness and species richness - computed with the diversity function of the
#'   vegan package. expShannon corresponds to the shannon index in terms of
#'   effective number of species (see Jost, 2007). If "all" is selected, all the
#'   implemented indices will be calculated.
#'
#'   If "all" is selected in the method argument, all the implemented methods
#'   will be applied.
#'
#' @return Ecofinder returns a list containing the original data, the value of
#'   the main arguments used in the function and the outcome of the selected
#'   analyses.
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
#'
#'
#'
EcotoneFinder <- function (data, dist, method = c("dca", "fanny", "vegclust", "diversity", "cmeans", "all"),
                           groups = NULL, m.exp = 2, standardize = NULL, seed = 1,
                           diversity = c("shannon", "richness", "expShannon", "pielou", "all"), na.rm = FALSE)
{
  # Checks:
  if (is.null(data) || all(sapply(data, is.numeric)) != TRUE) {
    stop("data must be provided and numerical.")
  }
  if (is.numeric(dist) == FALSE ) {
    stop("dist must be a numeric vector or column.")
  }
  if (any(method == c("fanny", "vegclust", "cmeans", "all")) && is.null(groups)) {
    stop("When methods include fuzzy clustering (i.e. fanny, vegclust or cmeans), groups must be provided.")
  }
  # List:
  list_output <- list()
  list_output[[paste("method")]] <- method
  list_output[[paste("data")]] <- data
  list_output[[paste("distance")]] <- dist
  list_output[[paste("groups")]] <- groups
  # Methods:
  if (any(method == "dca") || any(method == "all")) {
    dca=vegan::decorana(data)
    list_output[[paste("dca")]] <- dca
  }
  if (any(method == "fanny") || any(method == "all")) {
    if (is.null(standardize) == TRUE) {
      Fanny <- cluster::fanny(data, k = groups, memb.exp = m.exp)
    } else {
      if (standardize == "log") {
        log <- readline(prompt = "Enter the logarithm base to be used: ")
        StandardCommunity <- vegan::decostand(data, method = standardize, logbase = log, na.rm = na.rm)
      } else {
        StandardCommunity <- vegan::decostand(data, method = standardize, na.rm = na.rm)
      }
      Fanny <- cluster::fanny(StandardCommunity, k = groups, memb.exp = m.exp)
      list_output[[paste("stadardize")]] <- standardize
    }
    list_output[[paste("fanny")]] <- Fanny
  }
  if (any(method == "vegclust") || any(method == "all")) {
    if (is.null(seed)) {
      Vegclust <- vegclust::vegclust(data, mobileCenters = groups, method = "FCM", m = m.exp)
    } else {
    withr::with_seed(seed,
    Vegclust <- vegclust::vegclust(data, mobileCenters = groups, method = "FCM", m = m.exp)
    )
    }
    list_output[[paste("vegclust")]] <- Vegclust
  }
  if (any(method == "cmeans") || any(method == "all")) {
    if (is.null(standardize) == TRUE) {
      if (is.null(seed)) {
        Cmeans <- e1071::cmeans(data, centers = groups, method = "cmeans", m = m.exp)
      } else {
        withr::with_seed(seed,
                         Cmeans <- e1071::cmeans(data, centers = groups, method = "cmeans", m = m.exp)
        )
      }
    } else {
      if (standardize == "log") {
        log <- readline(prompt = "Enter the logarithm base to be used: ")
        StandardCommunity <- vegan::decostand(data, method = standardize, logbase = log, na.rm = na.rm)
      } else {
        StandardCommunity <- vegan::decostand(data, method = standardize, na.rm = na.rm)
      }
      if (is.null(seed)) {
        Cmeans <- e1071::cmeans(StandardCommunity, centers = groups, method = "cmeans", m = m.exp)
      } else {
        withr::with_seed(seed,
                         Cmeans <- e1071::cmeans(StandardCommunity, centers = groups, method = "cmeans", m = m.exp)
        )
      }
    }
    Cmeans$membership <- as.data.frame(Cmeans$membership)[,unique(Cmeans$cluster)]
    colnames(Cmeans$membership) <- sort(unique(Cmeans$cluster))
    Cmeans$centers <- as.data.frame(Cmeans$centers)[unique(Cmeans$cluster),]
    rownames(Cmeans$centers) <- sort(unique(Cmeans$cluster))
    list_output[[paste("cmeans")]] <- Cmeans
  }
  if (any(method == "diversity") || any(method == "all")) {
    list_diversity <- list()
    if (any(diversity == "richness") || any(diversity == "all")) {
      SpeciesRichness <- vegan::specnumber(data)
      list_diversity[[paste("SpeciesRichness")]] <- SpeciesRichness
    }
    if (any(diversity == "shannon") || any(diversity == "all")) {
      Shannon <- vegan::diversity(data, "shannon")
      list_diversity[[paste("Shannon")]] <- Shannon
    }
    if (any(diversity == "expShannon") || any(diversity == "all")) {
      ExpShannon <- vegan::diversity(data, "shannon")
      ExpShannon <- exp(ExpShannon)
      list_diversity[[paste("ExpShannon")]] <- ExpShannon
    }
    if (any(diversity == "pielou") || any(diversity == "all")) {
      Pielou <- vegan::diversity(data, "shannon")/vegan::specnumber(data)
      list_diversity[[paste("Pielou")]] <- Pielou
    }
    list_output[[paste("diversity")]] <- list_diversity
  }
  return(list_output)
}

################################ Analyses for data series ###############################################

#'Extension of EcotoneFinder for space/time series
#'
#'@param data A list of dataframes corresponding to the series or a single
#'  dataframe with a series factor. Species must appear as columns.
#'@param dist A vector or column containing the gradient along which the
#'  analysis will be done. Must be of the same lenght as data.
#'@param series If data is a single dataframe, must be the name or the number of
#'  the factor column identifying the series.
#'@param method One of c("dca", "fanny", "vegclust", "diversity", "cmeans", "all").
#'  Tell the function which analysis to perform. See details.
#'@param groups Interger. The desired number of clusters if any of the
#'  clustering method is selected.
#'@param m.exp Integer. The membership exponent for any of the clustering
#'  method.
#'@param standardize Standardize method to apply to the data before further
#'  analysis (for fanny and cmeans). Must be one of decostand methods (see
#'  decostand).
#'@param diversity diversity indice to be calculated. See details.
#'@param na.rm Logical. Should NAs be removed.
#'
#'@details EcotoneFinderSeries is a generalisation of the EcotoneFinder function
#'  to handle space/time series of data. If a dataframe is provided, it will
#'  convert it internally to a named list according to the factor provided by
#'  the series argument. The methods of analysis and standardizations, as well
#'  as the diversity indices are the same as those of the EcotoneFinder
#'  function.
#'
#'@return A list of lists containing the outcomes of the EcotoneFinder function
#'  for each series.
#'
#'@export
#'
#' @examples
#'  ############# Synthetic time series data:
#'  SyntheticTrialSeries <- SyntheticDataSeries(CommunityPool = 40,
#'                                              CommunityNum = 4, SpCo = NULL,
#'                                              Length = 500, SeriesNum = 5,
#'                                              Parameters = list(a=rep(60, 4),
#'                                                                b=c(0,200,350,500),
#'                                                                c=rep(0.03,4)),
#'                                              dev.c = .015,
#'                                              pal = c("#008585", "#B8CDAE", "#E6C186", "#C7522B"),
#'                                              replacement = FALSE,
#'                                              Parameters.repl = TRUE)
#'
#'  EcoTimeSeriesTrial <- EcotoneFinderSeries(data = SyntheticTrialSeries,
#'                                            dist = "Distance", method = "cmeans",
#'                                            series = "Time", groups = 4,
#'                                            standardize = "hellinger", na.rm = TRUE)
#'
#'
#'
EcotoneFinderSeries <- function (data, dist, series = NULL, method = c("dca", "fanny", "vegclust", "diversity", "cmeans", "all"),
                                 groups = NULL, m.exp = 2, standardize = NULL,
                                 diversity = c("shannon", "richness", "expShannon", "pielou", "all"), na.rm = FALSE)
{
  if (is.data.frame(data) == FALSE && is.list(data) == FALSE) {
    stop("data must be either a list of data.frame or a data.frame")
  }
  if (is.data.frame(data) == TRUE) {
    if (is.character(series) == FALSE && is.numeric(series) == FALSE) {
      stop("series must be the name, or number, of a factor column")
    }
    data <- split(data, f = as.data.frame(data)[,which(colnames(data) == series)])
    for (i in 1:length(data)) {
      data[[i]] <- as.data.frame(data[[i]])[, sapply(as.data.frame(data[[i]]), class) != "factor"]
      data[[i]] <- as.data.frame(data[[i]])[, sapply(as.data.frame(data[[i]]), class) != "Date"]
    }
  }
  if (is.character(dist) == FALSE && is.numeric(dist) == FALSE) {
    stop("dist must the name, or number, of the reference column along the sampling dimension")
  }
  TimeSeries <- list()
  for (i in 1:length(data)) {
    TimeSeries[[paste(names(data)[i])]] <- EcotoneFinder(data = data[[i]][-match(dist, names(data[[i]]))], dist = data[[i]][[dist]], method = method,
                                                         m.exp = m.exp, standardize = standardize, groups = groups,
                                                         diversity = diversity, na.rm = na.rm)
  }
  return(TimeSeries)
}

