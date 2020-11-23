########################## Plotting with ggplot grammar #######################################

#' GGplot method for EcotoneFinder
#'
#' @param ecotonefinder list containing elements named in the same way than the
#'   EcotoneFinder function outcomes.
#' @param slope list containing elements named in the same way than the Slope
#'   function outcomes.
#' @param plot.data Logical. Should the data be plotted? Default to FALSE.
#' @param method Analysis method to be plotted from the EcotoneFinder results or
#'   the Slope results. Must be one or several of
#'   "none","dca","fanny","vegclust", "cmeans","diversity",
#'   "dca_slope","fanny_slope","vegclust_slope", "cmeans_slope" or
#'   "diversity_slope".
#' @param axis.number Number of DCA axis to be plotted. Must be between 1 and 4.
#'   Default to 1.
#' @param diversity diversity indice to be plotted, if the method argument
#'   contains "diversity" or "diversity_slope". Must be one or several of "Shannon",
#'   "SpeciesRichness", "ExpShannon", "Pielou", "all"
#' @param facet Character vector of method names indicating how the plot should
#'   be facetted. Can be provided as a list if several methods are to be plotted
#'   on the same facets. Can contain "data" if plot.data = T. If NULL, no facets
#'   are returned. See details.
#' @param col Color palette to be used for plotting. Must be either of length 1,
#'   of the same lenght than the number of facets (when provided), or of the same
#'   length than the number of species (if plot.data = TRUE), or than the number
#'   of groups or axis plotted with the method argument. See details.
#' @param title Main title for the plot
#' @param xlab A title for the x-axis. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param return.plot Logical. If TRUE, the plot is directly plotted. If FALSE,
#'   the plot is stored as a ggplot object. Default to FALSE. See details.
#'
#' @details The ggEcotone function is intended to facilitate the plotting of
#'   EcotoneFinder lists with the use of the ggplot2 grammar. It either directly
#'   print its outputs (if plot = TRUE), or returns a ggplot object that can be
#'   further modified (if plot = FALSE). The latter allows for the addition of
#'   other ggplot2 layers to personalise graphical outputs (see examples).
#'
#'   Facetting options are implemented to allow for the separation of the
#'   different method outputs and facilitate comparisons. The facet parameter
#'   accepts lists, with each element of the list corresponding to a facet and
#'   consisting of the names of the methods to be plotted on that facet.
#'
#'   The col parameter allows for basic control over the colors of the lines.
#'   ggplot internally recycles colour vectors for each new facets, making it
#'   difficult to precisely control colours in facetted plots. Plotting the
#'   outputs on several graphs and arranging them on a grid is the best way to
#'   produce "facetted" plots with different coulour schemes. See examples.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' #### Artificial dataset:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                  SpCo = NULL ,Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.015,3)),
#'                                  pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  ## Analyses:
#'  EcoFinder <- EcotoneFinder(data = SyntheticTrial[,-1],
#'                             dist = SyntheticTrial$Distance,
#'                             method = "all", groups = 3,
#'                             standardize = "hellinger", diversity = "all")
#'
#'  ## Slope calculation:
#'  EcoSlope <- Slope(EcoFinder, method = "all", axis.number = 2,
#'                    diversity = "all")
#'
#'  ## Plots:
#'  \donttest{
#'  require(ggplot2)
#'  require(colorspace)
#'  # Species Distributions and Fuzzy clusters:
#'  Plot <- ggEcotone(EcoFinder, slope = EcoSlope, plot.data = TRUE,
#'                    method = c("cmeans", "fanny"),
#'                    col = c("#D33F6A", "#E99A2C", "#E2E6BD"),
#'                    facet = list(c("data"), c("cmeans", "fanny")),
#'                    title = "Species distribution and fuzzy clusters",
#'                    xlab = "Gradient", ylab = "Membership grades") +
#'    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
#'    theme_bw()
#'  Plot
#'
#'  # Fuzzy clusters & derivatives:
#'  Plot <- ggEcotone(EcoFinder, slope = EcoSlope, plot.data = FALSE,
#'                    method = c("cmeans", "cmeans_slope"),
#'                    col = c("#D33F6A", "#E99A2C", "#E2E6BD"),
#'                    facet = c("cmeans", "cmeans_slope"),
#'                    title = "fuzzy clusters and derivatives",
#'                    xlab = "Gradient", ylab = "Membership grades") +
#'    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
#'    theme_bw()
#'  Plot
#'
#'  # Multiplot layout:
#'  GG1 <- ggEcotone(EcoFinder, slope = EcoSlope, plot.data = TRUE,
#'                   method = c("none"), col = heat_hcl(21), facet = NULL,
#'                   title = "Species distributions", xlab = NULL,
#'                   ylab = "Abundances") +
#'    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
#'    theme_bw()
#'
#'  GG2 <- ggEcotone(EcoFinder, slope = EcoSlope, plot.data = FALSE,
#'                   method = c("cmeans"), col = c("#023FA5", "#BEC1D4", "#D6BCC0"),
#'                   facet = NULL, title = "Fuzzy clusters", xlab = NULL,
#'                   ylab = "Membership grades") +
#'    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
#'    theme_bw()
#'
#'  GG3 <- ggEcotone(EcoFinder, slope = EcoSlope, plot.data = FALSE,
#'                   method = c("diversity"),
#'                   col = c("#26A63A", "#B4B61A"), facet = NULL,
#'                   diversity=c("SpeciesRichness", "ExpShannon"),
#'                   title = "diversity indices", xlab = "Gradient",
#'                   ylab = "Index scores") +
#'    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
#'    theme_bw()
#'
#'  require(Rmisc)
#'  Rmisc::multiplot(GG1,GG2,GG3)
#'  }
#'
ggEcotone <- function (ecotonefinder, slope = NULL, plot.data = FALSE,
                       method = c("none", "dca", "fanny", "vegclust", "cmeans", "diversity",
                                  "dca_slope", "fanny_slope", "vegclust_slope", "cmeans_slope", "diversity_slope"),
                       axis.number = 1, diversity=c("Shannon", "SpeciesRichness", "ExpShannon", "Pielou",
                                                    "SpeciesRichness_slope","Shannon_slope","ExpShannon_slope","Pielou_slope"),
                       facet = NULL, col = "black",
                       title = NULL, xlab = NULL, ylab = NULL,
                       return.plot = TRUE)
{
  m <- c("dca", "fanny", "vegclust", "cmeans", "diversity")
  s <- c("dca_slope", "fanny_slope", "vegclust_slope", "cmeans_slope", "diversity_slope")
  if (any(is.na(names(ecotonefinder[intersect(m, method)])))) {
    stop("ecotonefinder must comntains the elements required in 'method'.")
  }
  if (length(intersect(s, method)) != 0 && any(is.na(names(slope[intersect(s, method)])))) {
    stop("slope must comntains the elements required in 'method'.")
  }
  GG <- list()
  GG[["distance"]] <- ecotonefinder$distance
  if (plot.data == TRUE) {
    GG[["data"]] <- ecotonefinder$data
  }
  if (any(method == "dca")) {
    GG[["dca"]] <- as.data.frame(ecotonefinder$dca$rproj[,1:axis.number])
  }
  if (any(method == "fanny")) {
    GG[["fanny"]] <- as.data.frame(ecotonefinder$fanny$membership)
  }
  if (any(method == "vegclust")) {
    GG[["vegclust"]] <- as.data.frame(ecotonefinder$vegclust$memb)
  }
  if (any(method == "cmeans")) {
    GG[["cmeans"]] <- as.data.frame(ecotonefinder$cmeans$membership)
  }
  if (any(method == "diversity")) {
    GG[["diversity"]] <- as.data.frame(ecotonefinder$diversity)[,intersect(diversity, names(ecotonefinder$diversity))]
  }
  if (any(method == "dca_slope")) {
    GG[["dca_slope"]] <- data.frame(rbind(matrix(NA, floor(slope$window/2), ncol(as.data.frame(slope$dca_slope[,1:axis.number]))),
                                          slope$dca_slope[,1:axis.number, drop = FALSE],
                                          matrix(NA, floor(slope$window/2), ncol(as.data.frame(slope$dca_slope[,1:axis.number])))))
  }
  if (any(method == "fanny_slope")) {
    GG[["fanny_slope"]] <- data.frame(rbind(matrix(NA, floor(slope$window/2), ncol(slope$fanny_slope)),
                                            slope$fanny_slope,
                                            matrix(NA, floor(slope$window/2), ncol(slope$fanny_slope))))
  }
  if (any(method == "vegclust_slope")) {
    GG[["vegclust_slope"]] <- data.frame(rbind(matrix(NA, floor(slope$window/2), ncol(slope$vegclust_slope)),
                                               slope$vegclust_slope,
                                               matrix(NA, floor(slope$window/2), ncol(slope$vegclust_slope))))
  }
  if (any(method == "cmeans_slope")) {
    GG[["cmeans_slope"]] <- data.frame(rbind(matrix(NA, floor(slope$window/2), ncol(slope$cmeans_slope)),
                                             slope$cmeans_slope,
                                             matrix(NA, floor(slope$window/2), ncol(slope$cmeans_slope))))
  }
  if (any(method == "diversity_slope")) {
    GG[["diversity_slope"]] <- data.frame(rbind(matrix(NA, floor(slope$window/2), length(intersect(diversity, names(slope$diversity_slope)))),
                                                as.matrix(do.call(cbind, slope$diversity_slope[intersect(diversity, names(slope$diversity_slope))])),
                                                matrix(NA, floor(slope$window/2), length(intersect(diversity, names(slope$diversity_slope))))))
    colnames(GG[["diversity_slope"]]) <- intersect(diversity, names(slope$diversity_slope))
  }
  if (is.null(facet)) {
    GGdata <- data.frame(do.call(cbind, GG))
  }
  if (is.null(facet) == FALSE) {
    if (plot.data == TRUE) {
      method = c("data", method)
    }
    if(all(method %in% do.call(base::c, as.list(facet))) == FALSE ||
       all(do.call(base::c, as.list(facet)) %in% method) == FALSE) {
      stop("facet and method must match.")
    }
    for (i in method) {
      GG[[i]]$Distance <- GG[["distance"]]
      GG[[i]]$Type <- i
      GG[[i]] <- arrange.vars(GG[[i]], c("Distance" = 1, "Type" = 2))
      if (ncol(GG[[i]]) < purrr::reduce(lapply(GG, ncol), max)) {
        GG[[i]][,(ncol(GG[[i]])+1):purrr::reduce(lapply(GG, ncol), max)] <- NA
      }
      colnames(GG[[i]]) <- c("Distance", "Type", 3:ncol(GG[[i]]))
    }
    GGdata <- data.frame(do.call(rbind, GG[method]))
  }
  if (is.null(facet)) {
    Plot <- ggplot2::ggplot(data=GGdata, ggplot2::aes(GGdata[,1]))
    if (length(col) == 1) {
      col <- rep(col, ncol(GGdata[,-1]))
    } else {
      if (length(col) != ncol(GGdata[,-1])) { stop("col must be either of length == 1,
                                                   or of the length as the number of desired curves") }
      }
    for (i in 2:ncol(GGdata)) {
      Plot <- Plot + ggplot2::geom_line(ggplot2::aes_string(y = GGdata[,i]), col = col[i-1], na.rm = TRUE)
    }
    }
  if (is.null(facet) == FALSE) {
    Fac <- NULL
    for (i in 1:length(facet)) {
      Fac[[i]] <- rep(paste(facet[[i]], collapse ="/"), nrow(GGdata[GGdata$Type %in% facet[[i]],]))
    }
    GGdata$Facet <- factor(do.call(base::c,Fac), levels = unique(do.call(base::c,Fac)))
    Plot <- ggplot2::ggplot(data = GGdata, ggplot2::aes(GGdata[,1]))
    Plot <- Plot + ggplot2::facet_wrap(~ Facet, nrow = length(facet), scales = "free_y")
    if (length(col) == 1) {
      col <- rep(col, length(method))
    }
    if (length(col) == length(method)) {
      names(col) <- method
      levels(names(col)) <- names(col)
      col <- col[order(names(col))]
      for (i in 3:(ncol(GGdata)-1)) {
        Plot <- Plot + ggplot2::geom_line(ggplot2::aes_string(y = GGdata[,i], colour = "Type"), na.rm = TRUE)
      }
    } else {
      if (length(col) == ncol(GGdata)-3) {
        warning("col will be reused in the different facets. Drawing separate plots is an easy way around.")
        for (j in method) {
          for (i in 3:ncol(Filter(function(x)!all(is.na(x)), GG[[j]]))) {
            Plot <- Plot + ggplot2::geom_line(ggplot2::aes_string(y=GGdata[,i]), col = col[i-2], na.rm = TRUE)
          }
        }
      }
    }
    Plot <- Plot + ggplot2::scale_colour_manual(values = purrr::reduce(col,base::c), breaks = method)
  }
  Plot <- Plot + ggplot2::ggtitle(paste(title)) +
    ggplot2::ylab(paste(ylab)) +
    ggplot2::xlab(paste(xlab))
  if (return.plot == FALSE) { print(Plot) } else { return(Plot) }
}
