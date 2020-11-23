################################## Plotting functions for EcotoneFinder ###################################

################################## Plot component for raw analyses #######################################

#' Plotting component for EcotoneFinder
#'
#' @param ecotonefinder A list containing elements named in the same way than
#'   EcotoneFinder function outcomes
#' @param plot.data Logical. Should the data be plotted.
#' @param plot.method Analysis method to be plotted from the EcotoneFinder
#'   analyses. Must be one or several of "none","dca","fanny","vegclust",
#'   "cmeans" or"diversity".
#' @param axis.number Number of axis to plot from the DCA.
#' @param magnification Magnification coefficient for the method. Usefull if the
#'   data are being plotted.
#' @param magnification.diversity Particular magnification for the diversity
#'   indices.
#' @param col.data Colors to be used for the data. See CommunityColor function.
#' @param col.method Colors to be used for the methods.
#' @param title An overall title for the plot. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param xlab A title for the x-axis. See plot.
#' @param na.rm Logical. Should NAs be removed.
#' @param alone Logical. If FALSE, lines are added to an existing plot.
#' @param ... Additional argument to be passed to the plot function.
#'
#' @return A plot with the EcotoneFinder results along the gradient, and
#'   optionally, the data.
#'
#' @export
#'
#' @details Internal component of the PlotEcotone function for the plotting of
#'   the EcotoneFinder analyses. Use PlotEcotone directly for more options.
#'
#' @examples
#'  ######## Artificial dataset & analysis:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 20, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.03,3)),
#'                                  dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  SyntheticEcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all", groups = 3,
#'                                      standardize = "hellinger",
#'                                      diversity = "all")
#'
#'  ### Plot:
#'  require(colorspace)
#'  plotEco(SyntheticEcoFinder, plot.data = FALSE,
#'          plot.method = c("cmeans", "dca"),
#'          axis.number = 2, col.method = terrain_hcl(3))
#'
plotEco <- function (ecotonefinder, plot.data = FALSE, plot.method = c("none","dca","fanny","vegclust","cmeans","diversity"),
                     axis.number = 1, magnification = 20, magnification.diversity = 5, col.data = "black", col.method = c("red","blue"),
                     title = NULL, ylab = "Species", xlab = "Gradient", na.rm = FALSE, alone = TRUE, ...)
{
  if (alone == TRUE) {
    graphics::plot(ecotonefinder$distance,
                   seq(0,max(ecotonefinder$data, na.rm = na.rm),
                       length.out = length(ecotonefinder$distance)),
                   pch = NA,
                   xlab = xlab,
                   ylab = ylab,
                   main = title,
                   ...)
  }
  if (plot.data == TRUE) {
    Names <- names(ecotonefinder$data)
    ecotonefinder$data <- as.data.frame(ecotonefinder$data)
    if (length(col.data) == 1) {
      col.data <- rep(col.data, times = length(Names))
    }
    names(col.data) <- Names
    for (i in Names) {
      graphics::lines(ecotonefinder$distance,
                      ecotonefinder$data[,i],
                      col = col.data[i],
                      ...)
    }
  }
  if (any(plot.method == "dca")) {
    for (i in 1:axis.number) {
      graphics::lines(ecotonefinder$distance,
                      ecotonefinder$dca$rproj[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "fanny")) {
    for (i in 1:ecotonefinder$fanny$k.crisp) {
      graphics::lines(ecotonefinder$distance,
                      ecotonefinder$fanny$membership[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "vegclust")) {
    for (i in 1:length(ecotonefinder$vegclust$memb)) {
      graphics::lines(ecotonefinder$distance,
                      ecotonefinder$vegclust$memb[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "cmeans")) {
    for (i in 1:length(as.data.frame(ecotonefinder$cmeans$membership))) {
      graphics::lines(ecotonefinder$distance,
                      ecotonefinder$cmeans$membership[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "diversity")) {
    for (i in 1:length(ecotonefinder$diversity)) {
      graphics::lines(ecotonefinder$distance,
                      ecotonefinder$diversity[[i]]*magnification.diversity,
                      col = col.method[i],
                      ...)
    }
  }
}

############################### Plot components for slopes ###############################################

#' Plotting component for Slope
#'
#' @param ecotoneslope A list containing elements named in the same way than
#'   Slope function outcomes
#' @param plot.data Logical. Should the data be plotted.
#' @param plot.method Analysis method to be plotted from the Slope results. Must
#'   be one or several of "none", "dca_slope","fanny_slope","vegclust_slope",
#'   "cmeans_slope" or "diversity_slope".
#' @param axis.number Number of axis to plot from the DCA.
#' @param magnification Magnification coefficient for the method. Usefull if the
#'   data are being plotted.
#' @param col.data Colors to be used for the data. See CommunityColor function.
#' @param col.method Colors to be used for the methods.
#' @param title An overall title for the plot. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param xlab A title for the x-axis. See plot.
#' @param na.rm Logical. Should NAs be removed.
#' @param alone Logical. If FALSE, lines are added to an existing plot.
#' @param ... Additional argument to be passed to the plot function.
#'
#' @return A plot with the Slope results along the gradient, and optionally, the
#'   data.
#'
#' @export
#'
#' @details Internal component of the PlotEcotone function for the plotting of
#'   the Slope analyses. Use PlotEcotone directly for more options.
#'
#' @examples
#'  ######## Artificial dataset & analysis:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 20, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.03,3)),
#'                                  dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  SyntheticEcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all", groups = 3,
#'                                      standardize = "hellinger",
#'                                      diversity = "all")
#'
#'  ### Derivatives:
#'  SyntheticSlope <- Slope(SyntheticEcoFinder, method = "all",
#'                          axis.number = 2, diversity = "all")
#'
#' ### Plot:
#' require(colorspace)
#' plotSlope(SyntheticSlope, plot.data = FALSE,
#'           plot.method = c("cmeans_slope", "vegclust_slope"),
#'           col.method = terrain_hcl(3))
#'
#'
plotSlope <- function (ecotoneslope, plot.data = FALSE, plot.method = c("none","dca_slope","fanny_slope","vegclust_slope", "cmeans_slope","diversity_slope"),
                       axis.number = 1, magnification = 500, col.data = "black", col.method = c("red","blue"),
                       title = NULL, ylab = "Species", xlab = "Gradient", na.rm = FALSE, alone = TRUE, ...)
{
  if (alone == TRUE) {
    graphics::plot(ecotoneslope$distance,
                   seq(0,max(ecotoneslope$data, na.rm = na.rm),
                       length.out = length(ecotoneslope$distance)),
                   pch = NA,
                   xlab = xlab,
                   ylab = ylab,
                   main = title,
                   ...)
  }
  if (plot.data == TRUE) {
    Names <- names(ecotoneslope$data)
    ecotoneslope$data <- as.data.frame(ecotoneslope$data)
    if (length(col.data) == 1) {
      col.data <- rep(col.data, times = length(Names))
    }
    names(col.data) <- Names
    for (i in Names) {
      graphics::lines(ecotoneslope$distance,
                      ecotoneslope$data[,i],
                      col = col.data[i],
                      ...)
    }
  }
  if (any(plot.method == "dca_slope")) {
    for (i in 1:axis.number) {
      graphics::lines(c((min(ecotoneslope$distance)+as.integer(ecotoneslope$window/2)):(max(ecotoneslope$distance)-as.integer(ecotoneslope$window/2))),
                      ecotoneslope$dca_slope[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "fanny_slope")) {
    for (i in 1:ecotoneslope$groups) {
      graphics::lines(c((min(ecotoneslope$distance)+as.integer(ecotoneslope$window/2)):(max(ecotoneslope$distance)-as.integer(ecotoneslope$window/2))),
                      ecotoneslope$fanny_slope[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "vegclust_slope")) {
    for (i in 1:ecotoneslope$groups) {
      graphics::lines(c((min(ecotoneslope$distance)+as.integer(ecotoneslope$window/2)):(max(ecotoneslope$distance)-as.integer(ecotoneslope$window/2))),
                      ecotoneslope$vegclust_slope[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "cmeans_slope")) {
    for (i in 1:ecotoneslope$groups) {
      graphics::lines(c((min(ecotoneslope$distance)+as.integer(ecotoneslope$window/2)):(max(ecotoneslope$distance)-as.integer(ecotoneslope$window/2))),
                      ecotoneslope$cmeans_slope[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "diversity_slope")) {
    for (i in 1:length(ecotoneslope$diversity)) {
      graphics::lines(c((min(ecotoneslope$distance)+as.integer(ecotoneslope$window/2)):(max(ecotoneslope$distance)-as.integer(ecotoneslope$window/2))),
                      ecotoneslope$diversity[[i]]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
}

################################## Plot component for Environment #######################################

#' Plotting component for EcotoneFinder when run on environmental data
#'
#' @param env A list containing elements named in the same way than
#'   EcotoneFinder function outcomes.
#' @param plot.data Logical. Should the data be plotted.
#' @param plot.method Analysis method to be plotted from the EcotoneFinder
#'   analyses. Must be one or several of "none","dca","fanny","vegclust",
#'   "cmeans" or"diversity".
#' @param axis.number Number of axis to plot from the DCA.
#' @param magnification Magnification coefficient for the method. Usefull if the
#'   data are being plotted.
#' @param magnification.diversity Particular magnification for the diversity
#'   indices.
#' @param col.data Colors to be used for the data. See CommunityColor function.
#' @param col.method Colors to be used for the methods.
#' @param title An overall title for the plot. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param xlab A title for the x-axis. See plot.
#' @param na.rm Logical. Should NAs be removed.
#' @param alone Logical. If FALSE, lines are added to an existing plot.
#' @param ... Additional argument to be passed to the plot function.
#'
#' @return A plot with the EcotoneFinder results along the gradient, and
#'   optionally, the data.
#'
#' @export
#'
#' @details Internal component of the PlotEcotone function for the plotting of
#'   the EcotoneFinder analyses. Use PlotEcotone directly for more options.
#'   The "diversity" method is still implemented, but will send a warning as it
#'   may not be relevant for environmental data.
#'
#' @examples
#'  ######## Artificial dataset & analysis:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 20, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.03,3)),
#'                                  dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  SyntheticEcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all", groups = 3,
#'                                      standardize = "hellinger",
#'                                      diversity = "all")
#'
#'  ### Plot:
#'  require(colorspace)
#'  plotEnv(SyntheticEcoFinder, plot.data = FALSE,
#'          plot.method = c("cmeans", "dca"),
#'          axis.number = 2, col.method = terrain_hcl(3))
#'
plotEnv <- function (env, plot.data = FALSE, plot.method = c("none","dca","fanny","vegclust","cmeans","diversity"),
                     axis.number = 1, magnification = 20, magnification.diversity = 5, col.data = "black", col.method = c("red","blue"),
                     title = NULL, ylab = "Species", xlab = "Gradient", na.rm = FALSE, alone = TRUE, ...)
{
  if (alone == TRUE) {
    graphics::plot(env$distance,
                   seq(0,max(env$data, na.rm = na.rm),
                       length.out = length(env$distance)),
                   pch = NA,
                   xlab = xlab,
                   ylab = ylab,
                   main = title,
                   ...)
  }
  if (plot.data == TRUE) {
    Names <- names(env$data)
    env$data <- as.data.frame(env$data)
    if (length(col.data) == 1) {
      col.data <- rep(col.data, times = length(Names))
    }
    names(col.data) <- Names
    for (i in Names) {
      graphics::lines(env$distance,
                      env$data[,i],
                      col = col.data[i],
                      ...)
    }
  }
  if (any(plot.method == "dca")) {
    for (i in 1:axis.number) {
      graphics::lines(env$distance,
                      env$dca$rproj[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "fanny")) {
    for (i in 1:env$fanny$k.crisp) {
      graphics::lines(env$distance,
                      env$fanny$membership[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "vegclust")) {
    for (i in 1:length(env$vegclust$memb)) {
      graphics::lines(env$distance,
                      env$vegclust$memb[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "cmeans")) {
    for (i in 1:length(as.data.frame(env$cmeans$membership))) {
      graphics::lines(env$distance,
                      env$cmeans$membership[,i]*magnification,
                      col = col.method[i],
                      ...)
    }
  }
  if (any(plot.method == "diversity")) {
    for (i in 1:length(env$diversity)) {
      message("diversity indices may not be relevants on environmental data.")
      graphics::lines(env$distance,
                      env$diversity[[i]]*magnification.diversity,
                      col = col.method[i],
                      ...)
    }
  }
}


##################################### Plot Ecotone #####################################################

#' Plot method for EcotoneFinder
#'
#' @param data A list containing elements named in the same way than
#'   EcotoneFinder function outcomes
#' @param slope A list containing elements named in the same way than Slope
#'   function outcomes
#' @param env A list containing elements named in the same way than
#'   EcotoneFinder function outcomes. Usefull if EcotoneFinder has been
#'   run on environmental data and the outcomes are to be compared with
#'   the outcomes from the community matrix. Can also be used to compare
#'   results from two different community matrices, if the x-axis are
#'   similar.
#' @param plot.data Logical. Should the data be plotted.
#' @param plot.method Analysis method to be plotted from the EcotoneFinder
#'   results or the Slope results. Must be one or several of
#'   "none","dca","fanny","vegclust", "cmeans","diversity",
#'   "dca_slope","fanny_slope","vegclust_slope", "cmeans_slope" or
#'   "diversity_slope".
#' @param axis.number Number of axis to plot from the DCA.
#' @param magnification Magnification coefficient for the method. Usefull if the
#'   data are being plotted.
#' @param magnification.diversity Particular magnification for the diversity
#'   indices.
#' @param magnification.slope Magnification coefficient for the Slope.
#' @param col.data Colors to be used for the data. See CommunityColor function.
#' @param col.method Colors to be used for the methods (for data).
#' @param col.slope Colors to be used for the methods (for slope).
#' @param col.env Colors to be used for the methods (for env).
#' @param title An overall title for the plot. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param xlab A title for the x-axis. See plot.
#' @param na.rm Logical. Should NAs be removed.
#' @param ... Additional argument to be passed to the plot function.
#'
#' @details The plotEcotone function is intended for easy visualisation of the
#'   results of the EcotoneFinder function and the Slope function along the
#'   sampling gradient. It also provide a way to plot the original species data
#'   - for comparison - along with magnification coefficients for the different
#'   method or slopes in order to facilitate visualization. Please note that
#'   large sets of species may be confusing to plot. The CommunityColor function
#'   is also provided to help with data coloring and easy visualisation.
#'
#'
#' @return A plot corresponding to the plotEcotone arguments.
#'
#' @export
#'
#' @examples
#'  ######## Artificial dataset & analysis:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 20, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.03,3)),
#'                                  dev.c = .015, pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  SyntheticEcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all", groups = 3,
#'                                      standardize = "hellinger",
#'                                      diversity = "all")
#'
#'  ##### Assigning colors to communities:
#'  SyntheticColor <- CommunityColor(SyntheticEcoFinder, pal = "diverge_hcl",
#'                                   method = "cmeans")
#'
#'  ###### Computing the derivatives:
#'  SyntheticSlope <- Slope(SyntheticEcoFinder, method = "all",
#'                          axis.number = 2, diversity = "all")
#'
#'  ####### Plot the derivative of the FCM with the synthetic species data:
#'  require(colorspace)
#'  plotEcotone(slope = SyntheticSlope, plot.data = TRUE,
#'              plot.method = c("cmeans_slope"), axis.number = 2,
#'              col.method = terrain_hcl(3), col.data = SyntheticColor)
#'
#'  ####### Plot the derivative and the FCM:
#'  require(colorspace)
#'  plotEcotone(data = SyntheticEcoFinder, slope = SyntheticSlope,
#'              plot.data = TRUE,
#'              plot.method = c("cmeans", "cmeans_slope"),
#'              axis.number = 2, col.method = terrain_hcl(3),
#'              col.data = SyntheticColor)
#'
#'
#'
plotEcotone <- function (data = NULL, slope = NULL, env = NULL, plot.data = FALSE, plot.method=c("none","dca","fanny","vegclust","cmeans","diversity", "dca_slope","fanny_slope","vegclust_slope", "cmeans_slope", "diversity_slope"),
                         axis.number = 1, magnification = 20, magnification.diversity = 5, magnification.slope = 500, col.data = "black", col.method = c("red","blue"), col.slope = c("darkgreen","green"), col.env = c("orange","gold"),
                         title = NULL, ylab = "Species", xlab = "Gradient", na.rm = FALSE, ...)
{
  alone <- FALSE
  if (is.null(data) && is.null(slope)) {
    stop("At least one of data or slope must be provided")
  }
  m <- c("dca", "fanny", "vegclust", "cmeans", "diversity")
  s <- c("dca_slope", "fanny_slope", "vegclust_slope", "cmeans_slope", "diversity_slope")
  if (any(is.na(names(data[intersect(m, plot.method)])))) {
    stop("data must comntains the elements required in 'method'.")
  }
  if (length(intersect(s, plot.method)) != 0 && any(is.na(names(slope[intersect(s, plot.method)])))) {
    stop("slope must comntains the elements required in 'method'.")
  }
  if (is.null(data) == FALSE) {
    dist <- data$distance
    max <- max(data$data, na.rm = na.rm)
    l <- length(data$distance)
  }
  if (is.null(slope) == FALSE) {
    dist <- slope$distance
    max <- max(slope$data, na.rm = na.rm)
    l <- length(slope$distance)
  }
  if (is.null(env) == FALSE) {
    dist <- env$distance
    max <- max(env$data, na.rm = na.rm)
    l <- length(env$distance)
  }
  graphics::plot(dist,
                 seq(0, max, length.out = l),
                 pch = NA,
                 xlab = xlab,
                 ylab = ylab,
                 main = title,
                 ...)
  if (is.null(data) == FALSE) {
    plotEco(ecotonefinder = data, plot.data = plot.data, plot.method = plot.method,
            axis.number = axis.number, magnification = magnification, magnification.diversity = magnification.diversity,
            col.data = col.data, col.method = col.method,
            title = title, ylab = ylab, xlab = xlab, na.rm = na.rm, alone = alone, ...)
  }
  if (is.null(slope) == FALSE) {
    plotSlope(ecotoneslope = slope, plot.data = plot.data, plot.method = plot.method,
              axis.number = axis.number, magnification = magnification.slope,
              col.data = col.data, col.method = col.slope,
              title = title, ylab = ylab, xlab = xlab, na.rm = na.rm, alone = alone, ...)
  }
  if (is.null(env) == FALSE) {
    plotEnv(env = env, plot.data = plot.data, plot.method = plot.method,
            axis.number = axis.number, magnification = magnification,
            col.data = col.data, col.method = col.env,
            title = title, ylab = ylab, xlab = xlab, na.rm = na.rm, alone = alone, ...)
  }
}

###################################### Colors for communities ############################################

#' Tool to assign color to species distribution plots given fuzzy clustering
#' results.
#'
#' @param ecotonefinder A list containing elements named in the same way than
#'   EcotoneFinder function outcomes. Must contain at least one of vegclust or
#'   cmeans results.
#' @param method Which fuzzy clustering to use. Either vegclust or cmeans.
#' @param pal Which palette to use for color picking. Chosen from the colorspace
#'   package.
#'
#' @details CommunityColor creates a color vector that can be used by plotting
#'   functions. It assigns colors to species of a community matrix given the
#'   results of vegclust or cmeans analyses (using cmeans$centers or
#'   vegclust$mobileCenters). Species are assigned to a color according to the
#'   cluster centroid for which they have their highest membership value (see
#'   Bandelj et al., 2012).
#'
#'   The palette must be one of the colorspace package.
#'
#'
#' @return A vector of color names from the palette in the pal argument, of the
#'   same length and in the same order than the species columns of the provided
#'   data.
#'
#' @export
#'
#' @examples
#'  ######## Artificial dataset & analysis:
#'  SyntheticTrial <- SyntheticData(SpeciesNum = 27, CommunityNum = 3,
#'                                  SpCo = NULL, Length = 500,
#'                                  Parameters = list(a=rep(60, 3),
#'                                                    b=c(0,250,500),
#'                                                    c=rep(0.03,3)),
#'                                  dev.c = .015,
#'                                  pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'  SyntheticEcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
#'                                      dist = SyntheticTrial$Distance,
#'                                      method = "all", groups = 3,
#'                                      standardize = "hellinger",
#'                                      diversity = "all")
#'
#'  ##### Assigning colors to communities:
#'  SyntheticColor <- CommunityColor(SyntheticEcoFinder, pal = "diverge_hcl",
#'                    method = "cmeans")
#'
#'  #### Plotting:
#'  plotEcotone(data = SyntheticEcoFinder, plot.data = TRUE, plot.method = "none",
#'              col.data = SyntheticColor)
#'
#'
#'
CommunityColor <- function (ecotonefinder, method = c("vegclust", "cmeans"),
                            pal = c("diverge_hcl", "terrain_hcl", "sequential_hcl", "rainbow_hcl"))
{
  # Warnings:
  if (method == "vegclust") {
    if (is.null(ecotonefinder$vegclust)) {
      stop("Community color require vegclust method in EcotoneFinder")
    }
    Cluster <- rep(NA, length(ecotonefinder$vegclust$mobileCenters))
    for (i in 1:length(ecotonefinder$vegclust$mobileCenters)) {
      Cluster[i] <- which.max(ecotonefinder$vegclust$mobileCenters[,i])
    }
  }
  if (method == "cmeans") {
    if (is.null(ecotonefinder$cmeans)) {
      stop("Community color require cmeans method in EcotoneFinder")
    }
    Cluster <- rep(NA, length(ecotonefinder$cmeans$centers))
    for (i in 1:length(ecotonefinder$cmeans$centers)) {
      Cluster[i] <- which.max(ecotonefinder$cmeans$centers[,i])
    }
  }

  # Color patterns:
  if (pal == "diverge_hcl") {
    colvec <- colorspace::diverge_hcl(max(Cluster))
  }
  if (pal == "terrain_hcl") {
    colvec <- colorspace::terrain_hcl(max(Cluster))
  }
  if (pal == "sequential_hcl") {
    colvec <- colorspace::sequential_hcl(max(Cluster))
  }
  if (pal == "rainbow_hcl") {
    colvec <- colorspace::rainbow_hcl(max(Cluster))
  }
  CommunityColor <- rep(NA, length(Cluster))
  for (i in 1:length(Cluster)) {
    CommunityColor[i] <- colvec[Cluster[i]]
  }
  return(CommunityColor)
}

