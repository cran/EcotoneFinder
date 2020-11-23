################################## Functions to build synthetic data ######################################

################################# Unique set of data ####################################################

#' Create synthetic gaussian-shaped species abundance data
#'
#' @param SpeciesNum An integer giving the total number of species in the
#'   synthetic data.
#' @param CommunityNum An integer giving the number of communities to be
#'   synthetised.
#' @param Length The lenght of the gradient. Corresponds to the x-axis in a
#'   plot.
#' @param SpCo The ratio of species per communities. If NULL, species will be
#'   spread evenly between communities with additional species in the last
#'   community if the quotient is not an integer. When specified, SpCo must be a
#'   vector of lenght equal to CommunityNum and whose sum is equal to SpeciesNum
#' @param Parameters A list containing the parameters (a, b and c) for the
#'   gaussians. Each parameter must be specified for each community. See
#'   Details.
#' @param dev.a The deviation around parameter a for the gaussian in a
#'   community. If 0 all species curve in the comunity will have the same a
#'   parameter.
#' @param dev.b The deviation around parameter b for the gaussian in a
#'   community. If 0 all species curve in the comunity will have the same b
#'   parameter.
#' @param dev.c The deviation around parameter a for the gaussian in c
#'   community. If 0 all species curve in the comunity will have the same c
#'   parameter.
#' @param down.limit The limit under which the gaussian curve will be rounded
#'   down to 0. The default is 1.
#' @param pal The color palette to be used. Species curves are colored according
#'   to communities. Either a colorspace palette or a vector of the same lenght
#'   as the number of species.
#' @param xlab A title for the x-axis. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param title An overall title for the plot. See plot.
#'
#' @details The SyntheticData function is intended for the creation of
#'   articficial dataset to test ecological patterns along gradients. The
#'   gaussian curves that it computes are of the form: \eqn{
#'   a*exp(-(((x-b)^2)/2*(c^2))) } The parameters can be interpreted as follow:
#'   a is the maximum height of the gaussian on the y-axis, b is the center of
#'   the gaussian on the x-axis and c is the steepness of the slopes on each
#'   side of the maximum. The gaussians create a set of continuous data that are
#'   akin to abundances. As gaussians of this type cannot reach 0, any value
#'   that is below the down.limit (default is 1) is rounded down to 0.
#'
#' @return SyntheticData returns a dataset with numbered species (sp.1, sp.2,
#'   ...) as columns. It also plot the obtained data.
#'
#' @export
#'
#' @examples
#' ### 3 distinct communities comprising a total of 21 species
#' SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                 SpCo = NULL, Length = 500,
#'                                 Parameters = list(a=rep(60, 3),
#'                                                   b=c(0,250,500),
#'                                                   c=rep(0.03,3)),
#'                                 pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#' ### 3 distinct communities with uneven species numbers
#' SyntheticTrial <- SyntheticData(SpeciesNum = 21, CommunityNum = 3,
#'                                 SpCo = c(5, 10, 6), Length = 500,
#'                                 Parameters = list(a=rep(60, 3),
#'                                                   b=c(0,250,500),
#'                                                   c=rep(0.03,3)),
#'                                 pal = c("#008585", "#FBF2C4", "#C7522B"))
#'
#'
#'
#'
SyntheticData <- function (SpeciesNum, CommunityNum, Length = 100, SpCo = NULL,
                           Parameters = list(a = NULL, b = NULL, c = NULL), dev.a = 10, dev.b = 10, dev.c = 0, down.limit = 1,
                           pal = NULL,
                           xlab = "Gradient", ylab = "Synthetic species", title = "Synthetic data")
{
  if (is.null(SpeciesNum) || is.null(CommunityNum)) {
    stop("Both SpeciesNum and CommuityNum must be specified.")
  }
  Specieslist <- rep(NA,SpeciesNum)
  for (i in 1:SpeciesNum) {
    Specieslist[i] <- paste("Sp",i)
  }
  Data <- data.frame(matrix(NA, nrow = Length, ncol = SpeciesNum+1, dimnames = list(1:Length,c("Distance",Specieslist))))
  # Distance:
  Data$Distance <- seq(1:Length)
  # Assigning species in communities:
  if (is.null(SpCo)) {
    ratio <- as.integer(SpeciesNum/CommunityNum)
    SpCo <- rep(ratio, CommunityNum)
    if ((SpeciesNum/CommunityNum) - (as.integer(SpeciesNum/CommunityNum)) != 0) {
      SpCo[CommunityNum] <- SpCo[CommunityNum]+(SpeciesNum-CommunityNum*ratio)
      # Warnings:
      message(paste("uneven munber of species / communities.
                    Computed with ", SpeciesNum-CommunityNum*ratio, "extra species in the last community \n
                    Otherwise you can specify the desired number of species / community with the SpCo argument"))
    }
    }
  if (length(SpCo) != CommunityNum) {
    stop("SpCo must be of the length of CommunityNum")
  }
  if (sum(SpCo) != SpeciesNum) {
    stop("missing species in SpCo: the sum of SpCo  must be equal to SpeciesNum")
  }
  SpList <- list()
  for (i in 1:CommunityNum) {
    SpList[[i]] <- c(1:SpCo[i])
  }
  SeqList <- list()
  for(i in 1:CommunityNum) {
    SeqList[[1]] <- 0
    SeqList[[i+1]] <- c(1:SpCo[i])
  }
  ## Curve Loop:
  if (length(Parameters) != 3) {
    stop("missing parameters for gaussian curve")
  }
  if (length(Parameters[[1]]) != length(Parameters[[2]]) ||
      length(Parameters[[1]]) != length(Parameters[[3]]) ||
      length(Parameters[[1]]) != CommunityNum) {
    stop("unequal parameter list for gaussian curve or parameter list is not the same length as CommunityNum")
  }
  gaussian <- function(x) a*exp(-(((x-b)^2)/2*(c^2)))
  for (k in 1:CommunityNum) {
    Aseq <- seq(Parameters[[1]][[k]]-dev.a, Parameters[[1]][[k]]+dev.a, by=1)
    Bseq <- seq(Parameters[[2]][[k]]-dev.b, Parameters[[2]][[k]]+dev.b, by=1)
    Cseq <- seq(Parameters[[3]][[k]]-dev.c, Parameters[[3]][[k]]+dev.c, by=0.005)
    for (i in SpList[[k]]) {
      if (is.null(dev.a) == T) { a <- Aseq } else { a <- sample(Aseq, 1) }
      if (is.null(dev.b) == T) { b <- Bseq } else { b <- sample(Bseq, 1) }
      if (is.null(dev.c) == T) { c <- Cseq } else { c <- sample(Cseq, 1) }
      Curve <- curveNoPlot(gaussian, from = 1, to = Length, n = Length)
      Curve$y[Curve$y<down.limit] <- 0
      Data[,i+purrr::reduce(lapply(SeqList[1:k], max), sum)+1] <- Curve$y
    }
  }
  # Plot loop:
  # Color pattern:
  if  (length(pal) != SpeciesNum && length(pal) != CommunityNum) {
    stop("pal should either of the same lenght as CommnityNum or of the same lenght as SpeciesNum.")
  }
  if (length(pal) == CommunityNum) {
    colvec <- rep(pal, times=SpCo)
  }
  if (length(pal) == SpeciesNum) {
    colvec <- pal
  }
  graphics::plot(Data$Distance,
       seq(0, max(Data[,-1]), length.out = length(Data$Distance)),
       pch = NA,
       xlab = xlab,
       ylab = ylab,
       main = title)
  for (i in 1:SpeciesNum) {
    graphics::lines(Data$Distance,
          Data[,i+1],
          col = colvec[i])
  }
  return(Data)
}

################################## Data Series ###########################################################

#' Synthetic data for Space/Time series
#'
#' @param CommunityPool Total number of species per series. must either be of
#'   length 1 (if the pool of species is of the same length for all communities)
#'   or of length equal to CommunityNum (specifying the species pool for each
#'   community)
#' @param CommunityNum An integer giving the number of communities.
#' @param Length The lenght of the gradient. Corresponds to the x-axis in a
#'   plot.
#' @param SpCo The ratio of species per communities. When replacement = T, SpCo
#'   is computed from CommunityPool and range.repl.
#' @param SeriesNum The number of series to be synthetised.
#' @param replacement Logical. Should the species of a community in a given
#'   series be a sample of the possible number of species in this community. See
#'   Details.
#' @param range.repl An integer. Gives the possible range of species turnover in
#'   a community.
#' @param Parameters A list containing the parameters (a, b and c) for the
#'   gaussians. Each parameter must be specified for each community. See
#'   Details.
#' @param Parameters.repl Logical. Sould the parameters vary between series.
#' @param dev.a The deviation around parameter a for the gaussian in a
#'   community. If 0 all species curve in the comunity will have the same a
#'   parameter.
#' @param dev.b The deviation around parameter b for the gaussian in a
#'   community. If 0 all species curve in the comunity will have the same b
#'   parameter.
#' @param dev.c The deviation around parameter a for the gaussian in c
#'   community. If 0 all species curve in the comunity will have the same c
#'   parameter.
#' @param Parameters.range An integer. Gives the possible range of the parameter
#'   variations between two series. The range of variation is the parameter
#'   diveded by the parameters range. Default is 10.
#' @param displacement Numeric matrix to control the direction of the changes in
#'   the mean values of the parameters over the series. See examples.
#' @param pal The color palette to be used. Species curves are colored according
#'   to communities.
#' @param xlab A title for the x-axis. See plot.
#' @param ylab A title for the y-axis. See plot.
#' @param title An overall title for the plot. See plot.
#'
#' @details SyntheticDataSeries is a extention of the SyntheticData function and
#'   intended to produce easy and consistent space/time series of artificial
#'   ecological community datasets. The series of dataframes are stored in an
#'   object of class list. The replacement and Parameters.repl arguments allow
#'   the user to choose wether or not the number of species and their
#'   distribution curves should vary among the different series. range.repl is
#'   an integer that define the boundaries  of the interval in which the number
#'   of species can vary - such as the number of species per community is :
#'   \emph{sample([[CommunityPool/CommunityNum - range.repl ;
#'   CommunityPool/CommunityNum + range.repl]])} The number it takes is then
#'   used to define SpCo. Given the difference of scales between the 3
#'   parameters a, b and c, the Parameters.range argument controls the
#'   variations of the parameters by divison so that they correspond to :
#'   \emph{sample([[Parameters - Parameters/Parameters.range ; Parameters +
#'   Parameters/Parameters.range]])} The obtained parameters are then used by
#'   the internal SyntheticData function as base parameters for a given series,
#'   on which dev.a, dev.b and dev.c will apply. The other arguments are
#'   equivalent to those of SyntheticData.
#'
#' @return SyntheticDataSeries returns a list of datasets with numbered species
#'   (sp.1, sp.2, ...) as columns. The list has length = SeriesNum. It also plot
#'   the obtained data.
#'
#' @export
#'
#' @examples
#' ##### 5 datasets of 40 species spread on 4 communities without turnover
#' ##### on the number of species nor variations in their distribution:
#'
#' SyntheticTrialSeries <- SyntheticDataSeries(CommunityPool = 40,
#'                                             CommunityNum = 4, SpCo = NULL,
#'                                             Length = 500, SeriesNum = 5,
#'                                             Parameters = list(a=rep(60, 4),
#'                                                               b=c(0,200,350,500),
#'                                                               c=rep(0.03,4)),
#'                                             pal = c("#008585", "#B8CDAE", "#E6C186", "#C7522B"),
#'                                             replacement = FALSE,
#'                                             Parameters.repl = FALSE)
#'
#' ##### 5 datasets of 40 species spread on 4 communities with species turnover
#' ##### and variations in their distributions along the gradient:
#'
#' SyntheticTrialSeries <- SyntheticDataSeries(CommunityPool = 40,
#'                                             CommunityNum = 4, SpCo = NULL,
#'                                             Length = 500, SeriesNum = 5,
#'                                             Parameters = list(a=rep(60, 4),
#'                                                               b=c(0,200,350,500),
#'                                                               c=rep(0.03,4)),
#'                                             pal = c("#008585", "#B8CDAE", "#E6C186", "#C7522B"),
#'                                             replacement = TRUE,
#'                                             Parameters.repl = TRUE)
#'
#' ##### With a displacement matrix to control the direction of the changes
#' ##### between series:
#'
#' # Dispalcement matrix (Parameters x Communities):
#' disp <- matrix(data=c(0,0,0,
#'                       0,35,-0.0007,
#'                       0,10,0), nrow = 3, ncol = 3)
#'
#' Series <- SyntheticDataSeries(CommunityPool = 60, CommunityNum = 3, Length = 500,
#'                                SeriesNum = 5, replacement = FALSE, SpCo = c(15,15,30),
#'                                Parameters = list(a = c(60,60,60),
#'                                                  b = c(-50,-50,400),
#'                                                  c = c(0.01, 0.01, 0.01)),
#'                                dev.a=30, dev.b=40, dev.c=0,
#'                                displacement = disp,
#'                                pal = c(rep("#008585",15), rep("#FBF2C4",15),
#'                                        rep("#C7522B",30)))
#'
SyntheticDataSeries <- function(CommunityPool, CommunityNum, Length = 100, SpCo = NULL,
                                SeriesNum, replacement = TRUE, range.repl = as.integer(CommunityPool/5),
                                Parameters = list(a = NULL,b = NULL,c = NULL), Parameters.repl = TRUE, dev.a = 10, dev.b = 10, dev.c = 0,
                                Parameters.range = 10, displacement = NULL,
                                pal = NULL,
                                xlab = "Gradient", ylab = "Synthetic species", title = "Synthetic data")
{
  Data <- list()
  Pool <- list()
  if (length(CommunityPool) !=1 && length(CommunityPool) != CommunityNum) {
    stop("CommunityPool must either be of length 1 (if the pool of species is of the same length for all communities)\n
         or of length equal to CommunityNum (specifying the species pool for each community)")
  }
  if (length(CommunityPool) == 1) {
    for (k in 1:CommunityNum) {
      for (i in 1:CommunityPool) {
        Pool[[paste("C",k, sep="")]][i] <- paste("Sp.",i+(CommunityPool*(k-1)), sep = "")
      }
    }
  }
  if (length(CommunityPool) == CommunityNum) {
    for (k in 1:CommunityNum) {
      for (i in 1:CommunityPool[k]) {
        Pool[[paste("C",k,sep = "")]][i] <- paste("Sp.", i+(CommunityPool[k]*(k-1)), sep = "")
      }
    }
  }
  for (j in 1:SeriesNum) {
    if (replacement == TRUE) {
      message("When replacement = T, the number of species per community is computed from CommunityPool and range.repl")
      SpCo <- rep(NA, CommunityNum)
      C <- list()
      for (k in 1:CommunityNum) {
        C[[paste("C",k,sep = "")]] <- sample(Pool[[paste("C",k,sep = "")]], sample(c(length(Pool[[paste("C",k,sep = "")]])-range.repl):length(Pool[[paste("C",k,sep = "")]]),1))
        C[[paste("C",k,sep = "")]] <- C[[paste("C",k,sep = "")]][order(as.numeric(gsub("\\D", "", C[[paste("C",k,sep = "")]])))]
        SpCo[k] <- length(C[[paste("C",k,sep = "")]])
      }
      SpeciesNum <- sum(SpCo)
    }
    if (replacement == FALSE) {
      C <- Pool
      if (is.null(SpCo)) {
        SpeciesNum <- CommunityPool*CommunityNum
        ratio <- as.integer(SpeciesNum/CommunityNum)
        SpCo <- rep(ratio, CommunityNum)
        if ((SpeciesNum/CommunityNum) - (as.integer(SpeciesNum/CommunityNum)) != 0) {
          SpCo[CommunityNum] <- SpCo[CommunityNum]+(SpeciesNum-CommunityNum*ratio)
          # Warnings:
          message(paste("uneven munber of species / communities.
                        Computed with ", CommunityPool-CommunityNum*ratio, "extra species in the last community \n
                        Otherwise you can specify the desired number of species / community with the SpCo argument"))
        }
        } else {
          SpeciesNum <- sum(SpCo)
        }
    }
    if (length(Parameters) != 3) {
      stop("missing parameters for gaussian curve")
    }
    if (length(Parameters[[1]]) != length(Parameters[[2]]) ||
        length(Parameters[[1]]) != length(Parameters[[3]]) ||
        length(Parameters[[1]]) != CommunityNum) {
      stop("unequal parameter list for gaussian curve or parameter list is not the same length as CommunityNum")
    }
    if (Parameters.repl == TRUE) {
      P <- list()
      Range <- lapply(Parameters, "/", Parameters.range)
      for (i in 1:CommunityNum) {
        P[["a"]][[i]] <- sample(c((Parameters[[1]][[i]]-Range[[1]][[i]]):(Parameters[[1]][[i]]+Range[[1]][[i]])),1)
        P[["b"]][[i]] <- sample(c((Parameters[[2]][[i]]-Range[[2]][[i]]):(Parameters[[2]][[i]]+Range[[2]][[i]])),1)
        P[["c"]][[i]] <- sample(c((Parameters[[3]][[i]]-Range[[3]][[i]]):(Parameters[[3]][[i]]+Range[[3]][[i]])),1)
      }
    }
    if (Parameters.repl == FALSE) {
      P <- Parameters
    }
    if (is.null(displacement) == FALSE) {
      P <- as.numeric(purrr::reduce(P, base::c)) + as.numeric(purrr::reduce(t(((j-1)*displacement)), base::c))
      dim(P) <- dim(displacement)
      P <- (split(t(P), c("a", "b", "c")))
    }
    Data[[paste("Time",j,sep = "")]] <- SyntheticData(SpeciesNum = SpeciesNum, CommunityNum = CommunityNum, SpCo = SpCo, Length = Length,
                                                      Parameters = P, dev.a = dev.a, dev.b = dev.b, dev.c = dev.c, pal = pal,
                                                      xlab = xlab, ylab = ylab, title = title)
    }
  return(Data)
}
