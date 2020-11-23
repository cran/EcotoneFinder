############################## Miscellenious functions #####################################

################################# cbind.na ##################################################

#' qpcR cbind.na method.
#'
#' @param ... (generalized) vectors or matrices. See base::cbind
#' @param deparse.level integer controlling the construction
#'   of labels in the case of non-matrix-like arguments. See base::cbind
#'
#' @return  a matrix combining the ... arguments column-wise.
#'
#' @export
#'
#' @examples
#' ### Vectors:
#' a <- c(rep(1, 5), NA, seq(1:5))
#' b <- c(rep(1, 4), NA, seq(1:7))
#'
#' # Complete shorter vector with NAs:
#' cbindna(a,b)
#'
cbindna <- function (..., deparse.level = 1)
{
  na <- nargs() - (!missing(deparse.level))
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na == 0)
    return(NULL)
  if (na == 1) {
    if (isS4(..1))
      return(methods::cbind2(..1))
    else return(matrix(...))
  }
  if (deparse.level) {
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level ==
            2)
          deparse(r)
      }
      else r
    }
  }
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, nrow))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    if (deparse.level) {
      if (fix.na)
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl)))
        iV <- iV & (nmi == "")
      ii <- if (fix.na)
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
          names(argl)[i] <- nmi
      }
    }
    nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
    maxRow <- max(nRow, na.rm = TRUE)
    argl <- lapply(argl, function(x) if (is.null(nrow(x)))
      c(x, rep(NA, maxRow - length(x)))
      else rbindna(x, matrix(, maxRow - nrow(x), ncol(x))))
    r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
  }
  d2 <- dim(r)
  r <- methods::cbind2(argl[[1]], r)
  if (deparse.level == 0)
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2)
    return(r)
  Ncol <- function(x) {
    d <- dim(x)
    if (length(d) == 2L)
      d[2L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(colnames(r)))
      colnames(r) <- rep.int("", ncol(r))
    setN <- function(i, nams) colnames(r)[i] <- if (is.null(nams))
      ""
    else nams
    if (nn1)
      setN(1, N1)
    if (nn2)
      setN(1 + l1, N2)
    if (fix.na)
      setN(ncol(r), Nna)
  }
  r
}

################################ rbind.na ###################################################

#' qpcR rbind.na method.
#'
#' @param ... (generalized) vectors or matrices. See base::rbind
#' @param deparse.level integer controlling the construction
#'   of labels in the case of non-matrix-like arguments. See base::rbind
#'
#' @return a matrix combining the ... arguments row-wise.
#'
#' @export
#'
#' @examples
#' ### Vectors:
#' a <- c(rep(1, 5), NA, seq(1:5))
#' b <- c(rep(1, 4), NA, seq(1:7))
#'
#' # Complete shorter vector with NAs:
#' rbindna(a,b)
#'
rbindna <- function (..., deparse.level = 1)
{
  na <- nargs() - (!missing(deparse.level))
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na == 0)
    return(NULL)
  if (na == 1) {
    if (isS4(..1))
      return(methods::rbind2(..1))
    else return(matrix(..., nrow = 1))
  }
  if (deparse.level) {
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level ==
            2)
          deparse(r)
      }
      else r
    }
  }
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, ncol))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    if (deparse.level) {
      if (fix.na)
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl)))
        iV <- iV & (nmi == "")
      ii <- if (fix.na)
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
          names(argl)[i] <- nmi
      }
    }
    nCol <- as.numeric(sapply(argl, function(x) if (is.null(ncol(x))) length(x) else ncol(x)))
    maxCol <- max(nCol, na.rm = TRUE)
    argl <- lapply(argl, function(x) if (is.null(ncol(x)))
      c(x, rep(NA, maxCol - length(x)))
      else cbind(x, matrix(, nrow(x), maxCol - ncol(x))))
    namesVEC <- rep(NA, maxCol)
    for (i in 1:length(argl)) {
      CN <- colnames(argl[[i]])
      m <- !(CN %in% namesVEC)
      namesVEC[m] <- CN[m]
    }
    for (j in 1:length(argl)) {
      if (!is.null(ncol(argl[[j]])))
        colnames(argl[[j]]) <- namesVEC
    }
    r <- do.call(rbind, c(argl[-1L], list(deparse.level = deparse.level)))
  }
  d2 <- dim(r)
  colnames(r) <- colnames(argl[[1]])
  r <- methods::rbind2(argl[[1]], r)
  if (deparse.level == 0)
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2)
    return(r)
  Nrow <- function(x) {
    d <- dim(x)
    if (length(d) == 2L)
      d[1L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Nrow(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Nrow(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(rownames(r)))
      rownames(r) <- rep.int("", nrow(r))
    setN <- function(i, nams) rownames(r)[i] <- if (is.null(nams))
      ""
    else nams
    if (nn1)
      setN(1, N1)
    if (nn2)
      setN(1 + l1, N2)
    if (fix.na)
      setN(nrow(r), Nna)
  }
  r
}

############################### Re-order columns ##########################################

#' Re-ordering columns in dataframes:
#'
#' @param data dataframe to be ordered
#' @param vars nammed vectors of new positions. See details.
#'
#' @details This function provides an easy way to re-order the columns of a
#'   dataframe. The "vars" parameter must be a nammed numeric vectors with the
#'   names corresponding to the targeted columns and the numbers corresponding
#'   to their desired new positions.
#'
#' @return The dataframe with the desired colum order.
#'
#' @export
#'
#' @examples
#'  #### Dummy data:
#'  dat <- data.frame("Fac1" = c(rep("A", 6), rep("B",6)),
#'                    "Var1" = rnorm(12, mean = 20, sd = 1),
#'                    "Fac2" = rep(c("Low","High","Low","High"),
#'                                 each=3),
#'                    "Var2" = c(rnorm(3,7), rnorm(3,9),
#'                               rnorm(3,12), rnorm(3,15)))
#'
#'  # factor columns at the begining.
#'  arrange.vars(dat, vars =c("Fac2" = 2))
#'  # factor columns at the end.
#'  arrange.vars(dat, vars =c("Fac1" = 3, "Fac2" = 4))
#'
arrange.vars <- function(data, vars){
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars

  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec) == var.nr )

  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}




