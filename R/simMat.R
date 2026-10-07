simMat <- function(data, method, diag = TRUE, upper = TRUE, na.rm = FALSE, verbosity = 2, plot = FALSE, ...) {

  # version 2.5 (7 Oct 2026)

  method <- match.arg(method, c("Baroni", "Jaccard", "Simpson", "Sorensen"))

  if (verbosity > 1)  start.time <- Sys.time()

  data <- as.data.frame(data)  # accommodates SpatRaster, tibble, etc.

  stopifnot(all(data >= 0 & data <= 1, na.rm = TRUE))

  n.subjects <- ncol(data)

  sim.mat <- matrix(nrow = n.subjects, ncol = n.subjects,
                    dimnames = list(colnames(data), colnames(data)))

  inds <- triMatInd(sim.mat, lower = TRUE, list = TRUE)

  n.pairs <- length(combn(n.subjects, m = 2)) / 2

  if (verbosity > 1)  message("Computing ", n.pairs, " pair-wise similarities...")

  if (verbosity > 0) {
    progbar <- txtProgressBar(min = 0, max = n.pairs, style = 3, char = "-")
  }

  if (!anyNA(data))  na.rm <- FALSE  # otherwise unnecessary slower pairwise checks

  pair <- 0
  for (ind in inds) {

    pair <- pair + 1
    if (verbosity > 0) setTxtProgressBar(progbar, pair)

    row <- rownames(sim.mat)[ind[1]]
    col <- colnames(sim.mat)[ind[2]]

    x <- data[ , row]
    y <- data[ , col]

    if (na.rm && anyNA(c(x, y))) {
      # here because slower if checked repeatedly by fuzSim()
      finite <- is.finite(x) & is.finite(y)
      x <- x[finite]
      y <- y[finite]
    }  # end if na.rm

    sim.mat[ind[1], ind[2]] <- fuzSim(x, y, method = method, simplif = TRUE)
  }  # end for ind

  if (diag) diag(sim.mat) <- 1
  if (upper) {
    # https://stat.ethz.ch/pipermail/r-help/2008-September/174475.html
    ind <- upper.tri(sim.mat)
    sim.mat[ind] <- t(sim.mat)[ind]
  }

  if (isTRUE(plot)) {
    graphics::image(x = 1:ncol(sim.mat), y = 1:nrow(sim.mat), z = sim.mat,
                    axes = FALSE, xlab = "", ylab = "", ...)
    axis(side = 1, at = 1:ncol(sim.mat), tick = FALSE, labels = colnames(sim.mat), las = 2, cex.axis = 0.6)
    axis(side = 2, at = 1:nrow(sim.mat), tick = FALSE, labels = rownames(sim.mat), las = 2, cex.axis = 0.6)
  }

  if (verbosity > 1) {
    message ("\nFinished!")
    timer(start.time)
  }
  return(sim.mat)
}
