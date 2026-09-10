rangemapSim <- function(rangemap.matrix,
                        total.area = NULL,
                        method = c("Jaccard", "Sorensen", "Simpson", "Baroni"),
                        diag = TRUE,
                        upper = TRUE,
                        verbosity = 0) {
  # version 2.3 (27 Aug 2026)

  # rangemap.matrix: a matrix resulting from the pairwiseRangemaps function, with the areas of the range maps in the diagonal, their pair-wise intersections in the lower triangle, and their pair-wise unions in the upper triangle
  # total.area: total size of the study area, or a polygon delimiting it

  method <- match.arg(method, c("Jaccard", "Sorensen", "Simpson", "Baroni"))

  if (!is.null(total.area)){
    if (inherits(total.area, "SpatVector") && terra::is.polygons(total.area))
      total.area <- sum(terra::expanse(total.area))
    if (!(is.numeric(total.area) && length(total.area) == 1))
      stop("'total.area' must be either a polygon 'SpatVector' or a numeric value")
  } else {
    if (method == "Baroni") stop ("'Baroni' method requires 'total.area' not NULL")
  }

  n.rangemaps <- ncol(rangemap.matrix)
  rangemap.names <- colnames(rangemap.matrix)
  sim.matrix <- matrix(nrow = n.rangemaps, ncol = n.rangemaps, dimnames = list(rangemap.names, rangemap.names))
  lower.inds <- triMatInd(rangemap.matrix, lower = TRUE, list = TRUE)

  for (i in lower.inds) {
    area1 <- rangemap.matrix[i[1], i[1]]
    area2 <- rangemap.matrix[i[2], i[2]]
    int <- rangemap.matrix[i[1], i[2]]
    uni <- rangemap.matrix[i[2], i[1]]
    sim.matrix[i[1], i[2]] <- simFromSetOps(size1 = area1, size2 = area2, intersection = int, union = uni, total.size = total.area, method = method, verbosity = verbosity)
  }; rm(i)

  if (diag) diag(sim.matrix) <- 1
  if (upper) {  # https://stat.ethz.ch/pipermail/r-help/2008-September/174475.html
    up <- upper.tri(sim.matrix)
    sim.matrix[up] <- t(sim.matrix)[up]
  }

  return(sim.matrix)
}
