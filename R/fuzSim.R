fuzSim <- function(x, y, method, na.rm = TRUE, simplif = FALSE) {

  # version 2.1 (7 Oct 2026)

  if (!simplif) {
    if (inherits(x, "SpatRaster"))
      x <- terra::values(x, mat = FALSE, dataframe = FALSE)
    if (inherits(y, "SpatRaster"))
      y <- terra::values(y, mat = FALSE, dataframe = FALSE)

    # for non-vector inputs:
    x <- unlist(x)
    y <- unlist(y)

    stopifnot(length(x) == length(y),
              # min(c(x, y, na.rm = TRUE)) >= 0,
              # max(c(x, y, na.rm = TRUE)) <= 1,
              all(c(x, y) >= 0 & c(x, y) <= 1, na.rm = TRUE)
    )

    if (na.rm && anyNA(c(x, y))) {
      finite <- is.finite(x) & is.finite(y)
      x <- x[finite]
      y <- y[finite]
    }

    method <- match.arg(method, c("Baroni", "Jaccard", "Simpson", "Sorensen"))
  }

  dab.methods <- c("Baroni")
  A <- sum(x)
  B <- sum(y)
  C <- sum(pmin(x, y))
  if (method %in% dab.methods) D <- sum(1 - pmax(x, y))

  if (method == "Baroni") return((sqrt(C * D) + C) / (sqrt(C * D) + A + B - C))
  if (method == "Jaccard") return(C / (A + B - C))
  if (method == "Simpson") return(C / min(A, B))
  if (method == "Sorensen") return(2 * C / (A + B))
}
