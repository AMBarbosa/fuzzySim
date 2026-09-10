clampVars <- function(vars,
                      ref,
                      # var.cols = NULL,
                      values = TRUE
) {
  # version 1.0 (4 Ago 2026)

  if (inherits(ref, "SpatVector")) {
    if (inherits(vars, "SpatRaster")) {
      if (nzchar(terra::crs(vars)) && nzchar(terra::crs(ref)) && !terra::same.crs(vars, ref)) {
        message("Projecting 'ref' to the same CRS as 'vars'")
        ref <- terra::project(ref, vars)
      }
      # } else {
      # if (ncol(ref) != 2)
      #   stop("If not a spatial object, 'ref' must have only two columns\nwith the X (longitude) and Y (latitude) coordinates, respectively.")
      # ref <- terra::vect(ref, geom = colnames(ref), crs = terra::crs(vars))
    } else stop("If 'ref' is a SpatVector, 'vars' must be a SpatRaster")

  } else {  # if 'ref' !SpatVector
    if (inherits(ref, "matrix") || inherits(ref, "data.frame")) {
      ref <- as.data.frame(ref)  # accommodates tibbles etc.
    } else if (inherits(ref, "maxnet")) {
        ref <- rbind(ref$varmin, ref$varmax)
    } else ref <- modEvA::mod2obspred(ref, x.only = TRUE)
  }

  if (inherits(vars, "SpatRaster")) {
    if (inherits(ref, "SpatVector"))
      ref <- terra::mask(vars, ref)
    # ref <- terra::extract(vars, ref)

    # if (!is.null(var.cols)) {
    #   vars <- vars[[var.cols]]
    #   if (inherits(ref, "SpatRaster")) ref <- ref[[var.cols]]
    #   if (inherits(ref, "data.frame")) ref <- ref[, var.cols]
    # }

    ranges <- apply(as.data.frame(ref), 2, range, na.rm = TRUE)

    for (v in colnames(ranges)) {
      vars[[v]] <- terra::clamp(vars[[v]],
                                lower = min(ranges[ , v]),
                                upper = max(ranges[ , v]),
                                values = values)
    }
    return(vars)

  } else {

    if (inherits(vars, "data.frame")) {
      vars <- as.data.frame(vars)  # accommodates matrices, tibbles, etc.
      ref <- as.data.frame(ref)

      # if (!is.null(var.cols)) {
      #   vars <- vars[ , var.cols]
      #   ref <- ref[ , var.cols]
      # }

      ranges <- apply(ref, 2, range, na.rm = TRUE)

      for (v in colnames(ref)) {
        vars[ , v] <- terra::clamp(vars[ , v],
                                   lower = ranges[1, v],
                                   upper = ranges[2, v],
                                   values = values)
      }

      return(vars)

    } else {

      stop("'vars' must be either a SpatRaster or an object inheriting class 'data.frame'")
    }
  }
}
