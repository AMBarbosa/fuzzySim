\name{clampVars}
\alias{clampVars}
\title{Clamp variables}

\description{
  This function is a wrapper (with additional input options) for \code{\link[terra:clamp]{terra::clamp}}, which caps each variable in a table or a (multilayer) raster to the extreme values observed in a reference (e.g. model training) dataset. Values below the reference minimum are set to that minimum (or to NA), and values above the reference maximum are set to that maximum (or NA). This may help prevent unreliable extrapolation when projecting a model beyond the training conditions.
}

\usage{
  clampVars(vars, ref, values = TRUE)
}

\arguments{
  \item{vars}{
    A \code{\link[terra:SpatRaster]{SpatRaster}} or an object inheriting class 'data.frame' containing the variables to be clamped, one in each layer or column, with names matching the corresponding variables in 'ref'.
  }
  \item{ref}{
    The reference (e.g. model training) dataset that defines the clamping limits for the (subset of) variables to be clamped in 'vars'. It can be:
      \itemize{
        \item a \code{\link[terra:SpatRaster]{SpatRaster}} or an object inheriting class \code{\link{data.frame}} containing, in separate layers or columns, the reference (e.g. model training) values for the variables in 'vars', or just two rows with the minimum and maximum value for each variable
        \item a \code{\link[terra:SpatVector]{SpatVector}} of points with the reference or training localities (both presences and pseudoabsences/background points), in which case 'vars' must be a \code{\link[terra:SpatRaster]{SpatRaster}} and should have values at all those localities
        \item a model object of class \code{\link[maxnet:maxnet]{maxnet}} or a class implemented in \code{\link[modEvA:mod2obspred]{modEvA::mod2obspred}} (except \code{\link[randomForest:randomForest]{randomForest}}, which does not store predictor values), from which the extreme training values of the predictor variables will be taken.
      }
  }
  \item{values}{logical value to pass to \code{\link[terra:clamp]{terra::clamp}}. If FALSE, 'vars' values outside the 'ref' range become NA; if TRUE (default), they get the extreme values found in 'ref'}
}

\details{
  If 'ref' is a \code{\link[terra:SpatVector]{SpatVector}} of points, it is used in a \code{\link[terra:mask]{terra::mask}} of the 'vars' SpatRaster, and values are then extracted from the raster predictors.

  Clamping of variables to stay within the range of model training values was originally implemented in Maxent (Phillips et al., 2006), but its rationale is applicable to other models. It can be recommended when projecting a model in space or time, to avoid extrapolating beyond the known responses.

  Future versions of this function might implement one-sided clamping (Pinilla-Buitrago et al., 2026).
}

\value{The function returns an object of the same general type as 'vars', i.e. a \code{\link[terra:SpatRaster]{SpatRaster}} or a \code{\link{data.frame}} of the clamped 'vars'.}

\author{A. Marcia Barbosa}

\seealso{
  \code{\link[terra:clamp]{terra::clamp}}, \code{\link{predict}}
}

\references{
Phillips S.J., Anderson R.P. & Schapire R.E. (2006) Maximum entropy modeling of species geographic distributions. Ecological Modelling 190, 231-259. https://doi.org/10.1016/j.ecolmodel.2005.03.026

Pinilla-Buitrago G.E., Kass J.M. & Anderson R.P. (2026) Extrapolation strategy matters when transferring ecological niche models: new visualization tools for informed decisions. Ecography n/a, e08590. https://doi.org/10.1002/ecog.08590
}

\examples{

\dontrun{
  # RASTER EXAMPLE:

  library(terra)

  elev <- terra::rast(system.file("ex/elev.tif", package = "terra"))
  slope <- terra::terrain(elev, v = "slope")
  vars <- c(elev, slope)
  terra::plot(vars)


  ## With SpatVector points as reference data:

  set.seed(1)
  pts <- terra::spatSample(vars, 100, na.rm = TRUE, as.points = TRUE)
  terra::plot(vars[[1]])
  terra::points(pts)

  vars_clamped <- clampVars(vars, pts)
  names(vars_clamped) <- paste0(names(vars_clamped), "_clamp")
  plot(c(vars, vars_clamped), axes = FALSE)


  ## With a data frame as reference data:

  terra::plot(vars)

  vals <- terra::extract(vars, pts)
  head(vals)

  vars_clamped <- clampVars(vars, vals)
  terra::plot(vars_clamped)
}

  # DATA FRAME EXAMPLE:

  cars
  ref <- data.frame(speed = c(6, 22), dist = c(12, 56))
  # you can use more rows, or fewer columns

  cars_clamp <- clampVars(cars, ref)
  cars_clamp


  # MODEL EXAMPLE:

  iris$Species <- ifelse(iris$Species == "virginica", NA,
                         ifelse(iris$Species == "setosa", 1, 0))
  mod <- glm(Species ~ Sepal.Length + Sepal.Width + Petal.Length,
             data = iris, family = binomial)

  iris_clamp <- clampVars(iris, mod)
  iris_clamp

  # try also with 'values = FALSE'
}
