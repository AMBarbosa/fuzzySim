\name{fuzzyRangeChange}
\alias{fuzzyRangeChange}
\title{
Range change based on continuous (fuzzy) values
}
\description{
This function quantifies overall range change (expansion, contraction, maintenance and balance) based on either presence-absence data or on the continuous predictions (presence probability, favourability, or suitability) of two species distribution / ecological niche models.
}
\usage{
fuzzyRangeChange(pred1, pred2, number = TRUE, prop = TRUE,
na.rm = TRUE, round.digits = 2,
measures = c("Gain", "Loss", "Stable positive", "Stable negative", "Balance"),
plot = TRUE, plot.type = "waterfall", x.lab = TRUE, ...)
}
\details{
No thresholds are involved in any of the calculations. See below for a simple illustration of how each metric works, using the special case when all favourabilities are either 0 or 1.

  \preformatted{
present <- c(0, 0, 1, 0, 1)
future  <- c(1, 1, 1, 0, 0)

fuzzyRangeChange(present, future)

#                 Number  Prop
# Gain                 2  1.00  # 2 new presences, 100\% gained presence area
                                # (new presences double the stable presences)
# Loss                -1 -0.50  # 1 (50\%) initial presences became absences
# Stable positive      1  0.50  # 1 (50\%) initial presences remained presences
# Stable negative      1  0.33  # 1/3 of the initial absences remained absences
# Balance              1  0.50  # 1 more total presence (50\% more) than initial


present <- c(0, 0, 1, 0, 1, 1)
future  <- c(1, 1, 1, 0, 0, 1)

fuzzyRangeChange(present, future)

#                 Number  Prop
# Gain                 2  0.67  # 2/3 of initial absences became presences
# Loss                -1 -0.33  # 1/3 of initial presences became absences
# Stable positive      2  0.67  # 2/3 of initial presences remained presences
# Stable negative      1  0.33  # 1/3 of initial absences remained absences
# Balance              1  0.33  # 1 more presence than initially, 33\% gain in
                                # total presence area (from 3 to 4 presences)
  }

With fuzzy input values (e.g. predicted favourabilities) instead of 1s and 0s, the mechanism is the same, but a fuzzy value (e.g. 0.35, 0.92) counts as a "partial" presence, and proportions are calculated accordingly - e.g., a 25\% increase is from 0.2 to 0.25, but from 0.5 to 0.625 and from 0.8 to 1. So, higher favourability values have more weight in the output, but small favourability values also chip in, as they also represent a bit of (potential) presence.
}
\arguments{
  \item{pred1}{numeric vector or SpatRaster layer containing the values (between 0 and 1) that will serve as reference.}
  \item{pred2}{numeric vector or SpatRaster layer containing the values (between 0 and 1) whose change will be computed relative to \code{pred1}. Must be of the same dimensions and in the same order.}
  \item{number}{logical value (default TRUE) indicating if results should include the (fuzzy) number of cases.}
  \item{prop}{logical value (default TRUE) indicating if results should include the proportion of the total number of cases.}
  \item{na.rm}{logical value (default TRUE) indicating whether NA values should be ignored.}
  \item{round.digits}{argument to pass to \code{\link{fuzzyOverlay}}, indicating the number of decimal places to which to round 'pred' for computing maintenance, i.e. stable positives and stable negatives. The default is 2.}
  \item{measures}{character vector listing the range change measures to compute. The default includes all available measures.}
  \item{plot}{logical value (default TRUE) indicating whether to plot the results.}
  \item{plot.type}{(if plot=TRUE) character value indicating the type of plot to produce. Can be "barplot" or "waterfall" (synonyms, default) or "lollipop".}
  \item{x.lab}{logical value indicating whether to add the x axis labels to the plot (i.e., the names below each lollipop or bar). The default is TRUE, but users may set it to FALSE and then add labels differently (e.g. with different names or rotations).}
  \item{\dots}{(if plot=TRUE) additional arguments to pass to \code{\link{barplot}} (if 'plot.type' is "barplot" or "waterfall") or to \code{\link[modEvA:lollipop]{modEvA::lollipop}} (if 'plot.type' is "lollipop").}
}
\value{
This function returns a data frame with the following values in different rows (among those included in 'measures'):

\item{Gain}{sum of the values that increase from 'pred1' to 'pred2' (fuzzy equivalent of the number of localities that gain presence)}
\item{Loss}{sum of the values that decrease from 'pred1' to 'pred2' (fuzzy equivalent of the number of localities that lose presence)}
\item{Stable positive}{number of localities that remain favourable from 'pred1' to 'pred2' (fuzzy equivalent of the number of presence localities that remain as such)}
\item{Stable negative}{number of localities that remain unfavourable from 'pred1' to 'pred2' (fuzzy equivalent of the number of absence localities that remain as such)}
\item{Balance}{sum of the change in values from 'pred1' to 'pred2' (fuzzy equivalent of the overall balance between gained and lost presences, or the proportional gain in potential presence area)}

If number=TRUE (default), the output has a column named "Number" with the number of localities in each of the above categories. If prop=TRUE (also default), the output has a column named "Proportion" where this number is divided by the total number of reference values (i.e., the fuzzy range or fuzzy non-range size). If plot=TRUE (default), a plot is also produced from the values in the last column of the output.
}
\author{A. Marcia Barbosa}
\seealso{
\code{\link{fuzzyOverlay}} for row-wise or pixel-wise fuzzy comparisons; \code{\link{fuzSim}}; \code{\link{modOverlap}} for other ways to compare model predictions}
\examples{
# get an environmental favourability model for a rotifer species:

data(rotif.env)

names(rotif.env)

fav_current <- multGLM(rotif.env, sp.cols = 18, var.cols = 5:17,
step = TRUE, FDR = TRUE, trim = TRUE, P = FALSE, Fav = TRUE) $
predictions


# imagine you have a model prediction for this species in a future time
# (here we will create one by randomly jittering the current predictions)

fav_imag <- jitter(fav_current, amount = 0.2)
fav_imag[fav_imag < 0] <- 0
fav_imag[fav_imag > 1] <- 1


# compute range change between current and imaginary future predictions:

fuzzyRangeChange(fav_current, fav_imag)

fuzzyRangeChange(fav_current, fav_imag, las = 2)  # vertical labels

fuzzyRangeChange(fav_current, fav_imag, prop = FALSE)  # different table and y axis

fuzzyRangeChange(fav_current, fav_imag, ylim = c(-0.3, 0.3))  # different axis limits

fuzzyRangeChange(fav_current, fav_imag, plot.type = "lollipop")
}
