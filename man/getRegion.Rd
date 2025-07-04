\name{getRegion}
\alias{getRegion}
\title{Get region}
\usage{
getRegion(
  pres.coords,
  type = "width",
  clust_dist = 100,
  dist_mult = 1,
  width_mult = 0.5,
  weight = FALSE,
  CRS = NULL,
  dist_mat = NULL,
  dist_method = "auto",
  verbosity = 2,
  plot = TRUE
)
}
\arguments{
\item{pres.coords}{a SpatVector of points, or an object inheriting class 'data.frame' with 2 columns containing, respectively, the x and y, or longitude and latitude coordinates (in this order!) of the points where species presence was recorded.
}
\item{type}{character indicating which procedure to use for defining the region around 'pres.coords'. Options are:
\itemize{
\item "width": a buffer whose radius is the minimum diameter of the 'pres.coords' spatial extent (computed with \code{\link[terra:width]{terra::width()}}), multiplied by 'width_mult';
\item "mean_dist": a buffer whose radius is the mean pairwise \code{\link[terra:distance]{terra::distance()}} among 'pres.coords', multiplied by 'dist_mult';
\item "inv_dist": a buffer whose radius is inversely proportional to the sum of the distances from each point to all other points in 'pres.coords' (a rough measure of how isolated each point is, possibly indicating an opportunistic record in a sparsely surveyed area);
\item "clust_mean_dist": a different buffer around each cluster of 'pres.coords' (clusters computed with \code{\link[stats:hclust]{stats::hclust()}}, method = "simple") and then \code{\link[stats:cutree]{stats::cutree()}} with h = clust_dist*1000), sized according to the mean pairwise distance of each cluster's 'pres.coords'.
\item "clust_width": a different buffer around each cluster of 'pres.coords' (clusters computed as described for 'clust_mean_dist'), sized according to the \code{\link[terra:width]{terra::width()}} of each cluster's 'pres.coords'.
}
}
\item{clust_dist}{if 'type' involves clusters, numeric value specifying the distance threshold (in km) within which points are clustered together. Default 100.}
\item{dist_mult}{if 'type' involves distance, multiplier of the mean pairwise point distance to use for the \code{\link[terra:buffer]{terra::buffer()}} radius around each cluster. Default 1.}
\item{width_mult}{if 'type' involves width, multiplier of the width to use for the \code{\link[terra:buffer]{terra::buffer()}} radius. Default 0.5.}
\item{weight}{logical (used only if 'type' includes clusters) indicating whether to weigh the radius of the buffer around each cluster proportionally to the number of points that it includes. Default FALSE; if set to TRUE, clusters with fewer points (possibly indicating more sparsely surveyed areas) get proportionally smaller buffers than the mean distances among them.}
\item{CRS}{coordinate reference system of 'pres.coords' (if it is not a SpatVector with a defined CRS already), in one of the following formats: WKT/WKT2, <authority>:<code>, or PROJ-string notation (see \code{\link[terra:crs]{terra::crs()}}).}
\item{dist_mat}{optional matrix of pairwise distances among 'pres.coords', to use (if 'type' includes a string "dist" or "clust") for efficiency instead of computing a new one. Should normally be computed with \code{\link[terra:distance]{terra::distance()}}, \code{\link[geodist:geodist]{geodist::geodist()}}, or another function that takes the Earth's curvature into account.}
\item{dist_method}{argument to pass to [distMat()] (if 'dist_mat' is NULL) specifying the method for distance calculation. The default is "auto", or "haversine" if 'type' includes the string "clust", to avoid different clusters generating a different automatic method selection.}
\item{verbosity}{integer indicating the amount of messages to display along the process. The default is 2, for all available messages.}
\item{plot}{logical (default TRUE) indicating whether to plot the resulting region (in yelow), together with the input 'pres.coords' (black points, or points coloured according to their cluster) and a label with the number of points in each cluster (if 'type' involves clusters).}
}
\value{
SpatVector polygon delimiting a region around 'pres.coords'
}
\description{
This function computes a polygon around a set of point coordinates under given criteria, which may be useful for delimiting background or (pseudo)absence regions for computing species distibution models. Some of the 'type' options, especially those involving clusters or inverse distance, attempt to address survey bias by making smaller polygons around areas with fewer or more isolated points.
}
\details{
Most methods for computing species distribution models require predictor values for regions beyond those with species occurrence records, i.e. background or (pseudo)absence areas. The extent (as well as the spatial resolution) of these regions has a strong effect on model predictions. Ideally, they should include the areas that are within the reach of the species AND were reasonably surveyed (though you can further refine the latter with \code{\link{selectAbsences}} and an optional \code{\link{biasLayer}}). While sometimes we have a large enough and delimited area that we can use (e.g. when modelling a region where a national or regional distribution atlas is available), often we need to approximate the areas that appear to be both reasonably surveyed and within the species' reach.

Mind that no automated procedure can properly address all possible issues related to uneven data collection, or properly conform to all possible species distribution and survey patterns. Mind also that the output region from this function does not consider geographical barriers, or other factors that should also be taken into account when delimiting a region for modelling.

It is thus recommended to try different values for 'type' and associated parameters; judge for yourself which one provides the most plausible approximation to the surveyed region accessible to your target species; and possibly post-process (i.e. further edit) the resulting region in light of the available knowledge of that species' distribution, survey patterns and study region.
}
\seealso{
\code{\link[terra:buffer]{terra::buffer()}}, \code{\link[terra:width]{terra::width()}}, \code{\link[terra:crop]{terra::crop()}}
}
\examples{
\dontrun{
# you can run these examples if you have 'terra' and 'geodata' installed

# download example data:

occs <- geodata::sp_occurrence("Triturus", "pygmaeus")

occs_sv <- terra::vect(occs, geom = c("lon", "lat"), crs = "EPSG:4326")

cntry <- geodata::world(path = tempdir())


terra::plot(occs_sv)

terra::plot(cntry, lwd = 0.2, add = TRUE)


# compute regions with some different methods:

reg1 <- fuzzySim::getRegion(occs_sv)

terra::plot(cntry, lwd = 0.2, add = TRUE)


reg2 <- fuzzySim::getRegion(occs_sv, type = "inv_dist")

terra::plot(cntry, lwd = 0.2, add = TRUE)

terra::plot(reg2, lwd = 4, border = "orange", add = TRUE)


reg3 <- fuzzySim::getRegion(occs_sv, type = "clust_width", weight = TRUE,
width_mult = 0.3)

terra::plot(cntry, lwd = 0.2, add = TRUE)

terra::plot(reg3, lwd = 4, border = "orange", add = TRUE)


# note it is up to the user to pre-process the data (e.g. by removing erroneous
# records) and/or post-process the region (e.g. by erasing islands, countries,
# or continents that are not accessible to the target species).
}
}
\author{
A. Marcia Barbosa
}
