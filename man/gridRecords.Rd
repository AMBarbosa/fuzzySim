\name{gridRecords}
\alias{gridRecords}
\title{
Get species occurrence records gridded onto a raster
}
\description{
This function takes a raster stack and a set of geographic coordinates of a species' presence (and optionally absence) records, and returns a table with the presences and absences and the values of the rasters in the grid of pixels. If absence coordinates are not supplied, all pixels without any presence point will be returned as absences.
}
\usage{
gridRecords(rst, pres.coords, abs.coords = NULL)
}
\arguments{
  \item{rst}{
a Raster* object with the desired spatial resolution and extent for the species presence-absence data. The raster should be masked so that pixels have values only in relevant and reasonably surveyed areas.
}
  \item{pres.coords}{
a matrix or data frame with two columns containing, respectively, the longitude and latitude coordinates (in this order, and in the same coordinate reference system as 'rst') of the points where species presence was detected.
}
  \item{abs.coords}{
same as 'pres.coords' but for points where the species was not detected. If kept NULL, all pixels with a finite value that are not intersected by 'pres.coords' will be output as absence grid cells.
}
  \item{keep.n}{logical, whether to output the number of presence points that fall within each pixel, as opposed to 1 for every cell with at least one presence point. The default is FALSE (the latter).
}
}
\details{
}
\value{
This function returns a data frame with the following columns:
- 'presence': 1 for the cells with at least one presence point, and 0 for the cells with no presence points
- 'x', 'y': centroid coordinates of each pixel
- 'cellnumber': 
- one column for each variable in 'rst', with the value of each pixel for that variable.
}
\note{
}
\references{
}
\author{
A. Marcia Barbosa
}
\seealso{
'extract' in package \pkg{raster}
}
\examples{
}
