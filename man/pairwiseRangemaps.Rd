\name{pairwiseRangemaps}
\alias{pairwiseRangemaps}
\title{
Pairwise intersection (and union) of range maps
}
\description{
This function takes a set of rangemaps and returns a matrix containing the areas (in square meters) of their pairwise intersections; optionally, also their individual areas and/or the areas of their pairwise unions.
}
\usage{
\special{
pairwiseRangemaps(rangemaps, projection = NULL, diag = TRUE, unions = TRUE,
region = NULL, buffer = FALSE, verbosity = 1, Ncpu = 1, Nchunks = 1,
subchunks = NULL, filename = paste0(rangemap_matrix_, format(Sys.time(),
"\%Y-\%m-\%d_\%Hh\%Mm\%Ss"), ".csv")
}
}
\arguments{
  \item{rangemaps}{
named list of \code{\link[terra:SpatVector]{terra::SpatVector}} polygon layers; or a character vector of file names (in which case they will be imported with \code{\link[terra:vect]{terra::vect}}), including the extension (e.g. ".shp" or ".gpkg") and the folder paths if not in the working directory.
}
  \item{projection}{
DEPRECATED argument, previously required by function 'PBSmapping::importShapefile', which is now replaced with \code{\link[terra:vect]{terra::vect}}. Ignored with a message if provided. Note that area computations are now more accurate after the implementation of geodesic algorithms, which account for the curvature of the Earth (see \code{\link[terra:expanse]{terra::expanse}}).
}
  \item{diag}{
logical, whether to fill the diagonal of the resulting matrix with the individual area of each rangemap. The default is TRUE, and it is also automatically set to TRUE (as it is necessary) if unions = TRUE.
}
  \item{unions}{
logical (default TRUE), whether to fill the upper triangle of the resulting matrix with the area of union of each pair of rangemaps. It is not as computationally intensive as the intersection, as it is calculated not with spatial but with algebraic operations within the matrix (union = area1 + area2 - intersection).
}
  \item{region}{
optional \code{\link[terra:SpatVector]{terra::SpatVector}} polygon layer delimiting the study region. If provided, 'rangemaps' will be clipped with it using \code{\link[terra:crop]{terra::crop}} before computing areas.
}
  \item{buffer}{
optional logical value to pass to \code{\link[terra:makeValid]{terra::makeValid}}. The default is FALSE; try setting it to TRUE if you get topological errors, but see the \code{\link[terra:makeValid]{terra::makeValid}} help file for caveats!
}
  \item{verbosity}{
integer number (default 1, possible values 0 to 2) indicating the volume of progress messages to display.
}
  \item{Ncpu}{
integer indicating the number of CPUs (central processing units) to employ if parallel computing is to be used. The default is 1, which means no parallel computing, and this is currently the recommended value, as parallelizing doesn't necessarily speed things up and spatial operations are now substantially more efficient than when this argument was introduced. You might want to increase 'Ncpu' if you have a very large number of rangemaps, the computation is taking way too long, and your machine has several cores that can be used simultaneously.
}
  \item{Nchunks}{
either an integer indicating the number of chunks of rows in which to divide the results matrix for calculations, or character "decreasing" to indicate that the matrix should be divided into chunks of decreasing number of rows (as intersections are calculated in the lower triangle, rows further down the matrix have an increasing number of intersections to compute). Note, however, that rangemap size, not rangemap number, is the main determinant of computation time. The default is 1 (no division of the matrix), which is currently the recommended value, as spatial operations are now substantially more efficient than when this argument was introduced, but you might want to increase this if you have a very large number of rangemaps and your machine is getting clogged. With chunks, each set of rows of the matrix is computed and saved to disk, and the memory is cleaned before the next chunk begins.
}
  \item{subchunks}{
  optional integer vector specifying which chunks to actually compute. This is useful if a previous, time-consuming run of pairwiseRangemaps() with chunks was interrupted (e.g. by a power outage) and you want to compute only the remaining chunks.
}
  \item{filename}{
optional character vector indicating the name of the .csv file to save the resulting matrix to. The default is e.g. "rangemap_matrix_2026-07-23_09h21m05s.csv", but with the user's current timestamp. Can be NULL for no file export.
}
}
\value{
This function returns a square matrix containing, in the lower triangle, the areas (in square meters) of the pairwise intersections among the input 'rangemaps'; in the diagonal (if diag = TRUE or union = TRUE), the area of each rangemap; and in the upper triangle (if union = TRUE), the areas of the pairwise unions among the rangemaps.
}
\references{
Barbosa A.M. & Estrada A. (2016) Calcular corotipos sin dividir el territorio en OGUs: una adaptacion de los indices de similitud para su utilizacion directa sobre areas de distribucion. In: Gomez Zotano J., Arias Garcia J., Olmedo Cobo J.A. & Serrano Montes J.L. (eds.), Avances en Biogeografia. Areas de Distribucion: Entre Puentes y Barreras, pp. 157-163. Editorial Universidad de Granada & Tundra Ediciones, Granada (Spain)
}
\author{
A. Marcia Barbosa
}
\note{
This function previously used the \pkg{PBSmapping} package to import and intersect the rangemaps and to calculate areas. Now it uses the \pkg{terra} package instead and it is considerably faster, following recent developments therein. Mind that, following the implementation of geodesic algorithms, area computations are more accurate with unprojected maps, which are now used by default by \code{\link[terra:expanse]{terra::expanse}}. Some differences (for better) can thus arise between the results of the previous version and the current version (from \pkg{fuzzySim} 4.9.4).
}
\seealso{
\code{\link{rangemapSim}}
}
