\name{rangemapSim}
\alias{rangemapSim}
\title{
Pairwise similarity between rangemaps
}
\description{
Calculate pairwise similarity among rangemaps from a matrix of their areas of intersection and union
}
\usage{
rangemapSim(rangemap.matrix, total.area,
method = c("Jaccard", "Sorensen", "Simpson", "Baroni"),
diag = TRUE, upper = TRUE, verbosity = 0)
}
\arguments{
  \item{rangemap.matrix}{
a matrix like the one produced by function \code{\link{pairwiseRangemaps}}, containing the areas of pairwise intersection among rangemaps in the lower triangle, individual rangemap areas in the diagonal, and pairwise union areas in the upper diagonal.
}
  \item{total.area}{
[if 'method' uses shared absences (as is the case of "Baroni")] \code{\link[terra:SpatVector]{terra::SpatVector}} polygon delimiting the study region (same as the 'region' argument provided to \code{\link{pairwiseRangemaps}}); or a numeric value indicating the total size of the study area, in the same units as \code{rangemap.matrix}.
}
  \item{method}{
character value indicating the similarity index to use. Currently implemented indices are "Jaccard" (default), "Sorensen", "Simpson" and "Baroni".
}
  \item{diag}{
logical value (default TRUE) indicating whether the diagonal of the output matrix should be filled and not NA.
}
  \item{upper}{
logical value (default TRUE) indicating whether the upper triangle of the output matrix should be filled (symmetrical to the lower triangle) and not NA.
}
  \item{verbosity}{
integer number (default 0) indicating the volume of messages to display.
}
}
\details{
Distributional relationships among species are commonly determined based on pair-wise (dis)similarities in species' occurrence patterns. Some of the most commonly employed similarity indices are those of Jaccard (1901), Sorensen (1948), Simpson (1960) and Baroni-Urbani & Buser (1976), which are here implemented for comparing rangemaps based on their areas of intersection and union (Barbosa & Estrada, 2016).
}
\value{
This function returns a square matrix of pairwise similarities between the rangemaps in 'rangemap.matrix', calculated with the (first) similarity index specified in 'method'.
}
\references{
Barbosa A.M. & Estrada A. (2016) Calcular corotipos sin dividir el territorio en OGUs: una adaptacion de los indices de similitud para su utilizacion directa sobre areas de distribucion. In: Areas de distribucion: entre puentes y barreras. Universidad de Granada, Spain.

Baroni-Urbani C. & Buser M.W. (1976) Similarity of Binary Data. Systematic Zoology, 25: 251-259

Jaccard P. (1901) Etude comparative de la distribution florale dans une portion des Alpes et des Jura. Memoires de la Societe Vaudoise des Sciences Naturelles, 37: 547-579

Simpson G.G. (1960) Notes on the measurement of faunal resemblance. Amer. J. Sci. 258A, 300-311

Sorensen T. (1948) A method of establishing groups of equal amplitude in plant sociology based on similarity of species and its application to analyses of the vegetation on Danish commons. Kongelige Danske Videnskabernes Selskab, 5(4): 1-34
}
\author{
A. Marcia Barbosa
}
\seealso{
\code{\link{pairwiseRangemaps}}; \code{\link{simFromSetOps}}; \code{\link{simMat}}
}
