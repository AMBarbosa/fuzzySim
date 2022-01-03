# fuzzySim
Fuzzy Similarity in Species Distributions

`fuzzySim` can convert binary presence-absence to fuzzy occurrence data, using e.g. trend surface analysis, inverse distance interpolation or prevalence-independent environmental favourability modelling, for one or more species simultaneously.

It can then calculate e.g. fuzzy change and fuzzy similarity among (fuzzy) species distributions and/or among (fuzzy) regional species compositions, avoiding the use of thresholds and instead using fuzzy logic versions of known similarity indices such as Jaccard, Sørensen, Simpson, and Baroni-Urbani & Buser.

# install the latest version of `fuzzySim`:
install.packages("fuzzySim", repos="http://R-Forge.R-project.org")

`fuzzySim` is on CRAN, but the latest version (maintained on R-Forge) includes recent developments, such as new fuzzy metrics and a transition from 'raster' and 'sp' to the `terra` package, which is much faster.
