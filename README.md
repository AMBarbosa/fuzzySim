# `fuzzySim`: Fuzzy Similarity in Species Distributions

<p align="center">
<img width="180" height="180" alt="image" title="logo composed with terra, geovizr, magick and hexSticker" src="https://github.com/user-attachments/assets/eec03c5e-80de-4048-8dbd-bdd6235062c4" />
</p>

`fuzzySim` can convert binary presence-absence to fuzzy occurrence data, using e.g. trend surface analysis, inverse distance interpolation or prevalence-independent environmental favourability modelling, for one or more species simultaneously.

It can then calculate e.g. fuzzy change and fuzzy similarity among (fuzzy) species distributions and/or among (fuzzy) regional species compositions, avoiding the use of thresholds and instead using fuzzy logic versions of known similarity indices.


# Install the package:

Package `fuzzySim` can be installed from CRAN, where it is updated relatively often:

`install.packages("fuzzySim")`

The development version may include more recent enhancements and small bug fixes. `fuzzySim` was created long before I knew GitHub, and it's still normally maintained on its original platform, [R-Forge](https://fuzzySim.r-forge.r-project.org/). So, you can normally install the latest version from there:

`install.packages("fuzzySim", repos="http://R-Forge.R-project.org")`

Sadly, [R-Forge](https://fuzzySim.r-forge.r-project.org/) servers are sometimes offline, so I try to also keep an updated version of the development version here on GitHub:

`remotes::install_github("AMBarbosa/fuzzySim")`


# Package homepage:

Find out more about the package, along with usage manuals and literature examples, at https://fuzzysim.r-forge.r-project.org/
