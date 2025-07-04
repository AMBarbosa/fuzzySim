# Version 4.34
### (Committed 2025-  - )

### New functions:

* distMat
- now used by distPres() and getRegion()


### Modified functions:

* getRegion
- added option dist_method="auto" (new default), as in distPres()
- renamed argument 'dist_measure' to 'dist_method', to match the terra::distance() new argument that it will call, as well as the corresponding argument to distPres()


* getRegion
- before geodist(), project to EPSG:4326 if !is.lonlat( perhaps=TRUE)


* getRegion
- use geodist::geodist instead of terra::distance (faster)
- added 'dist_measure' argument for geodist


* getRegion
- rows with no coordinates removed (with message) to avoid errors



# Version 4.33 
### (Committed 2025-03-21) -> CRAN


### Modified functions:
* Fav
    - added 'inv' argument to compute instead probability from favourability


### Other modified files:

* NEWS
    - reformatted and converted to .md to work with news()



# Version 4.32 
### (Committed 2025-02-24)


### New functions:

* biasLayer


### Modified functions:

* cleanCoords
    - fixed bug when input data is SpatVector
    - 'coord.cols' ignored with warning if provided and input is SpatVector



# Version 4.31 
### (Committed 2024-12-27)


### Modified functions:

* partialResp
    - fixed bug when only one variable in model
    - polygon instead of lines for confidence interval
    - new argument se.mult, default 1.96 for the 95% confidence interval, but can be changed e.g. to 1 (for just the SE) or 0 (for no CI)
    - initial plot with type="n" and points added if plot.points=TRUE, removing boilerplate code
    - added args 'reset.par' and '...' for plot()

* selectAbsences
    - moved "Absences not enough..." section upwards, avoiding unnecessary remaining computations

* distPres
    - removed double loop to increase efficiency
    - new default method "auto"
    - option to use faster methods like haversine and cosine distance if packageVersion("terra") >= "1.8.7"



# Version 4.30 
### (Committed 2024-12-16)


### New functions:

* partialResp



# Version 4.29 
### (Committed 2024-12-12) -> CRAN


### Modified functions:

* getRegion
    - added 'dist_mat' argument
    - changed default 'weight' to FALSE, as it looked too restrictive



# Version 4.28 
### (Committed 2024-12-05)


### Modified functions:

* cleanCoords
    - added 'extend' argument, default 0.1 to increase the output plot extent by 10% of the input coordinates range



# Version 4.27 
### (Committed 2024-11-22)


### Modified functions:

* fuzzyOverlay
    - fixed bug for raster inputs with NA pixels (reported by Jakub Nowosad: https://github.com/AMBarbosa/fuzzySim/issues/2#issue-2679366909)


### Other modified files:

* getRegion.Rd
    - added examples

* selectAbsences.Rd
    - added example with bias layer



# Version 4.26 
### (Committed 2024-10-29) -> CRAN


### Modified functions:

* gridRecords
    - output is SpatVector only if input is too

* selectAbsences
    - added data.frame() to fix error when input SpatVector



# Version 4.25 
### (Committed 2024-10-24)


### Modified functions:

* selectAbsences
    - 'bias' can now be a raster layer

* modelTrim
    - implement categorical variables in phylolm models



# Version 4.24 
### (Committed 2024-10-03)


### Modified functions:

* fuzSim, simMat, fuzzyOverlay, modOverlap, fuzzyRangeChange
    - inputs can be SpatRaster

* simMat
    - added 'plot' argument


### Other modified files:

* fuzzyOverlay.Rd
    - 'op' argument itemized list of options

* NAMESPACE
    - removed import of utils::installed.packages() [now using .packages(all.available = TRUE) instead]



# Version 4.23 
### (Committed 2024-09-27)


### Modified functions:

* getRegion
    - labels now appear for all clusters in plot

* gridRecords, selectAbsences
    - output plot SpatVector when 'terra' available



# Version 4.22 
### (Committed 2024-09-26)


### Modified functions:

* getRegion
    - hclust() now with method = "single"
    - cutree() with h = clust_dist (new argument)
    - fixed bug in 'weight', which wasn't correctly ordered (replaced merge with match)
    - 'weight' argument now also applies to type="width"



# Version 4.21 
### (Committed 2024-09-25)


### Modified functions:

* distPres, selectAbsences
    - add 'CRS' argument; if provided, terra::distance() is used

* distPres
    - add 'verbosity' argument

* getRegion
    - fix plotting error, now specifying 'terra::' before plot() and text()



# Version 4.20 
### (Committed 2024-09-24)


### New functions:

* getRegion



# Version 4.13 
### (Committed 2024-09-17)


### Modified functions:

* modelTrim
    - implemented for 'phylolm' models
    - added 'phy' and 'verbosity' arguments
    - changed 'if (class(model) %in% c("glm", "lm"))' to 'if (all(class(model) %in% c("glm", "lm")))', otherwise CRAN check error


### Other modified files:

* DESCRIPTION
    - added contributors



# Version 4.12 
### (Committed 2024-06-03)


### New functions:

* dms2dec



# Version 4.11.1 
### (Committed 2024-03-07)


### Modified functions:

* FDR, multGLM, corSelect
    - added argument 'test' (default "Chisq" for back-compatibility; future default probably "Rao" score as in SPSS)

* stepByStep
    - message and NA outputs if model has no variables



# Version 4.11 
### (Committed 2024-03-06)


### Modified functions:

* stepByStep
    - added option select="p.value" and associated arguments


### Other modified files:

* stepwise.Rd
    - changed 'sp.col' in Examples, for a model where a variable is included and then dropped



# Version 4.10.8 
### (Committed 2024-03-05)


### Modified functions:

* appendData
    - 'data1' and 'data2' coerced to data.frame (to work e.g. for SpatVectors)

* stepByStep
    - 'data' can now be a 'glm' model object, from which the order of the variables will be taken



# Version 4.10.7 
### (Committed 2024-01-24) -> CRAN


### Modified functions:

* selectAbsences
    - added 'dist.mat' argument for 'distPres'
    - updated help file to reflect this


### Other modified files:

* distPres.Rd
    - corrected typo "euclidian" -> "euclidean"
    - mention that (external) spherical distance is recommended for large extents
    - mention that 'method' is only used if 'dist.mat' is not provided



# Version 4.10.6 
### (Committed 2023-11-02)


### Modified functions:

* corSelect
    - fixed incorrect use of all.equal() to compare results between different 'select' criteria (thanks to bug report by Jose Carlos Guerrero)


### Other modified files:

* corSelect.Rd
    - mention {collinear} pkg under "See also"
    - mention that method 'pearson' is recommended for >30 rows, 'spearman' for >10 rows (as in ?collinear::collinear)

* fuzSim.Rd
    - mention fuzzy Jaccard = weighted Jaccard (Ioffe 2010), overlap, coincidence, consistence (Real et al. 2010)



# Version 4.10.5 
### (Committed 2023-10-04) -> CRAN


### Modified functions:

* corSelect
    - added 'cor' to 'select' criteria

* getPreds
    - replaced 'requireNamespace(raster)' with 'raster %in% rownames(installed.packages())', to avoid loading the pkg if unnecessary


==
# Version 4.10.4 
### (Committed 2023-09-28)
==

### Modified functions:

* cleanCoords
    - added 'year.min' and 'year.na.pass' arguments


==
# Version 4.10.3 
### (Committed 2023-09-25)
==

### Modified functions:

* appendData
    - added 'add.source' argument
    - x and y arguments renamed to data1 and data2
    - warning when no matching columns
    - out data frame created empty, instead of taking data1 values and then converting to NA

* gridRecords
    - replaced 'requireNamespace(raster)' with 'raster %in% rownames(installed.packages())', to avoid loading 'sp' and triggering all the deprecation warnings


### Other modified files:

* cleanCoords.Rd
    - now states that row names are carried from input dataframe


==
# Version 4.10.2 
### (Committed 2023-09-23)
==

### New functions:

* appendData


==
# Version 4.10.1 
### (Committed 2023-07-04)
==

### Modified functions:

* cleanCoords
    - added 'terra::' before plot() on SpatVector objects, to prevent newly appeared error


### Other modified files:

* Fav.Rd
    - now states that 'model' or 'pred' must have been obtained with weights=NULL

* corSelect.Rd
    - added "however" references against dropping correlated variables


==
# Version 4.10 
### (Committed 2023-05-24)
==

### Modified functions:

* multGLM
    - added 'coeff' argument to pass to 'corSelect', and changed 'cor.thresh' default accordingly
    - precluded error when no variables pass corSelect


==
# Version 4.9.13 
### (Committed 2023-05-22)
==

### Modified functions:

* corSelect
    - fixed bug that didn't allow correctly reporting the selected variables when coeff=FALSE


==
# Version 4.9.12 
### (Committed 2023-05-21)
==

### Modified functions:

* corSelect
    - added 'coeff' argument; if FALSE, high correlations are identified based on p-value (not coefficient) of correlation
    - 'p.value' column added to 'high.correlations' output table
    - fixed bug when there's >1 categorical variable to exclude a priori


==
# Version 4.9.11 
### (Committed 2023-05-17)
==

### Modified functions:

* multGLM
    - warning emitted if select="p.value" and trim=TRUE
    - started implementation of 'block.cols' argument (still not useable)


==
# Version 4.9.10 
### (Committed 2023-04-21)
==

### Modified functions:

* corSelect
    - 'select' default now VIF if 'sp.cols' is null
    - 'select' can now be null
    - message displayed saying what 'select' criterion was used
    - help file updated accordingly


==
# Version 4.9.9 
### (Committed 2023-04-14) -> CRAN
==

### Modified functions:

* cleanCoords
    - added arguments 'rm.abs' and 'abs.col'
    - when all records are removed, return 'data' with zero rows instead of error
    - input 'data' can also be a SpatVector of points, in which case output is too
    - 'coord.cols' default now NULL, can be taken from 'data' if class 'SpatVector'
    - 'coord.cols' can be NULL, in which case only the other removal criteria (absences, coordinate uncertainty) are applied

* gridRecords
    - input 'pres.coords' and 'abs.coords' can also be SpatVector of points


==
# Version 4.9.8 
### (Committed 2023-03-06) -> CRAN
==

### Modified files:

* CITATION
    - removed old-style personList(), and replaced citEntry() with bibentry(), as per new CRAN requirements



# Version 4.9.8 
### (Committed 2023-02-22)


### Modified functions:

* selectAbsences
    - added 'df' argument
    - output df now ordered by as.integer(rownames) rather than character rownames (1, 10, 100, 2, ...)
    - points now plotted in this order: unselected absences, selected absences, presences on top

* gridRecords
    - points now plotted in this order: absences, presences on top



# Version 4.9.7 
### (Committed 2023-01-25)


### Modified functions:

* cleanCoords
    - added 'plot' argument
    - changed default 'uncert.limit' from Inf to 50000



# Version 4.9.6 
### (Committed 2023-01-24)


### New functions:

* selectAbsences


### Modified functions:

* cleanCoords
    - fixed bug in 'rm.imprecise'



# Version 4.9.5 
### (Committed 2023-01-20)


### New functions:

* cleanCoords


### Modified functions:

* fuzzyRangeChange
    - added argument 'plot.type'
    - default 'plot.type' now "lollipop"

* entropy
    - added argument 'plot.type'
    - default 'plot.type' now "lollipop"



# Version 4.9.4 
### (Committed 2023-01-12)


### Modified functions:

* pairwiseRangemaps
    - replaced 'PBSmapping' (and consequently 'sp' and 'maptools') with 'terra' (future-proof and much faster)


### Other modified files:

* DESCRIPTION
    - packages 'PBSmapping', 'sp' and 'maptools' removed from 'Suggests'



# Version 4.9.3 
### (Committed 2023-01-09)


### Modified functions:

* gridRecords
    - added 'plot' argument

* corSelect
    - non-numeric variables excluded with warning (instead of causing error)

* multGLM
    - corrected bug that miscounted the variables excluded by FDR when FDR.first=FALSE (affecting only console messages)
    - total N variables now reported when TSA=TRUE (i.e. "... 46 with the spatial trend variable" instead of just "...plus the spatial trend variable")

* modelTrim
    - added family <- family(model) when method = "summary", as per bug report by J.C.Guerrero email 30/12/2022



# Version 4.9.2 
### (Committed 2022-10-27)


### Modified functions:

* multGLM
    - added 'FDR.first' argument


### Other modified files:

* multGLM.Rd
    - added recommendation to set trim=FALSE when select="p.value"
    - added JCG funding acknowledgment to Note

* fuzSim.Rd
    - help file now mentions Jaccard and Sorensen also as recommended metrics for model evaluation



# Version 4.9.1 
### (Committed 2022-10-13)


### Modified functions:

* corSelect
    - added 'verbosity' argument
    - message emitted saying how many vars selected and excluded (like in 'FDR') if verbosity > 0
    - message emitted saying which vars selected and excluded if verbosity > 1

* FDR
    - argument 'verbose' replaced with 'verbosity' (with warnings)
    - message emitted saying which vars selected and excluded if verbosity > 1

* multGLM
    - argument 'verbose' for 'FDR' replaced with 'verbosity' (as per 'FDR' change above)
    - argument 'select' now accepts additional criterion "p.value", in which case 'stepwise' used instead of 'step'
    - added arguments test.in, test.out, p.in and p.out to pass to 'stepwise'
    - trim=TRUE reverts back to using 'modelTrim' (rather than 'stepwise' with direcction='backward' as in previous experimental versions)
    - argument trim.fun (added in previous experimental versions) removed



# Version 4.9 
### (Committed 2022-10-11)


### Modified functions:

* summaryWald
    - fixed bug when model has no variables

* stepwise
    - break loop when a variable simultanously meets the criteria for inclusion and exclusion


### Other modified files:

* stepwise
    - help file with new references (two against AIC, and one more against stepwise)



# Version 4.8.1 
### (Committed 2022-08-01)


### Modified functions:

* gridRecords
    - moved if(na.rm) uppward to go back to using 5:ncol(result) instead of names(result), which was causing problems with some raster names (bug report by Alba Estrada)
    - ordered result by 'cell' and removed (confusing) row names
    - help file Examples now includes mapping multi-species gridded presences


### Other modified files:

* multGLM
    - help file notes that 'TSA' uses "type="Y" and is included in FDR, corSelect, etc.



# Version 4.8 
### (Committed 2022-07-21)


### Modified functions:

* gridRecords
    - warning if 'species' (not character or) has leading/trailing spaces
    - added check.names=FALSE to data.frame call, so that output column names are exactly as the unique values of 'species' (even if they are numbers), and the first species name doesn't come out different from the rest
    - added last argument in abs.coords <- terra::crds(rst, df=TRUE, na.rm=FALSE), otherwise some rows with missing data in some variables would be missing from the result (bug report by Alba Estrada)


### Other modified files:

* multGLM
    - help file notes that 'modelTrim' and 'stepwise' have different default significance thresholds, to explain why 'stepwise' may leave more variables in the model



# Version 4.7 
### (Committed 2022-07-07)


### Modified functions:

* multGLM
    - fixed message about N variables excluded by 'modelTrim' (now by trim.fun)
    - set trace=0 for 'stepwise'

* fuzzyConsensus
    - changed default biplot to FALSE, as it makes function much slower

* fuzzyRangeChange
    - changed "Stable presence" and "Stable absence" to "Stable positive" and "Stable negative" to avoid erroneous interpretations
    - added y axis label



# Version 4.6 
### (Committed 2022-07-05)


### New functions:

* stepwise

* summaryWald


### Modified functions:

* fuzzyRangeChange
    - added 'x.lab=TRUE' argument, so user can set to FALSE and add labels differently

* stepByStep
    - added warning("'Favourability' is only applicable when family=binomial(link='logit'), so it was automatically set to FALSE.")

* multGLM
    - added argument trim.fun="modelTrim", which can be changed to "stepwise"



# Version 4.5 
### (Committed 2022-06-15)


### Modified functions:

* stepByStep
    - added 'k' argument to pass to 'step' (allows approaching p<0.05; https://www.researchgate.net/post/Why-stepAIC-gives-a-model-with-insignificant-variables-in-the-summarymodel)
    - added 'direction' argument to pass to 'step' (implementation by Alba Estrada)



# Version 4.4 
### (Committed 2022-06-14)


### Modified functions:

* fuzzyConsensus
    - 'cat' message now mentions the number of input vectors

* multGLM
    - added 'cor.method' argument (for 'cor' via 'corSelect'; default "pearson")

* stepByStep
    - 'sp.col' and 'var.cols' can now be column names, not just index numbers



# Version 4.3 
### (Committed 2022-06-01) -> CRAN


### Modified functions:

* fuzzyConsensus
    - 'data' can also be SpatRaster
    - removed unnecessary argument 'fav.cols' (as no other '*.cols' argument exists)



# Version 4.2 
### (Committed 2022-05-31)


### New functions:

* fuzzyConsensus


### Modified functions:

* Fav
    - stops with error if any 'pred' is outside the [0,1] interval (matching 'fuzzyConsensus')


### Other modified files:

* DESCRIPTION
    - added reference about the method, as suggested by Uwe Ligges after previous CRAN submission



# Version 4.1 
### (Committed 2022-05-15)


### Modified functions:

* modOverlap
    - improved help file and removed reference to requiring GLM predictions



# Version 4.0 
### (Committed 2022-05-02) -> CRAN


### New functions:

* prevalence (copied from 'modEvA' pkg)


### Modified functions:

* bioThreat
    - "white" changed to "grey" for 0-values to be visible on the map

* corSelect, FDR, fuzzyOverlay, integerCols, multConvert, multGLM, multicol, multTSA, simMat, splist2presabs, stepByStep, transpose
    - added data=as.data.frame(data), to avoid errors when input is tibble

* bioThreat, Fav, favClass, fuzSim, fuzzyRangeChange, modOverlap, sharedFav, spCodes
    - added 'unlist' to avoid obscure error when input is one-column tibble instead of vector



# Version 3.96 
### (Committed 2022-03-22)


### Modified functions:

* sharedFav
    - improved help file and corrected a mistake (bug report by Alba Estrada)
    - replaced 'main' argument with '...' for barplot (to accept also 'las')
    - plot now produced in phases as computations are made, to have something plotted even in case of error
    - 'conf' can be set to NA, to avoid error and still get values and plot (without CIs)
    - returned value now includes 'FOvI' and 'bins_table'



# Version 3.95 
### (Committed 2022-03-22)


### New functions:

* rarity

* vulnerability

* entropy


### Modified functions:

* spCodes
    - added 'verbosity' argument
    - 'cat' now used instead of red 'message' when OK no dupl codes

* corSelect
    - replaced '...' with arguments 'use' and 'method', the former defaulting to "pairwise.complete.obs"

* simMat
    - added 'na.omit' to 'data' in 'stopifnot', to avoid error when NAs

* sharedFav
    - added 'bin_interval' argument, currently accepting '0.1' (the default, for back-compatibility) or 'quantiles'
    - added Note to help file about possible error when overly small bins



# Version 3.9 
### (Committed 2022-03-13)


### Modified functions:

* gridRecords
    - result now ordered by "cells" column, instead of presences then absences
    - help file "Examples" changed to "elev" raster for more clarity
    - included also example with multiple species

* corSelect
    - added 'family' argument (to pass to FDR when 'sp.col' is not NULL)
    - fixed NA 'excluded.vars' output when input 'var.cols' is character

* FDR
    - added 'Gamma' for positive non-integer responses when family="auto"



# Version 3.8 
### (Committed 2022-02-05)


### Modified functions:

* gridRecords
    - added 'species' argument for multi-species gridding
    - 'pres.coords' (and 'abs.coords' if not NULL) automatically converted with 'as.data.frame'
    - changed 'raster' to 'terra' in help file "Examples"
    - fixed bug that did not grid abs.coords (if not NULL)



# Version 3.7 
### (Committed 2022-01-21) -> CRAN


### Modified functions:

* corSelect
    - fixed new CRAN check error (with 'all.equal' in 'corSelect') on r-devel-linux-x86_64-fedora-clang and r-devel-linux-x86_64-fedora-gcc (by manually setting the 'tolerance' argument for 'all.equal')

* Fav
    - added import modEvA::mod2obspred


### Other modified files:

* fuzsim.Rd
    - updated Sorensen wikipedia link with 'https'



# Version 3.6 
### (Committed 2021-09-29)


### Modified functions:

* Fav
    - 'model' argument can now be of class 'glm', 'gam', 'gbm', 'randomForest' or 'bart'
    - 'pred' can now also be a 'SpatRaster' ('terra' pkg) object

* gridRecords
    - 'rst' can now also be a 'SpatRaster' ('terra' pkg) object
    - Raster* 'rst' coerced to SpatRaster if 'terra' pkg is installed



# Version 3.5 
### (Committed 2021-09-04)


### Modified functions:

* FDR
    - fixed bug when using 'pvalues' as input (following bug report by Stephen via modTools contact form)



# Version 3.4 
### (Committed 2021-09-02)


### Modified functions:

* gridRecords
    - fixed bug when only one raster layer [ , drop = FALSE]



# Version 3.3 
### (Committed 2021-04-24)


### Modified functions:

* simMat
    - added 'verbosity' argument
    - replaced "25%-50%-75% done" with 'utils::txtProgressBar'

* multTSA
    - added 'verbosity' argument

* getPreds
    - added 'verbosity' argument
    - 'data' can now also be a RasterBrick (not only dataframe or a RasterStack)

* Fav
    - added 'verbosity' argument
    - 'pred' can now also be a RasterLayer (not only a numeric vector)



# Version 3.2 
### (Committed 2020-12-12)


### Modified functions:

* favClass (benefitting also bioThreat)
    - added 'is.na(fav) | ' to 'stopifnot' to avoid NA fail


### Other changes:

    - added Linero et al. (2020) to references of papers using fuzzySim



# Version 3.1 
### (Committed 2020-09-18)


### Modified functions:

* gridRecords
    - added 'absences' logical argument



# Version 3.0 
### (Committed 2020-02-03) -> CRAN


### Modified functions:

* gridRecords
    - added 'na.rm' argument (default TRUE)
    - fixed bug when no gridded absences (i.e. no absence AND not presence)


### Other changes:

    -  clarifications and typo/format corrections in the manual



# Version 2.5 
### (Committed 2020-01-31)


### Modified files:

* gridRecords.Rd
    - added examples



# Version 2.2.4 
### (Committed 2020-01-30)


### New functions:

* gridRecords


### Modified functions:

* getPreds:
    - replaced 'if("raster" %in% class(data))' (no longer working) with 'if (is(data, "RasterStack"))'
    - added 'Imports: methods' to DESCRIPTION and 'importFrom("methods", "is")' to NAMESPACE
    - if 'data' are raster, added 'raster::' before 'stack'



# Version 2.2.3 
### (Committed 2020-01-06)


### Modified functions:

* multGLM:
    - fixed id.col name in 'predictions' when id.col input as character

### Other changes:

  - added 'inst' folder with article citation information
  - removed 'onAttach' function with citation message on load
  - added package URLs to DESCRIPTION file



# Version 2.2.2 
### (Committed 2020-01-03)


### Modified functions:

* multGLM:
    - output now includes list of selected variables per model
    - sp.cols, var.cols and id.col can now optionally be input as column names rather than index numbers


### Other changes:

- fixed length of some lines along the PDF manual



# Version 2.2.1 
### (Committed 2019-10-18)


### Modified functions:

* modOverlap:
    - corrected typo in na.rm closing parentheses (reported by Heidi K. Mod)

* Fav:
    - slightly reduced probabilities of exactly 1, which would cause division by zero (resulting Fav is still 1)



# Version 2.2 
### (Committed 2019-03-10)


### New functions:

* sharedFav



# Version 2.1 
### (Committed 2019-03-07)


### New functions:

* favClass
* bioThreat


### Modified functions:

* multGLM:
    - spatial_trend variable in models (when TSA=TRUE and the spatial trend is selected) now named after the response variable (e.g. 'sptrend_giraffe')

* getPreds:
    - fixed new bug by replacing 'if (class(data) == "RasterStack")' with 'if ("RasterStack" %in% class(data))'


### Modified .Rd files:

* multGLM:
    - added example with TSA=TRUE


### Other changes:

- updated maintainer e-mail address



# Version 2.0 
### (Committed 2018-12-05) -> CRAN


### Modified functions:

* multTSA:
    - resulting model object now named 'models' rather than 'TSA.models'
    - corrected bug in coordinate polynomial names (if save.models = TRUE)


# Version 1.9 
### (Committed 2018-11-20)


### Modified functions:

* distPres
    - 'inv' now subtracts from 1 (after standardizing) rather than dividing by 1

* multTSA:
    - added 'criterion' argument, which can be 'AIC' (to use 'step') or 'significance' (to use 'modelTrim')
    - added '...' argument to pass to modelTrim
    - introduced more informative coordinate polynomial names (visible if 'save.models = TRUE')
    - added 'simple = TRUE' to 'poly' for speedup

* modelTrim:
    - added 'data' argument to 'update', to avoid 'attach'

* multGLM:
    - added arguments 'verbosity', 'TSA' and 'coord.cols'
    - added 'data = train.data' to 'step', and 'with(train.data' to 'model.formula', so finally got rid of 'attach'
    - "Building model 1..." 'message' instances replaced with 'cat' so that they are saved if 'sink' is used



# Version 1.8.3 
### (Committed 2018-07-06)


### Modified functions:

* multTSA:
    - added argument 'type' which can be "Y", "P" (the default, for back-compatibility) or "F" (which substitutes the deprecated argument below)
    - deprecated argument 'Favourability = FALSE'


### Modified .Rd files:

* multTSA:
    - documented argument changes described above

* pairwiseRangemaps:
    - updated literature reference


## webpage:

- added additional article citing fuzzySim



# Version 1.8.2 
### (Committed 2018-05-23)


### Modified functions:

* getPreds:
    - 'data' can now be a RasterStack



# Version 1.8.1 
### (Committed 2018-05-15)


### Modified functions:

* multTSA:
    - fixed bug when only one species was used


## webpage:
- added articles citing fuzzySim



# Version 1.8.0 
### (Committed 2017-07-07)


### Modified functions:

* corSelect:
    - included 'VIF' criterion



# Version 1.7.9 
### (Committed 2017-03-27)


### Modified functions:

* multGLM:
    - fixed bug when only one variabe passed to corSelect

* corSelect:
    - added check that 'data' is not missing and is a data frame
    - help file now states that corSelect is included as option in multGLM
    - added 'use = "pairwise.complete.obs"' to examples

* PDF manual:
    - reduced several code line lengths to avoid overboard



# Version 1.7.8 
### (Committed 2016-09-15)


## Removed empty sections from .Rd files

### Modified functions:

* fuzzyRangeChange:
    - underscore replaced with space in measure names
    - x axis labels plotted instead of legend
    - 'col' no longer supplied by default



# Version 1.7.7 
### (Committed 2016-08-01)


### Modified functions:

* pairwiseRangemaps:
    - argument 'chunks' replaced with 'nchunks' for operative reasons
    - added 'subchunks' argument for continuing interrupted runs



# Version 1.7.6 
### (Committed 2016-05-05)


### Modified functions:

* corSelect:
    - removed error when no corrs above threshold
    - sp.cols now NULL by default (to avoid error when missing)

* multGLM:
    - FDR correction reverted to "fdr" by default



# Version 1.7.5 
### (Committed 2016-04-26)


### Modified functions:

* FDR:
    - BIC now also provided

* corSelect:
    - BIC now included as selection criterion

* multGLM:
    - added "select" argument (for 'step') - AIC or BIC



# Version 1.7.4 
### (Committed 2016-04-15)


### Modified functions:

* multGLM:
    - FDR correction now "BY" by default



# Version 1.7.3 
### (Committed 2016-04-12)


### Modified functions:

* corSelect:
    - error message when length(sp.cols) > 1
    - message about excluded missing data
    - added option sp.cols=NULL to get only high.cor.mat

* multGLM:
    - added 'correction' argument to pass to 'FDR'



# Version 1.7.2 
### (Committed 2016-03-22)


### Modified functions:

* fuzzyRangeChange:
    - user can now choose which measures to calculate
    - 'loss' now reported as originally negative
    - results provided also as a barplot



# Version 1.7.1 
### (Committed 2016-03-17)


### Modified functions:

* multTSA:
    - eliminated call to 'attach'



# Version 1.7 
### (Committed 2016-02-12)


### New functions:

* pairwiseRangemaps (calculate area of pairwise intersection and union between rangemaps)

* rangemapSim (calculate  rangemap similarity using common similarity indices)



# Version 1.6.3 
### (Committed 2015-12-02)


### Modified functions:

* FDR:
    - now uses only finite sp.col values



# Version 1.6.2 
### (Committed 2015-11-23)


### Modified functions:

* fuzzyRangeChange:
    - corrected proportional changes to be relative to reference range size, not total study area
    - replaced "stable" with "stable presence" and added (fuzzy equivalent of) "stable absence"
    - replaced "change" with "balance" (overall loss/gain, not amount of changed cells)
    - result is now data frame, not named vector

* distPres:
    - now allows NA values



# Version 1.6.1 
### (Committed 2015-11-13)


### Modified functions:

* fuzzyOverlay:
    - added arguments 'prop = TRUE' and 'overlay.cols = 1:ncol(data)'
    - result now includes 4 (not 8) values, either sum or proportion
    - result is now a named vector instead of a list
    - 'na.rm' now TRUE by default

* fuzSim:
    - added argument 'na.rm = TRUE'

* multGLM:
    - bug corrected in corSelect (was reporting but not really eliminating variables)
    - suppressMessages in corSelect


### Modified help files:

* fuzzyOverlay: modifications reflecting function changes

* modOverlap: example now provided

* fuzzyOverlay, fuzzyRangeChange, modOverlap:
    - examples containing 'jitter' corrected to avoid values outside the [0, 1] interval

* fuzSim:
    - examples now provided also for similarity between fuzzy data
    - tables of significance for Jaccard and Baroni's indices now referred



# Version 1.6 
### (Committed 2015-11-03)


### New functions:

* fuzzyOverlay (calculate row-wise intersection, union, expansion, contraction or consensus among continuous model predictions)

* fuzzyRangeChange (calculate overal loss, gain, and maintenance of favourability between models)



# Version 1.5 
### (Committed 2015-10-29)


### New functions:

* corSelect (select among correlated variables based on their bivariate relationship with the response)

* modOverlap (asses the total overlap between model predictions using niche comparison metrics)


### Modified functions:

* multGLM:
    - 'corSelect' now included as additional option for variable selection

* FDR:
    - AIC now also calculated
    - 'model.type' deprecated
    - 'family = "auto"' by default
    - 'simplif' argument added



# Version 1.4
## ... and previous (edits I can remember)


### Modified functions:

* multicol:
    - 'model' argument added (user can provide a model object instead of a set of variables)
    - variables in output are now ordered according to VIF

* simFromSetOps:
    - the similarity index used is now mentioned in a message

* FDR:
    - data input format changed, with former parameters 'response' and 'predictors' replaced with 'data', 'sp.cols' and 'var.cols' (for coherence and compatibility with 'multGLM' function)
