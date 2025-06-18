getRegion <- function (pres.coords, type = "width", clust_dist = 100, dist_mult = 1,
          width_mult = 0.5, weight = FALSE, CRS = NULL, dist_mat = NULL,
          verbosity = 2, plot = TRUE)
{
  if (!(type %in% c("width", "mean_dist", "inv_dist", "clust_mean_dist",
                    "clust_width")))
    stop("Invalid 'type'. See help file for options.")
  stopifnot(dist_mult > 0, width_mult > 0, is.numeric(verbosity),
            is.logical(plot))
  if (!inherits(pres.coords, "SpatVector")) {
    pres.coords <- as.data.frame(pres.coords)
    if (ncol(pres.coords) > 2)
      stop("If not a SpatVector, 'pres.coords' must have only two columns")
    pres.coords <- terra::vect(pres.coords, geom = colnames(pres.coords))
  }
  if (inherits(pres.coords, "SpatVector")) {
    if (isFALSE(terra::is.points(pres.coords)))
      stop("If 'pres.coords' is of class 'SpatVector', its 'geomtype' must be 'points'.")
    if (!is.null(CRS)) {
      if (terra::crs(pres.coords) == "") {
        terra::set.crs(pres.coords, CRS)
      }
      else {
        message("'CRS' argument ignored, as 'pres.coords' already has one.")
      }
    }
  }
  if (grepl("dist", type)) {
    if (is.null(dist_mat)) {
      if (verbosity > 0)
        message("Computing pairwise distance between points...")
      dist_mat <- terra::distance(pres.coords)
    }
    else {
      if (verbosity > 0)
        message("Using supplied pairwise distance between points...")
    }
    dist_mean <- mean(dist_mat, na.rm = TRUE)
  }
  if (grepl("clust", type)) {
    if (verbosity > 1)
      message("Computing point clusters...")
    # tree <- stats::hclust(d = dist_mat, method = "single")
    # pres.coords$clust <- stats::cutree(tree, h = clust_dist *
    #                                      1000)
    buffers <- terra::buffer(pres.coords, width = clust_dist * 1000)
    groups <- terra::disagg(terra::aggregate(buffers))
    groups$id <- 1:nrow(groups)
    pres.coords$clust <- terra::extract(groups, pres.coords)$id
    clusters <- unique(pres.coords$clust)
  }
  if (type == "mean_dist") {
    reg <- terra::buffer(pres.coords, width = dist_mean *
                           dist_mult)
  }
  else if (type == "inv_dist") {
    dist_df <- as.data.frame(as.matrix(dist_mat))
    dist_sums <- sapply(dist_mat, sum, na.rm = TRUE)
    range01 <- function(x) {
      (x - min(x))/(max(x) - min(x))
    }
    dist_sums_01 <- range01(dist_sums)
    dist_sums_01[dist_sums_01 == 0] <- 0.001
    reg <- terra::buffer(pres.coords, width = dist_mean *
                           rev(dist_sums_01) * dist_mult)
  }
  else if (type == "clust_mean_dist") {
    if (verbosity > 0)
      message("Computing pairwise distance within clusters...")
    for (i in clusters) {
      clust_pts <- pres.coords[pres.coords$clust == i,
      ]
      buff_radius <- mean(terra::distance(clust_pts)) *
        dist_mult
      if (!is.finite(buff_radius) || buff_radius <= 0)
        buff_radius <- 0.001
      pres.coords[pres.coords$clust == i, "buff_radius"] <- buff_radius
    }
    if (isTRUE(weight)) {
      counts <- table(pres.coords$clust)
      clust_n <- counts[match(pres.coords$clust, names(counts))]
      clust_abund <- (clust_n - min(clust_n))/(max(clust_n) -
                                                 min(clust_n))
      clust_abund[clust_abund == 0] <- 0.001
      pres.coords$buff_radius <- pres.coords$buff_radius *
        clust_abund
    }
    reg <- terra::buffer(pres.coords, width = "buff_radius")
  }
  else if (type == "clust_width") {
    if (verbosity > 0)
      message("Computing cluster widths...")
    for (i in clusters) {
      clust_pts <- pres.coords[pres.coords$clust == i,
      ]
      clust_width <- terra::width(terra::aggregate(clust_pts)) *
        width_mult
      if (clust_width <= 0)
        clust_width <- 0.001
      pres.coords[pres.coords$clust == i, "clust_width"] <- clust_width
    }
    if (isTRUE(weight)) {
      counts <- table(pres.coords$clust)
      clust_n <- counts[match(pres.coords$clust, names(counts))]
      clust_abund <- (clust_n - min(clust_n))/(max(clust_n) -
                                                 min(clust_n))
      clust_abund[clust_abund == 0] <- 0.001
      pres.coords$clust_width <- pres.coords$clust_width *
        clust_abund
    }
    reg <- terra::buffer(pres.coords, width = "clust_width")
  }
  else if (type == "width") {
    wdth <- terra::width(terra::aggregate(pres.coords))
    reg <- terra::buffer(pres.coords, width = wdth * width_mult)
  }
  reg <- terra::aggregate(reg)
  if (plot) {
    terra::plot(reg, col = "yellow", border = NA, main = paste("type =",
                                                               type))
    if (!grepl("clust", type)) {
      terra::plot(pres.coords, cex = 0.3, add = TRUE)
    }
    else {
      clust_aggregates <- terra::aggregate(pres.coords,
                                           "clust")
      clust_centroids <- terra::centroids(clust_aggregates)
      terra::plot(clust_aggregates, cex = 0.3, add = TRUE,
                  col = grDevices::hcl.colors(length(clusters),
                                              palette = "dark2"))
      terra::text(clust_centroids, "agg_n", cex = 0.5)
    }
  }
  if (verbosity > 1 && grepl("dist|clust", type))
    message("Finished!")
  return(reg)
}
