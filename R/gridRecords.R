gridRecords <- function(rst,
                        pres.coords,
                        abs.coords = NULL,
                        keep.n = FALSE
) {
  
  if (is.null(abs.coords)) abs.coords <- coordinates(rst)
  
  p <- pres.coords
  a <- abs.coords
  
  p_extract <- raster::extract(rst, p, cellnumbers = TRUE, df = TRUE)[ , -1]
  a_extract <- raster::extract(rst, a, cellnumbers = TRUE, df = TRUE)[ , -1]
  
  if (!keep.n) p_extract <- unique(p_extract)
  a_extract <- unique(a_extract)
  
  a_extract <- a_extract[!(a_extract$cells %in% p_extract$cells), ]
  
  p_centroids <- xyFromCell(rst, p_extract$cells)
  a_centroids <- xyFromCell(rst, a_extract$cells)
  
  p_extract <- data.frame(presence = 1, p_centroids, p_extract)
  a_extract <- data.frame(presence = 0, a_centroids, a_extract)
  
  result <- rbind(p_extract, a_extract)
}
