require(raster)
require(rgdal)
require(MASS)
require(pROC)
require(ggplot2)

setwd('~/Documents/CSRM/JAS/Data/')

wet <- raster('./CD_AOI_All_Wetlands/Irwin_CD_AOI_All_Wetlands.tif')
dem <- raster('./DEM_AOI/Irwin_DEM_AOI.tif')
lrm <- raster('./LRM_AOI/Irwin_LRM_AOI.tif')
sbs <- raster('./Site_Boundaries_Lithic_Scatter/Irwin_Prehistoric_NREI_Lithic_Scatter_AOI.tif')
slp <- raster('./Slope_AOI/Irwin_Slope_AOI.tif')
wv3 <- raster::stack('./Resampled_WV3/15SEP21184609-M3DS_R01C1-057029324010_01_P001_radiance_resample_AOI.tif')
# smi <- raster::MoranLocal(slp, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3))
# wmi <- raster::MoranLocal(wet, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3))
# dmi <- raster::MoranLocal(dem, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3))
# lmi <- raster::MoranLocal(lrm, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3))
smi <- raster('./Slope_Moran_I/slope_moran_i.tif')
wmi <- raster('./Wet_Moran_I/wet_moran_i.tif')
dmi <- raster('./Distance_Moran_I/distance_moran_i.tif')
lmi <- raster('./LRM_Moran_I/LRM_moran_i.tif')

# Calculate Moran's I of first principal component on the WorldView-2 bands
# wv3.df <- as.data.frame(wv3)
# colnames(wv3.df) <- c('red', 'redge', 'coastal', 'blue', 'green', 'yellow', 'nir1', 'nir2')
# wv3.pcs <- prcomp(~ red + redge + coastal + blue + green + yellow + nir1 +  nir2, data = wv3.df,
#                                center = TRUE, scale = TRUE)
# wv3.coords <- predict(wv3.pcs, newdata = wv3.df)
# wv3.df$pc.1 <- wv3.coords[, 1]
# wv.pc.1 <- raster(nrow = nrow(wv3), ncol = ncol(wv3), vals = wv3.df$pc.1, crs = crs(wv3))
# wvi <- raster::MoranLocal(wv.pc.1, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3))
# extent(wvi) <- extent(wv3)

# with.indices <- sbs
# with.indices[Which(with.indices == 128, cells = TRUE)] <- 0  # changes non-site value from 128 to 0

# data.set <- raster::stack(with.indices, wet, dem, lrm, slp, wv3, smi, wmi, dmi, lmi, wvi)
# df <- as.data.frame(data.set, na.rm = TRUE)
# colnames(df) <- c('ID', 'wet', 'dem', 'lrm', 'slp',
#                         'red', 'redge', 'coastal', 'blue', 'green',
#                       'yellow', 'nir1', 'nir2', 'slope.moran', 'wet.moran',
#                         'distance.moran', 'lrm.moran', 'wv3.moran')
# pcs <- prcomp(~ red + redge + coastal + blue + green + yellow + nir1 +  nir2, data = df,
              center = TRUE, scale = TRUE)
# df.coords <- predict(pcs, newdata = df[, 6:13])
# df$pc.1 <- df.coords[, 1]
# df$pc.2 <- df.coords[, 2]
# df$pc.3 <- df.coords[, 3]

# save(df, file = './site_df_indices.RData')

load('./site_df_indices.RData')  # simply loads in the data frame created by the ~25 lines just above

df$ID <- as.factor(df$ID)
print(summary(df$ID))
site.ids <- as.integer(names(summary(df$ID)[-1]))

moran.m <- matrix(nrow = length(site.ids), ncol = 6)
moran.m[, 1] <- site.ids

morans <- c('slope.moran', 'wet.moran', 'distance.moran', 'lrm.moran', 'wv3.moran')

i <- 1
j <- 2

for(k in morans) {
  for(l in site.ids) {
    moran.m[i, j] <- mean(df[which(df$ID == l), k])
    i <- i + 1
  }
  i <- 1
  j <- j + 1
}



get.site.vals <- function(site.1, site.2, var.name, total.size) {
  site.1.vals <- numeric(length = total.size)
  site.2.vals <- numeric(length = total.size)
  
  site.1.in <- df[which(df$ID == site.1), var.name]
  site.2.in <- df[which(df$ID == site.2), var.name]
  
  for(i in 1:total.size) {
    s.1.samp <- sample(site.1.in, 20)
    s.2.samp <- sample(site.2.in, 20)
    
    s.1.mean <- mean(s.1.samp)
    s.2.mean <- mean(s.2.samp)
    
    site.1.vals[i] <- s.1.mean
    site.2.vals[i] <- s.2.mean
  }
  out.df <- data.frame(site.1.vals, site.2.vals)
  t <- t.test(out.df$site.1.vals, out.df$site.2.vals, alternative = 'two.sided')
  return(t$p.value)
}

pairwise.matrix <- matrix(nrow = 12, ncol = 12, dimnames = list(Site.1 = c(site.ids),
                                                                Site.2 = c(site.ids)))

for(i in 1:nrow(pairwise.matrix)) {
  for(j in i:ncol(pairwise.matrix)) {
    pairwise.matrix[i, j] <- get.site.vals(site.ids[i], site.ids[j], 'slp', 100)
  }
}

slope.matrix <- pairwise.matrix

pairwise.matrix <- matrix(nrow = 12, ncol = 12, dimnames = list(Site.1 = c(site.ids),
                                                                Site.2 = c(site.ids)))

for(i in 1:nrow(pairwise.matrix)) {
  for(j in i:ncol(pairwise.matrix)) {
    pairwise.matrix[i, j] <- get.site.vals(site.ids[i], site.ids[j], 'pc.1', 100)
  }
}

pc.1.matrix <- pairwise.matrix


wet.matrix <- matrix(nrow = 12, ncol = 12, dimnames = list(Site.1 = c(site.ids),
                                                           Site.2 = c(site.ids)))

for(i in 1:nrow(wet.matrix)) {
  for(j in i:ncol(wet.matrix)) {
    wet.matrix[i, j] <- get.site.vals(site.ids[i], site.ids[j], 'wet', 100)
  }
}


dem.matrix <- matrix(nrow = 12, ncol = 12, dimnames = list(Site.1 = c(site.ids),
                                                           Site.2 = c(site.ids)))

for(i in 1:nrow(dem.matrix)) {
  for(j in i:ncol(dem.matrix)) {
    dem.matrix[i, j] <- get.site.vals(site.ids[i], site.ids[j], 'dem', 100)
  }
}


lrm.matrix <- matrix(nrow = 12, ncol = 12, dimnames = list(Site.1 = c(site.ids),
                                                           Site.2 = c(site.ids)))

for(i in 1:nrow(lrm.matrix)) {
  for(j in i:ncol(lrm.matrix)) {
    lrm.matrix[i, j] <- get.site.vals(site.ids[i], site.ids[j], 'lrm', 100)
  }
}



