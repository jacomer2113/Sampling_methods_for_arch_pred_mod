require(raster)
require(rgdal)
require(MASS)
require(pROC)

setwd('~/Documents/CSRM/JAS/Data/')

# load in raster data layers
wet <- raster('./CD_AOI_All_Wetlands/Irwin_CD_AOI_All_Wetlands.tif')
dem <- raster('./DEM_AOI/Irwin_DEM_AOI.tif')
lrm <- raster('./LRM_AOI/Irwin_LRM_AOI.tif')
sbs <- raster('./Site_Boundaries_Lithic_Scatter/Irwin_Prehistoric_NREI_Lithic_Scatter_AOI.tif')
slp <- raster('./Slope_AOI/Irwin_Slope_AOI.tif')
wv3 <- raster::stack('./Resampled_WV3/15SEP21184609-M3DS_R01C1-057029324010_01_P001_radiance_resample_AOI.tif')

# In the original site layer, the individual sites have ID numbers. The cells composing each site have
# values equal to their site's ID. The following lines turn the site layer into a binary, where all
# site cells, regardless of which site they belong to, equal 1, and all else is NA
# Non-site value in the original layer is 128
sbs.binary <- Which(sbs != 128)

# This layer has cell values = 1 at sites and within a 200 m buffer around them, and NA elsewhere. We can sample
# from outside this buffer (from the cells that are now NA) to get non-site cells that are a safe distance from
# site cells
# I'm loading in the ready-made raster here, because the buffer function takes a while. The code that
# actually creates the buffer is the next line, commented out
# big.buff.sites <- buffer(sbs.binary, width = 200, doEdge = TRUE)
big.buff.sites <- raster('./200m_buffered_lith_scatt.tif')

# In the site buffer raster layer, switch the cell values -- that is, change 1s to NAs and NAs to 1s
big.buff.sites[Which(big.buff.sites == 1, cells = TRUE)] <- 2
big.buff.sites[Which(is.na(big.buff.sites), cells = TRUE)] <- 1
big.buff.sites[Which(big.buff.sites == 2, cells = TRUE)] <- NA

# Define a raster extent that is 50 meters inward from the borders of the full-size rasters, so that
# none of the non-site areas we select are too close to the edge
ext.100m.edge <- extent(big.buff.sites, 50, (nrow(big.buff.sites) - 50),
                                                    50, (ncol(big.buff.sites) - 50))

# Summary of sites:
site.summary <- summary(as.factor(values(sbs)))
site.summary <- site.summary[1:(length(site.summary) - 1)]
site.summary <- matrix(ncol = 2, data = as.numeric(c((names(site.summary)), (unname(site.summary)))),
                       dimnames = list(1:length(site.summary), c('Site number', 'Number of cells')))
print(site.summary)

# Selecting non-site cells at random from within the smaller extent defined above
ns.cells.matrix <- matrix(nrow = nrow(site.summary), ncol = 3)
s <- 1
set.seed(2113)
while(s <= nrow(site.summary)) {
  cur.cell.samp <- sampleRandom(big.buff.sites, size = 1, ext = ext.100m.edge, na.rm = FALSE, cells = TRUE)
  if(!is.na(cur.cell.samp[2])) {
    ns.cells.matrix[s, 1] <- cur.cell.samp[1]
    ns.cells.matrix[s, 2] <- cur.cell.samp[2]
    ns.cells.matrix[s, 3] <- 1000 + s
    s <- s + 1
    print('success')
  } else {
    print('site cell')
    next
  }
}

# Get the mean size of sites
mean.site.size <- round(mean(site.summary[, 2]))
ns.radius <- round(sqrt(((mean.site.size * 4) / 3.14)))

# Buffer the randomly selected non-site cells so that the non-sites are approximately equal to
# the mean site size
ns.layer <- big.buff.sites
ns.layer[] <- NA
ns.layer[ns.cells.matrix[, 1]] <- 1
ns.layer <- buffer(ns.layer, width = ns.radius)

# Add non-sites to the binary site layer created up top, so that all non-site cells have value 0
sbs.binary[Which(sbs.binary == 0, cells = TRUE)] <- NA
sbs.binary[Which(!is.na(ns.layer), cells = TRUE)] <- 0

# Site and non-site raster:
raster::plot(sbs.binary)

# Assign individual non-sites individual IDs, then make all cell values within each non-site equal that
# ID (and all other cells equal 0)
ns.stack <- sbs
ns.stack[] <- NA
ns.stack[ns.cells.matrix[1, 1]] <- 1
ns.stack <- buffer(ns.stack, width = ns.radius)
ns.stack[Which(!is.na(ns.stack), cells = TRUE)] <- ns.cells.matrix[1, 3]
ns.stack[Which(is.na(ns.stack), cells = TRUE)] <- 0

# for(i in 2:nrow(ns.cells.matrix)) {
#   cur.rast <- sbs
#   cur.rast[] <- NA
#   cur.rast[ns.cells.matrix[i, 1]] <- 1
#   cur.rast <- buffer(cur.rast, width = ns.radius)
#   cur.rast[Which(!is.na(cur.rast), cells = TRUE)] <- ns.cells.matrix[i, 3]
#   cur.rast[Which(is.na(cur.rast), cells = TRUE)] <- 0
#   
#   ns.stack <- raster::stack(ns.stack, cur.rast)
# }

# ns.layer.w.ids <- sum(ns.stack)
# writeRaster(ns.layer.w.ids, filename = './ns_layer_w_ids.tif')

ns.layer.w.ids <- raster('./ns_layer_w_ids.tif')
ns.layer.w.ids[Which(is.na(ns.layer.w.ids), cells = TRUE)] <- 0

# Make site layer where all site cells have values that equal their site's ID, and all non-site cells are
# 0
site.layer.w.ids <- sbs
site.layer.w.ids[Which(site.layer.w.ids == 128, cells = TRUE)] <- 0

# Make single layer for all sites and non-sites, including individual site and non-site IDs
combined.site.ns <- ns.layer.w.ids + site.layer.w.ids

site.layer.w.ids[Which(site.layer.w.ids == 0, cells = TRUE)] <- NA
site.stack <- raster::stack(site.layer.w.ids, wet, dem, lrm, slp, wv3)
site.layer.df <- as.data.frame(site.stack, na.rm = TRUE)

ns.layer.w.ids[Which(ns.layer.w.ids == 0, cells = TRUE)] <- NA
ns.stack <- raster::stack(ns.layer.w.ids, wet, dem, lrm, slp, wv3)
ns.layer.df <- as.data.frame(ns.stack, na.rm = TRUE)

colnames(site.layer.df) <- colnames(ns.layer.df) <- c('ID', 'wet', 'dem', 'lrm', 'slp',
                                                      'red', 'redge', 'coastal', 'blue', 'green',
                                                      'yellow', 'nir1', 'nir2')

site.layer.df$site <- factor(1)
ns.layer.df$site <- factor(0)
full.df <- rbind(site.layer.df, ns.layer.df)

# Calculate principal components of spectral measurements
pcs <- prcomp(~ red + redge + coastal + blue + green + yellow + nir1 +  nir2, data = full.df,
              center = TRUE, scale = TRUE)
summary(pcs)
rand.samp.coords <- predict(pcs, newdata = full.df[, 6:13])
full.df$pc.1 <- rand.samp.coords[, 1]

# Summarize non-sites
site.ns.summary <- summary(as.factor(values(combined.site.ns)))
site.ns.summary <- site.ns.summary[1:(length(site.ns.summary))]
site.ns.summary <- matrix(ncol = 2, data = as.numeric(c((names(site.ns.summary)), (unname(site.ns.summary)))),
                     dimnames = list(1:length(site.ns.summary), c('Site number', 'Number of cells')))

site.ns.summary <- site.ns.summary[2:nrow(site.ns.summary), ]
print(site.ns.summary)

# Make a matrix for LOOCV lda results
loocv.results.lda <- matrix(c(0, 0, 0, 0), nrow = 2, dimnames = list(Sitehood = c('Site', 'Non-site'),
                                                                 Prediction = c('Site', 'Non-site')))

set.seed(2113)
class.fun <- formula(site ~ wet + dem + lrm + slp + pc.1)

loocv.res.table.lda <- matrix(nrow = nrow(site.ns.summary), ncol = 2,
                           dimnames = list(NULL, c('True', 'Predicted')))

for(i in 1:nrow(site.ns.summary)) {
  left.out.site <- site.ns.summary[i, 1]
  print(paste0('left out site: ', left.out.site))
  
  training.data <- full.df[which(full.df$ID != left.out.site), ]  # randomly sample half of these for each iteration of LOOCV model?
  testing.data <- full.df[which(full.df$ID == left.out.site), ]
  
  class.mod <- lda(class.fun, data = training.data)
  res <- predict(class.mod, testing.data)$class
  print(summary(res))
  
  summary <- unname(summary(res))
  
  if(summary[1] < summary[2]) {
    site.prediction <- FALSE
    loocv.res.table.lda[i, 2] <- 0
  } else {
    site.prediction <- TRUE
    loocv.res.table.lda[i, 2] <- 1
  }
  
  if(left.out.site < 1000) {
    true.sitehood <- TRUE
    loocv.res.table.lda[i, 1] <- 1
  } else {
    true.sitehood <- FALSE
    loocv.res.table.lda[i, 1] <- 0
  }
  
  if(site.prediction == TRUE & true.sitehood == TRUE) {
    loocv.results.lda[1, 1] <- loocv.results.lda[1, 1] + 1
    print('true positive')
  } else if(site.prediction == TRUE & true.sitehood == FALSE) {
    loocv.results.lda[2, 1] <- loocv.results.lda[2, 1] + 1
    print('false positive')
  } else if(site.prediction == FALSE & true.sitehood == TRUE) {
    loocv.results.lda[1, 2] <- loocv.results.lda[1, 2] + 1
    print('false negative')
  } else {
    loocv.results.lda[2, 2] <- loocv.results.lda[2, 2] + 1
    print('true negative')
  }
  
}

loocv.res.df.lda <- as.data.frame(loocv.res.table.lda)
loocv.roc.lda <- roc(loocv.res.df.lda$True, loocv.res.df.lda$Predicted, auc = TRUE, ci = TRUE, plot = TRUE)

# Make a matrix for LOOCV log reg results
loocv.results.log <- matrix(c(0, 0, 0, 0), nrow = 2, dimnames = list(Sitehood = c('Site', 'Non-site'),
                                                                     Prediction = c('Site', 'Non-site')))

set.seed(2113)

loocv.res.table.log <- matrix(nrow = nrow(site.ns.summary), ncol = 2,
                              dimnames = list(NULL, c('True', 'Predicted')))

for(i in 1:nrow(site.ns.summary)) {
  left.out.site <- site.ns.summary[i, 1]
  print(paste0('left out site: ', left.out.site))
  
  training.data <- full.df[which(full.df$ID != left.out.site), ]
  testing.data <- full.df[which(full.df$ID == left.out.site), ]
  
  class.mod <- glm(class.fun, family = binomial, data = training.data)
  res <- predict(class.mod, testing.data, type = 'response')
  pred.classes.log <- res
  n <- length(pred.classes.log[which(pred.classes.log >= 0.5)])
  s <- length(pred.classes.log[which(pred.classes.log <= 0.5)])
  
  print(paste0('nonsite votes: ', n))
  print(paste0('site votes: ', s))

  if(n < s) {
    site.prediction <- TRUE
    loocv.res.table.log[i, 2] <- 1
  } else {
    site.prediction <- FALSE
    loocv.res.table.log[i, 2] <- 0
  }
  
  if(left.out.site < 1000) {
    true.sitehood <- TRUE
    loocv.res.table.log[i, 1] <- 1
  } else {
    true.sitehood <- FALSE
    loocv.res.table.log[i, 1] <- 0
  }
  
  if(site.prediction == TRUE & true.sitehood == TRUE) {
    loocv.results.log[1, 1] <- loocv.results.log[1, 1] + 1
    print('true positive')
  } else if(site.prediction == TRUE & true.sitehood == FALSE) {
    loocv.results.log[2, 1] <- loocv.results.log[2, 1] + 1
    print('false positive')
  } else if(site.prediction == FALSE & true.sitehood == TRUE) {
    loocv.results.log[1, 2] <- loocv.results.log[1, 2] + 1
    print('false negative')
  } else {
    loocv.results.log[2, 2] <- loocv.results.log[2, 2] + 1
    print('true negative')
  }
  
}

loocv.res.df.log <- as.data.frame(loocv.res.table.log)

save.image('./LOOCV.RData')
# Remove the big objects we used for LOOCV, just to save space
rm(training.data, testing.data, full.df, site.layer.df, ns.layer.df)

# Start random cell sampling strategy
site.1.ns.0 <- Which(sbs != 128)
cellnums <- site.1.ns
cellnums[] <- 1:ncell(cellnums)
binary.stack <- raster::stack(site.1.ns.0, wet, dem, lrm, slp, wv3, cellnums)
binary.stack.df <- as.data.frame(binary.stack)
binary.stack.df <- binary.stack.df[complete.cases(binary.stack.df), ]
colnames(binary.stack.df) <- c('site', 'wet', 'dem', 'lrm', 'slp',
                            'red', 'redge', 'coastal', 'blue', 'green',
                            'yellow', 'nir1', 'nir2', 'index')
binary.stack.df$site <- as.factor(binary.stack.df$site)

pcs <- prcomp(~ red + redge + coastal + blue + green + yellow + nir1 +  nir2, data = binary.stack.df,
              center = TRUE, scale = TRUE)
summary(pcs)
binary.coords <- predict(pcs, newdata = binary.stack.df[, 6:13])
binary.stack.df$pc.1 <- binary.coords[, 1]

site.cells <- which(binary.stack.df$site == 1)
ns.cells <- which(binary.stack.df$site == 0)
samp.size <- round(length(site.cells) / 2)

training.site.cell.samp <- sample(1:length(site.cells), samp.size)
training.ns.cell.samp <- sample(1:length(ns.cells), samp.size)

training.site.cells <- site.cells[training.site.cell.samp]
training.ns.cells <- ns.cells[training.ns.cell.samp]

testing.site.cells <- site.cells[-training.site.cell.samp]
testing.ns.cells <- sample(ns.cells[-training.ns.cell.samp], samp.size)

training.df <- rbind(binary.stack.df[training.site.cells, ], binary.stack.df[training.ns.cells, ])
testing.df <- rbind(binary.stack.df[testing.site.cells, ], binary.stack.df[testing.ns.cells, ])

binary.mod.lda <- lda(class.fun, data = training.df)
binary.res.lda <- predict(binary.mod.lda, testing.df)

binary.mod.log <- glm(class.fun, data = training.df, family = binomial)
binary.res.log <- predict(binary.mod.log, testing.df, type = 'response')

pred.classes.lda <- binary.res.lda$class
pred.classes.log <- binary.res.log
pred.classes.log[which(pred.classes.log <= 0.5)] <- 0
pred.classes.log[which(pred.classes.log >= 0.5)] <- 1

rand.samp.res.lda.df <- data.frame(testing.df$site, pred.classes.lda)
colnames(rand.samp.res.lda.df) <- c('True', 'Predicted')
rand.samp.res.lda.df$True <- as.numeric(rand.samp.res.lda.df$True) - 1
rand.samp.res.lda.df$Predicted <- as.numeric(rand.samp.res.lda.df$Predicted) - 1

rand.samp.res.log.df <- data.frame(testing.df$site, pred.classes.log)
colnames(rand.samp.res.log.df) <- c('True', 'Predicted')
rand.samp.res.log.df$True <- as.numeric(rand.samp.res.log.df$True) - 1

# Make a matrix for random cell sampling results lda
rand.samp.results.lda <- matrix(c(0, 0, 0, 0), nrow = 2, dimnames = list(Sitehood = c('Site', 'Non-site'),
                                                                     Prediction = c('Site', 'Non-site')))
for(i in 1:nrow(rand.samp.res.lda.df)) {
  if(rand.samp.res.lda.df[i, 'True'] == 1 & rand.samp.res.lda.df[i, 'Predicted'] == 1) {
    rand.samp.results.lda[1, 1] <- rand.samp.results.lda[1, 1] + 1
  } else if(rand.samp.res.lda.df[i, 'True'] == 0 & rand.samp.res.lda.df[i, 'Predicted'] == 1) {
    rand.samp.results.lda[2, 1] <- rand.samp.results.lda[2, 1] + 1
  } else if(rand.samp.res.lda.df[i, 'True'] == 1 & rand.samp.res.lda.df[i, 'Predicted'] == 0) {
    rand.samp.results.lda[1, 2] <- rand.samp.results.lda[1, 2] + 1
  } else {
    rand.samp.results.lda[2, 2] <- rand.samp.results.lda[2, 2] + 1
  }
}

# Make a matrix for random cell sampling results log reg
rand.samp.results.log <- matrix(c(0, 0, 0, 0), nrow = 2, dimnames = list(Sitehood = c('Site', 'Non-site'),
                                                                         Prediction = c('Site', 'Non-site')))

for(i in 1:nrow(rand.samp.res.log.df)) {
  if(rand.samp.res.log.df[i, 'True'] == 1 & rand.samp.res.log.df[i, 'Predicted'] == 1) {
    rand.samp.results.log[1, 1] <- rand.samp.results.log[1, 1] + 1
  } else if(rand.samp.res.log.df[i, 'True'] == 0 & rand.samp.res.log.df[i, 'Predicted'] == 1) {
    rand.samp.results.log[2, 1] <- rand.samp.results.log[2, 1] + 1
  } else if(rand.samp.res.log.df[i, 'True'] == 1 & rand.samp.res.log.df[i, 'Predicted'] == 0) {
    rand.samp.results.log[1, 2] <- rand.samp.results.log[1, 2] + 1
  } else {
    rand.samp.results.log[2, 2] <- rand.samp.results.log[2, 2] + 1
  }
}

loocv.roc.lda <- roc(loocv.res.df.lda$True, loocv.res.df.lda$Predicted, auc = TRUE, ci = TRUE, plot = TRUE)
loocv.roc.log <- roc(loocv.res.df.log$True, loocv.res.df.log$Predicted, auc = TRUE, ci = TRUE, plot = TRUE)

rand.samp.roc.lda <- roc(rand.samp.res.lda.df$True, rand.samp.res.lda.df$Predicted, auc = TRUE, ci = TRUE, plot = TRUE)
rand.samp.roc.log <- roc(rand.samp.res.log.df$True, rand.samp.res.log.df$Predicted, auc = TRUE, ci = TRUE, plot = TRUE)

calc.gain <- function(m) {
  tp <- m[1, 1]
  fp <- m[2, 1]
  fn <- m[1, 2]
  all <- m[1, 1] + m[1, 2] + m[2, 1] + m[2, 2]
  
  num <- ((tp + fp) / all)
  den <- (tp / (tp + fn))
  
  gain <- 1 - (num / den)
}

loocv.lda.gain <- calc.gain(loocv.results.lda)
loocv.log.gain <- calc.gain(loocv.results.log)
rand.samp.lda.gain <- calc.gain(rand.samp.results.lda)
rand.samp.log.gain <- calc.gain(rand.samp.results.log)


vis.training <- site.1.ns.0
vis.site.cells <- Which(vis.training == 1, cells = TRUE)
vis.ns.cells <- Which(vis.training == 0, cells = TRUE)

vis.site.cell.samp <- sample(vis.site.cells, samp.size)
vis.ns.cell.samp <- sample(vis.ns.cells, samp.size)

vis.training[c(vis.site.cell.samp, vis.ns.cell.samp)] <- 2


# Our Fort Irwin original sampling method: spatially random sample across all cells
