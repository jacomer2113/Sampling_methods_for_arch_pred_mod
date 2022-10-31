load('~/Documents/CSRM/JAS/Data/LOOCV.RData')

require(pROC)
require(raster)
require(rgdal)
require(MASS)

set.seed(2113)

runs <- 200

aucs <- numeric(length = runs)

class.fun <- formula(site ~ wet + dem + lrm + slp + pc.1)

for(j in 1:runs) {
  
  Status <- numeric()
  Pred <- numeric()
  res.df <- data.frame(Status, Pred)
  
  for(i in 1:nrow(site.ns.summary)) {
    left.out.site <- site.ns.summary[i, 1]
    print(paste0('left out site: ', left.out.site))
    
    training.data <- full.df[which(full.df$ID != left.out.site), ]
    training.data <- training.data[sample((1:nrow(training.data)), size = round((nrow(training.data) / 4)), replace = FALSE), ]
    testing.data <- full.df[which(full.df$ID == left.out.site), ]
    testing.data <- testing.data[sample((1:nrow(testing.data)), size = round((nrow(testing.data) / 4)), replace = FALSE), ]
    
    class.mod <- glm(class.fun, family = binomial, data = training.data)
    pred.classes.log <- predict(class.mod, testing.data, type = 'response')
    pred.classes.log[which(pred.classes.log >= 0.5)] <- 1
    pred.classes.log[which(pred.classes.log <= 0.5)] <- 0
    
    cur <- data.frame(as.numeric(testing.data$site), pred.classes.log)
    colnames(cur) <- c('Status', 'Pred')
    
    res.df <- rbind(res.df, cur)
  }
  
  
  a <- roc(res.df$Status, res.df$Pred, auc = TRUE)
  aucs[j] <- a$auc[1]
  
  print(paste0('AUC: ', aucs[j]))
  
  rm(res.df)
  
  print(paste0('Run ', j, ' of ', runs, ' finished'))
  
}

loocv.log.aucs <- aucs

save(loocv.log.aucs, file = 'loocv_log_aucs.RData')

rm(list = ls())

load('~/Documents/CSRM/JAS/Data/LOOCV.RData')

set.seed(2113)

runs <- 200

aucs <- numeric(length = runs)

class.fun <- formula(site ~ wet + dem + lrm + slp + pc.1)

for(j in 1:runs) {
  
  Status <- numeric()
  Pred <- numeric()
  res.df <- data.frame(Status, Pred)
  
  for(i in 1:nrow(site.ns.summary)) {
    left.out.site <- site.ns.summary[i, 1]
    print(paste0('left out site: ', left.out.site))
    
    training.data <- full.df[which(full.df$ID != left.out.site), ]
    training.data <- training.data[sample((1:nrow(training.data)), size = round((nrow(training.data) / 4)), replace = FALSE), ]
    testing.data <- full.df[which(full.df$ID == left.out.site), ]
    testing.data <- testing.data[sample((1:nrow(testing.data)), size = round((nrow(testing.data) / 4)), replace = FALSE), ]
    
    class.mod <- lda(class.fun, data = training.data)
    res <- as.numeric(predict(class.mod, testing.data)$class)
    
    cur <- data.frame(as.numeric(testing.data$site), res)
    colnames(cur) <- c('Status', 'Pred')
    
    res.df <- rbind(res.df, cur)
    
  }
  
  
  a <- roc(res.df$Status, res.df$Pred, auc = TRUE)
  aucs[j] <- a$auc[1]
  
  print(paste0('AUC: ', aucs[j]))
  
  rm(res.df)
  
  print(paste0('Run ', j, ' of ', runs, ' finished'))
  
}

loocv.lda.aucs <- aucs

save(loocv.lda.aucs, file = 'loocv_lda_aucs.RData')

rm(list = ls())

wet <- raster('./CD_AOI_All_Wetlands/Irwin_CD_AOI_All_Wetlands.tif')
dem <- raster('./DEM_AOI/Irwin_DEM_AOI.tif')
lrm <- raster('./LRM_AOI/Irwin_LRM_AOI.tif')
sbs <- raster('./Site_Boundaries_Lithic_Scatter/Irwin_Prehistoric_NREI_Lithic_Scatter_AOI.tif')
slp <- raster('./Slope_AOI/Irwin_Slope_AOI.tif')
wv3 <- raster::stack('./Resampled_WV3/15SEP21184609-M3DS_R01C1-057029324010_01_P001_radiance_resample_AOI.tif')

site.1.ns.0 <- Which(sbs != 128)
cellnums <- site.1.ns.0
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
binary.coords <- predict(pcs, newdata = binary.stack.df[, 6:13])
binary.stack.df$pc.1 <- binary.coords[, 1]

site.cells <- which(binary.stack.df$site == 1)
ns.cells <- which(binary.stack.df$site == 0)

samp.size <- round(length(site.cells) / 4)

set.seed(2113)

aucs <- numeric(length = runs)

class.fun <- formula(site ~ wet + dem + lrm + slp + pc.1)

for(i in 1:runs) {
  
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
  
  pred.classes.lda <- binary.res.lda$class
  
  rand.samp.res.lda.df <- data.frame(testing.df$site, pred.classes.lda)
  colnames(rand.samp.res.lda.df) <- c('True', 'Predicted')
  rand.samp.res.lda.df$True <- as.numeric(rand.samp.res.lda.df$True) - 1
  rand.samp.res.lda.df$Predicted <- as.numeric(rand.samp.res.lda.df$Predicted) - 1
  
  rand.samp.roc.lda <- roc(rand.samp.res.lda.df$True, rand.samp.res.lda.df$Predicted, auc = TRUE, ci = TRUE)
  
  aucs[i] <- rand.samp.roc.lda$auc[1]
  
  print(aucs[i])
  
  print(paste0('Run ', i, ' of ', runs, ' finished'))
  
  
}

rs.lda.aucs <- aucs

save(rs.lda.aucs, file = 'rs_lda_aucs.RData')

rm(list = ls())

wet <- raster('./CD_AOI_All_Wetlands/Irwin_CD_AOI_All_Wetlands.tif')
dem <- raster('./DEM_AOI/Irwin_DEM_AOI.tif')
lrm <- raster('./LRM_AOI/Irwin_LRM_AOI.tif')
sbs <- raster('./Site_Boundaries_Lithic_Scatter/Irwin_Prehistoric_NREI_Lithic_Scatter_AOI.tif')
slp <- raster('./Slope_AOI/Irwin_Slope_AOI.tif')
wv3 <- raster::stack('./Resampled_WV3/15SEP21184609-M3DS_R01C1-057029324010_01_P001_radiance_resample_AOI.tif')

site.1.ns.0 <- Which(sbs != 128)
cellnums <- site.1.ns.0
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
binary.coords <- predict(pcs, newdata = binary.stack.df[, 6:13])
binary.stack.df$pc.1 <- binary.coords[, 1]

site.cells <- which(binary.stack.df$site == 1)
ns.cells <- which(binary.stack.df$site == 0)

samp.size <- round(length(site.cells) / 4)

set.seed(2113)

aucs <- numeric(length = 200)

class.fun <- formula(site ~ wet + dem + lrm + slp + pc.1)

runs <- 200

for(i in 1:runs) {
  
  training.site.cell.samp <- sample(1:length(site.cells), samp.size)
  training.ns.cell.samp <- sample(1:length(ns.cells), samp.size)
  
  training.site.cells <- site.cells[training.site.cell.samp]
  training.ns.cells <- ns.cells[training.ns.cell.samp]
  
  testing.site.cells <- site.cells[-training.site.cell.samp]
  testing.ns.cells <- sample(ns.cells[-training.ns.cell.samp], samp.size)
  
  training.df <- rbind(binary.stack.df[training.site.cells, ], binary.stack.df[training.ns.cells, ])
  testing.df <- rbind(binary.stack.df[testing.site.cells, ], binary.stack.df[testing.ns.cells, ])
  
  class.mod <- glm(class.fun, family = binomial, data = training.df)
  pred.classes.log <- predict(class.mod, testing.df, type = 'response')
  pred.classes.log[which(pred.classes.log >= 0.5)] <- 1
  pred.classes.log[which(pred.classes.log <= 0.5)] <- 0
  
  cur <- data.frame(as.numeric(testing.df$site), pred.classes.log)
  colnames(cur) <- c('Status', 'Pred')
  
  a <- roc(cur$Status, cur$Pred, auc = TRUE)
  
  aucs[i] <- a$auc[1]
  
  print(aucs[i])
  
  print(paste0('Run ', i, ' of ', runs, ' finished'))
  
}

rs.log.aucs <- aucs

save(rs.log.aucs, file = 'rs_log_aucs.RData')

