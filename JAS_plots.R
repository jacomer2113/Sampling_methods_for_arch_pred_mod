require(raster)
require(rgdal)
require(MASS)
require(pROC)
require(ggplot2)

setwd('~/Documents/CSRM/JAS/Data/')

load('loocv_lda_aucs.RData')
load('loocv_log_aucs.RData')
load('rs_lda_aucs.RData')
load('rs_log_aucs.RData')

df <- data.frame(loocv.lda.aucs, loocv.log.aucs, rs.lda.aucs, rs.log.aucs)

d.loocv.lda <- data.frame(loocv.lda.aucs)
d.loocv.log <- data.frame(loocv.log.aucs)
d.rs.lda <- data.frame(rs.lda.aucs)
d.rs.log <- data.frame(rs.log.aucs)

p1 <- ggplot(d.loocv.lda, aes(loocv.lda.aucs)) + geom_histogram(bins = 200) +
  geom_histogram(data = d.loocv.log, aes(loocv.log.aucs), bins = 200) +
  geom_histogram(data = d.rs.lda, aes(rs.lda.aucs), bins = 200) +
  geom_histogram(data = d.rs.log, aes(rs.log.aucs), bins = 200)

loocv.plot <- ggplot(d.loocv.lda, aes(loocv.lda.aucs)) + geom_histogram(bins = 75, fill = 'blue', alpha = 0.6) +
  geom_histogram(data = d.loocv.log, aes(loocv.log.aucs), bins = 75, fill = 'red', alpha = 0.6) +
  geom_vline(xintercept = mean(d.loocv.lda$loocv.lda.aucs), color = 'blue') +
  geom_vline(xintercept = mean(d.loocv.log$loocv.log.aucs), color = 'red') +
  labs(list(title = 'LOOCV', x = 'Area under the curve (AUC)', y = 'Count'))

rs.plot <- ggplot(d.rs.lda, aes(rs.lda.aucs)) + geom_histogram(bins = 75, fill = 'blue', alpha = 0.6) +
  geom_histogram(data = d.rs.log, aes(rs.log.aucs), bins = 75, fill = 'red', alpha = 0.6) +
  geom_vline(xintercept = mean(d.rs.lda$rs.lda.aucs), color = 'blue') +
  geom_vline(xintercept = mean(d.rs.log$rs.log.aucs), color = 'red') +
  labs(list(title = 'Random cell sample', x = 'Area under the curve (AUC)', y = 'Count'))
