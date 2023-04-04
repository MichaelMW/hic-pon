## Pre-requisite
# https://rdrr.io/bioc/HiTC/f/inst/doc/HiTC.pdf
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("HiTC")
# BiocManager::install("HiCDataHumanIMR90")

# install.packages("RCurl", type = "binary")
# install.packages("ggplot2")

library(HiTC)
require(Matrix)

## convert hic to pca A/B compartment:
hic_to_pca = function(data) {
  # parameters
  bin_size = 50000
  num_pc = 1
  
  # process hi-c data
  data_chrX = extractRegion(data$chrXchrX, chr="chrX", from=0, to=1000000000)
  data_chrX = binningC(data_chrX, binsize=bin_size, method="median", step=3)
  
  # get PCA
  pc = pca.hic(data_chrX, normPerExpected=TRUE, method="loess", npc=num_pc)
  score = score(pc$PC1)
  start = start(pc$PC1)
  score[is.na(score)] = 0
  return(list(score=score, start = start, data=data_chrX))
}

## data:
data(Nora_5C)
# E14     HiTC - 5C data    male undifferentiated ES cells -> use this as our "normal signal"
# MEF     HiTC - 5C data.   male embryonic fibroblasts  -> use this as our "tumor signal"

## E14, this is our base case, the underlying "normal" signal
result_E14 = hic_to_pca(E14)
score_normal = result_E14$score
start = result_E14$start

## build a panel of normal with E14 with random noise
n_sample = 10
noise_level = 0.00001
noise_matrix = matrix(rnorm(n_sample * length(score_normal), mean = 0, sd = noise_level), ncol = length(score_normal), byrow = TRUE)
panel_of_normal = matrix(NA, ncol = length(score_normal), nrow = n_sample)
for (i in 1:n_sample) {
  panel_of_normal[i,] = score_normal + noise_matrix[i,]
}

## MEF, this is our "pure tumor" case. 
result_MEF = hic_to_pca(MEF)
score_tumor = result_MEF$score

## mixture case: E14 + MEF * spike-in ratio
tumor_fraction = 0.01 # 1% 
score_mixture = score_normal + score_tumor * tumor_fraction

## Build a linear regression model from panel of normal:
pca = prcomp(panel_of_normal, rank.=5)
data_pca = data.matrix(pca$x)
reg = lm(panel_of_normal ~ data_pca)

# Project the new sample onto the principal components
new_sample = score_mixture
new_sample_pca = predict(pca, newdata=t(new_sample))
predicted_score = predict(reg, newdata=list(data_pca=new_sample_pca))

# Calculate residuals as final normalized values for the given sample
residuals = new_sample - predicted_score

# visualizing the original HiC data
mapC(HTClist(result_MEF$data), maxrange=100)
mapC(HTClist(result_E14$data), maxrange=100)

# visualizing the PCA normalization
dev.off()
par(mfrow = c(2,1))
par(mar=c(5,5,2,2))
plot(start, score_tumor, type="h", xlab ="", ylab="", main="Tumor", col = "red")
plot(start, score_normal, type="h", xlab ="",  ylab="", main="Normal", col = "blue")
plot(start, score_mixture, type="h", xlab ="", ylab="", main="Normal + 1% Tumor", col = "purple")
plot(start, residuals, type="h", xlab ="", ylab="", main="Residuals", col = "red")
