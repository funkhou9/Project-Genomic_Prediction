#' **Script: 1-test_w_bglr.R**
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20151102*
#' - *Project: [Genomic_prediction](../../../README.md)*
#' - *Sub Folder: [wheat_dataset](../../wheat_dataset.md)*
#'
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 	- [Cluster and estimate composition](#cluster-and-estimate-composition)
#'	- [Fit gBLUP model](#fit-gblup-model)
#'
#' <br />
#' ## Objectives
#' Test methodology for the analysis of heterogeneous populations using genomic predictions.
#' Methods are intended to be a modification of those in
#' "Incorporating Genetic Heterogeneity in Whole-Genome Regressions Using Interactions".
#' by Gustavo de los Compos et al. This will involve:
#'
#' 1. Clustering wheat strains by PCA (as in the paper)
#' 2. Sampling individuals that represent the center of each cluster to build a reference panel of
#'	allele frequencies used in breed/strain composition estimation.
#' 3. Estimate composition of all samples using these references.
#' 4. Figure out a way to implement these breed composition estimates into the model proposed in the
#'	aforementioned paper.

setwd("/mnt/research/pigsnp/NSR/Genomic_prediction/wheat_dataset/scripts")

#' <br />
#' ## Install libraries
library(breedTools)
library(magrittr)
# devtools::install_github("gdlc/BGLR-R")
library(BGLR)
library(ggplot2)

sessionInfo()

#' <br />
#' ## Load data
#' Provided by BGLR package, a dataset consisting of 599 pure lines of wheat genotyped
#'	at 1279 DArT markers.
data(wheat)

#' Inspect wheat dataset components.
ls()
dim(wheat.X)

#' <br />
#' ## Analysis
#' ### Cluster and estimate composition
#' Center and scale X.
X <- scale(wheat.X, center = TRUE, scale = TRUE)
y <- wheat.Y[, 1]

#' Compute G.
G <- X %*% t(X) / ncol(X)

#' Use PCA to separate wheat groups as done previously.
evd <- eigen(G)
plot(evd$vectors[, 1],
	 evd$vectors[, 2],
	 xlab = "PC1",
	 ylab = "PC2")

#' By eye, arbitrarily define the individuals we want to isolate for "breed composition" reference.
#' Might there be a more precise way to "grap the peaks" from each?
group1_idx <- (evd$vectors[, 1] > -0.05 & evd$vectors[, 1] < -0.035 & evd$vectors[, 2] > -0.01 & evd$vectors[, 2] < 0.03)
group2_idx <- (evd$vectors[, 1] > 0.018 & evd$vectors[, 1] < 0.04 & evd$vectors[, 2] > -0.055 & evd$vectors[, 2] < -0.03)
plot(evd$vectors[, 1],
	 evd$vectors[, 2],
	 xlab = "PC1",
	 ylab = "PC2",
	 col = c("black", "red")[(group1_idx | group2_idx) + 1])


#' Allele frequency calculations with breedTools v0.1 requires rownames for genotype matrix X
rownames(wheat.X) <- seq(1, nrow(wheat.X))

#' Solve "breed composition" for each line.
bc <- breedTools::allele_freq(wheat.X, list(group1 = rownames(wheat.X[group1_idx, ]),
									  	  	group2 = rownames(wheat.X[group2_idx, ]))) %>%
		breedTools::solve_composition(wheat.X, .)

#' Obtain "purity" of group 1.
prob <- bc[, 1]

#' Re-plot and color according to "group purity".
pcs <- as.data.frame(evd$vectors)
colnames(pcs) <- paste(rep("PC", ncol(pcs)), seq(1:ncol(pcs)), sep = '')
pcs <- cbind(prob, pcs)

ggplot(pcs, aes(x = PC1, y = PC2, color = prob)) +
	geom_point(size = 4, alpha = 0.7) + scale_color_gradient(low = "red", high = "blue")

#' ### Fit gBLUP model
X0 <- X
X1 <- X
X2 <- X

for (i in 1:nrow(X1)) { 
  X1[i, ] <- X1[i, ] * prob[i]
  X2[i, ] <- X2[i, ] * (1 - prob[i])
}

#' Standardize each incidence matrix.
X0 <- scale(X0) / sqrt(ncol(X0))
X1 <- scale(X1) / sqrt(ncol(X1))
X2 <- scale(X2) / sqrt(ncol(X2))

#' Create linear predictor.
ETA <- list(main = list(X = X0, model = 'BRR', df0 = .1),
  		    int1 = list(X = X1, model = 'BRR', df0 = .1),
  		    int2 = list(X = X2, model = 'BRR', df0 = .1))

#' Fit the model.
#+ results = 'hide'
fm <- BGLR(y = y,
		   ETA = ETA,
		   nIter = 12000,
		   burnIn = 2000,
		   groups = ifelse(prob > .5, 1, 2), df0 = .1)

varU0 <- fm$ETA$main$varB
varU1 <- fm$ETA$int1$varB
varU2 <- fm$ETA$int2$varB
COR <- varU0 / sqrt((varU0 + varU1) * (varU0 + varU2))
fm$varE



