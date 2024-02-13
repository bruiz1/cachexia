setwd("C://Users//bryan//Kong Lab Dropbox//Bryan Ivan Ruiz//biruiz@uci.edu’s files//ct26 vs c26 all data//cholsoon lcms//R code for ct26 vs c26 figures")

### Install package dependencies  

metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()
install.packages("pacman")
pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", 
                 "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", 
                 "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest",
                 "RBGL","edgeR","fgsea","httr","qs"), 
               force = TRUE,
               character.only = TRUE)
library(pacman)
pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","httr","qs"))
# Install the package
# Step 1: Install devtools
install.packages("devtools")
library(devtools)
# Step 2: Install MetaboAnalystR with documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)


###Univariate Methods

# Load MetaboAnalystR
# Clean global environment
rm(list = ls())
library(MetaboAnalystR)
#InitDataObjects - Initializes an empty mSet object to hold the data
mSet<-InitDataObjects("conc", "stat", FALSE);
#Read.TextData - Reads in the data from the CSV file at the given URL. rowu indicates rows are sample IDs, disc means discrete outcome variables.
mSet<-Read.TextData(mSet, "https://www.xialab.ca/api/download/metaboanalyst/human_cachexia.csv", "rowu", "disc");
#SanityCheckData - Checks data integrity, removes entries with too many missing values, etc. Gets it ready for preprocessing.
mSet<-SanityCheckData(mSet);
#ReplaceMin - Replaces exceptionally small values (set as fraction of minimum positive value) with a small positive number to avoid taking log of zero or negative values later.
mSet<-ReplaceMin(mSet);
#PreparePrenormData - Additional data check after sanitizing.
mSet<-PreparePrenormData(mSet);
#Normalization - Actually normalizes the data using a log transform, mean-centering, and Pareto scaling. Specific preprocessing parameters are specified
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20);
#Normalization Plot
mSet<-PlotNormSummary(mSet, "norm_0_", format ="png", dpi=72, width=NA);
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "png", dpi=72, width=NA);



###Fold-change analysis

# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
# Plot fold-change analysis
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
# To view fold-change 
mSet$analSet$fc$fc.log



###T-Test

# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, "fdr", FALSE)
# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)



###Volcano Plot

# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")
# Create the volcano plot
mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="png", dpi=72, width=NA)



### One-way Analysis of Variance (ANOVA) - ONLY FOR MULTI-GROUP ANALYSIS!

# Perform ANOVA
mSet <- ANOVA.Anal(mSet, F, 0.05, "fisher")
# Plot ANOVA
mSet <- PlotANOVA(mSet, "aov_0_", "png", 72, width=NA)



### Correlation Analysis NON FUNCTIONAL AS OF 2.12.24

# OPTION 1 - Heatmap specifying pearson distance and an overview
mSet<-PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, 0.0)
# OPTION 2 - Heatmap specifying pearson correlation and a detailed view
mSet<-PlotCorrHeatMap(mSet, "corr_1_", format = "png", dpi=72, width=NA, "col", "spearman", "bwm", "detail", F, F, 999)



### Pattern Searching

# Perform correlation analysis on a pattern (a feature of interest in this case)
mSet<-FeatureCorrelation(mSet, "pearson", "1,6-Anhydro-beta-D-glucose")
# Plot the correlation analysis on a pattern
mSet<-PlotCorr(mSet, "ptn_3_", format="png", dpi=72, width=NA)



### Principal Component Analysis (PCA)

# Dependencies
install.packages("ellipse")
# Perform PCA analysis
mSet<-PCA.Anal(mSet)
# Create PCA overview
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "png", dpi = 72, width=NA, 5)
# Create PCA scree plot
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", dpi = 72, width=NA, 5)
# Create a 2D PCA score plot
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
# Create a 3D PCA score plot
mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
# Create a PCA loadings Plots
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
# Create a PCA Biplot
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "png", dpi = 72, width=NA, 1, 2)



### Partial Least Squares - Discriminant Analysis (PLS-DA)

mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2);
mSet<-PLSDA.CV(mSet, "5", 5,5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15, FALSE)
mSet<-PLSDA.Permut(mSet, 100, "accu")
mSet<-PlotPLS.Permutation(mSet, "pls_perm_1_", "png", 72, width=NA)
# View the 3D interactive PLS-DA score plot
mSet$imgSet$plsda.3d



### Sparse Partial Least Squares - Discriminant Analysis (sPLS-DA)

# Perform sPLS-DA analysis
mSet<-SPLSR.Anal(mSet, 5, 10, "same", "Mfold", 5, T)
# Plot sPLS-DA overview
mSet<-PlotSPLSPairSummary(mSet, "spls_pair_0_", format = "png", dpi=72, width=NA, 5)
# Create 2D sPLS-DA Score Plot
mSet<-PlotSPLS2DScore(mSet, "spls_score2d_0_", format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
# Create 3D sPLS-DA Score Plot
mSet<-PlotSPLS3DScoreImg(mSet, "spls_score3d_0_", format = "png", 72, width=NA, 1, 2, 3, 40)
# Create sPLS-DA loadings plot
mSet<-PlotSPLSLoading(mSet, "spls_loading_0_", format = "png", dpi=72, width=NA, 1,"overview")
# Perform cross-validation and plot sPLS-DA classification
mSet<-PlotSPLSDA.Classification(mSet, "spls_cv_0_", format = "png", dpi=72, width=NA)
# View the 3D interactive PLS-DA score plot
mSet$imgSet$splsda.3d



### Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)

# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)
# Create a 2D oPLS-DA score plot
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "png", dpi=72, width=NA, 1,2,0.95,1,0)
# Create a significant features plot
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
# Create a plot of features ranked by importance
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE)
# Create a plot of the model overview
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "png", dpi=72, width=NA)
# Perform and plot oPLS-DA permutation 
mSet<-OPLSDA.Permut(mSet, 100)
## [1] "time taken for 10 permutations:  0.0561895370483398"
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "png", dpi=72, width=NA)
## [1] "Empirical p-values Q2:  p < 0.01 (0/100)  and R2Y:  p = 0.01 (1/100)"



### Significance Analysis of Microarrary (and Metabolites) (SAM)

# Perform SAM analysis
mSet<-SAM.Anal(mSet, "d.stat", FALSE, TRUE, 0.0, "sam_imp_0_")
# Create a SAM plot of FDR values
mSet<-PlotSAM.FDR(mSet, "sam_view_0_", format = "png", dpi=72, width=NA)
# Create a SAM plot of results
mSet<-PlotSAM.Cmpd(mSet, "sam_imp_0_", format = "png", dpi=72, width=NA)



### Empirical Bayesian Analysis of Microarray (and Metabolites) (EBAM)

# Perform EBAM analysis, plot EBAM analysis and create the EBAM matrix of significant features
mSet<-EBAM.Init(mSet, FALSE, TRUE, FALSE, -99.0, 0.9, "ebam_view_0_", "ebam_imp_0_")
# Create a EBAM plot of results
PlotEBAM.Cmpd(mSet, "ebam_imp_0_", "png", 72, width=NA)



### Hierarchical Clustering: Dendogram

# Perform hierarchical clustering and plot dendogram
mSet<-PlotHCTree(mSet, "tree_0_", format = "png", dpi=72, width=NA, "euclidean", "ward.D")



### Partitional Clustering: K-Means

# Perform K-means analysis
mSet<-Kmeans.Anal(mSet, 3)
# Plot K-means analysis 
mSet<-PlotKmeans(mSet, "km_0_", format = "png", dpi=72, width=NA)



### Partitional Clustering: Self Organizing Maps (SOM)

# Perform SOM analysis
mSet<-SOM.Anal(mSet, 1, 3,"linear","gaussian")
# Plot SOM analysis
mSet<-PlotSOM(mSet, "som_0_", format = "png", dpi=72, width=NA)



### Random Forest

# Perform random forest analysis
mSet<-RF.Anal(mSet, 500, 7, 1)
# Plot random forest classification
mSet<-PlotRF.Classify(mSet, "rf_cls_0_", format = "png", dpi=72, width=NA)
# Plot random forest variables of importance
mSet<-PlotRF.VIP(mSet, "rf_imp_0_", format = "png", dpi=72, width=NA)
# Plot random forest outliers 
mSet<-PlotRF.Outlier(mSet, "rf_outlier_0_", format = "png", dpi=72, width=NA)



### Support Vector Machine (SVM)

# Perform SVM 
mSet<-RSVM.Anal(mSet, 10)
mSet<-PlotRSVM.Classification(mSet, "svm_cls_0_", format = "png", dpi=72, width=NA)
mSet<-PlotRSVM.Cmpd(mSet, "svm_imp_0_", format = "png", dpi=72, width=NA)


