# 
# Drug response prediction
# This analysis predict IC50 values of individual cancer cell lines based on four single data type: 
# Gene Expression, Mutations, CNV, and Methylations
# 
# Author: Saeid Parvandeh Dec 2022
library(parallel)
library(foreach)
library(doParallel)
library(cvTools)
library(CORElearn)
library(glmnet)
library(randomForest)
library(e1071)
library(rpart)

# Calculate the number of cores
no_cores <- 4

# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# setwd("...")
# Gene Expression
load("rnaseq_data.RData")

# Mutation
load("mut_data.RData")
# drug_data <- mutation_data

# CNV
load("cnv_data.RData")

# Methylation
load("meth_data.RData")

pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
all_individual_tumors <- unique(pheno_data$TCGA_DESC[pheno_data$TCGA_DESC != ""])

cancers <- NULL
for (ca in all_individual_tumors){
  pheno_data_ca <- pheno_data[pheno_data$TCGA_DESC == ca, ]
  if (length(unique(pheno_data_ca$SANGER_MODEL_ID)) >= 19 &
      length(unique(pheno_data_ca$SANGER_MODEL_ID)) <= 30){
    cancers <- c(cancers, ca)
  }
}
cancers <- cancers[which(cancers!="UNCLASSIFIED")]

cancers <- c("COREAD", "HNSC", "LUAD", "SKCM")
ca <- cancers[3]
pheno_data_ca <- pheno_data[pheno_data$TCGA_DESC %in% ca, ]

# Prioritizing the drugs
pheno_data_ca$Z_SCORE <- as.numeric(pheno_data_ca$Z_SCORE)
avg_dr <- NULL
for (dr in unique(pheno_data_ca$DRUG_NAME)){
  pheno_data_ca_dr <- pheno_data_ca[pheno_data_ca$DRUG_NAME == dr, ]
  pheno_data_ca_dr <- pheno_data_ca_dr[order(pheno_data_ca_dr$Z_SCORE), ]
  avg_dr <- c(avg_dr, mean(pheno_data_ca_dr$Z_SCORE[1:(nrow(pheno_data_ca_dr)/3)]))
}
names(avg_dr) <- unique(pheno_data_ca$DRUG_NAME)
# top sensitive drugs based on average drug response across cell lines
drugs_name <- names(sort(avg_dr))[1:50]

nfolds=c(10,10)
data_types <- c("rnaseq", "mut", "cnv", "meth")
foreach(dr=1:length(drugs_name), .errorhandling = "pass") %dopar% {
  pheno_data_ca_dr <- pheno_data_ca[pheno_data_ca$DRUG_NAME %in% drugs_name[dr], ]
  drug_vs <- unique(pheno_data_ca_dr$DRUG_ID)
  if (length(drug_vs) >= 2) {# if there is two versions of a drug
    pheno_data_ca_dr$LN_IC50 <- as.numeric(pheno_data_ca_dr$LN_IC50)
    pheno_data_ca_dr <- aggregate(pheno_data_ca_dr$LN_IC50, by=list(pheno_data_ca_dr$SANGER_MODEL_ID), FUN=mean)
    colnames(pheno_data_ca_dr) <- c("SANGER_MODEL_ID", "LN_IC50")
  }

  if (nrow(pheno_data_ca_dr)>=30){
    for (dt in data_types){
      setwd("...")
      drug_data <- get(load(paste0(dt, "_data.RData")))
      drug_data_fltr <- drug_data[rownames(drug_data) %in% pheno_data_ca_dr$SANGER_MODEL_ID, ]
      # Pre-filtering:  
      # (1) At least 75% of cell lines > 0
      # (2) RNA average > 1
      # pre_fltr_idx <- which(colSums(drug_data_fltr==0)/nrow(drug_data_fltr)*100 > 75 |
      #                         colMeans(drug_data_fltr, na.rm = T) < 1)
      # drug_data_fltr <- drug_data_fltr[, -pre_fltr_idx]
      pheno_data_ca_dr_fltr <- pheno_data_ca_dr[match(rownames(drug_data_fltr),  pheno_data_ca_dr$SANGER_MODEL_ID),]
      
      fold_data <- data.frame(drug_data_fltr, label=as.numeric(pheno_data_ca_dr_fltr$LN_IC50))
      outer_folds <- cvTools::cvFolds(nrow(fold_data), K=nfolds[1])
      
      nCV_atts <- list()
      en_pred_r2 <- NULL
      rf_pred_r2 <- NULL
      svm_pred_r2 <- NULL
      dt_pred_r2 <- NULL
      pred_r2 <- NULL
      en_obs_pred <- NULL
      rf_obs_pred <- NULL
      svm_obs_pred <- NULL
      dt_obs_pred <- NULL
      nCV_atts_sig <- NULL
      for (i in 1:nfolds[1]){
        outer_train <- fold_data[which(outer_folds[[5]]!=i), ]
        outer_test <- fold_data[which(outer_folds[[5]]==i), ]
        inner_folds <- cvTools::cvFolds(nrow(outer_train), K=nfolds[2])
        atts <- list()
        atts_sig <- NULL
        folds_list <- list()
        for (j in 1:nfolds[1]){
          inner_train <- outer_train[which(inner_folds[[5]]!=j), ]
          # Percentage of missing values for each gene
          pMiss <- function(x){sum(is.na(x))/length(x)*100}
          # removing columns with missing value
          inner_train_fltr <- inner_train[ ,which(apply(inner_train, 2, pMiss) != 100)]
          inner_train_fltr <- inner_train_fltr[which(apply(t(inner_train_fltr), 2, pMiss) != 100), ]
          # Relief
          vars_scores <- CORElearn::attrEval("label",
                                             inner_train_fltr,
                                             estimator = "RReliefFequalK",
                                             kNearestEqual = floor((dim(inner_train_fltr)[1]-1)*0.154))
          
          # Statistical Significants
          vars_scores_fltr <- sort(vars_scores[which(vars_scores>0)], T)
          
          # Back-Transformation (square root)
          trans_fun <- function(x){x^(1/2)}
          
          vars_scores_trans <- trans_fun(vars_scores_fltr)
          
          # calculate the p-value for each gene
          Pvals <- NULL
          for (p in vars_scores_trans){
            Pvals <- c(Pvals, pnorm(p, mean(vars_scores_trans), sd(vars_scores_trans), lower.tail = F))
          }
          
          genelist <- data.frame(gene=names(vars_scores_trans), pvalue=Pvals)
          
          genelist_sorted <- genelist[order(genelist$pvalue), ]
          atts[[j]] <-  genelist_sorted$gene#[genelist_sorted$pvalue<0.05]
          atts_sig <-  rbind(atts_sig, genelist_sorted)#[genelist_sorted$pvalue<0.05, ])
        }
        nCV_atts[[i]] <- Reduce(intersect, atts)
        nCV_atts_sig <- rbind(nCV_atts_sig, atts_sig[which(atts_sig$gene%in%nCV_atts[[i]]), ])
        
        # Classification with five classifiers
        # reduce columns to consensus features
        outer_train_fltr <- outer_train[, nCV_atts[[i]], drop=FALSE]
        outer_test_fltr <- outer_test[, nCV_atts[[i]], drop=FALSE]
        
        # replace NAs with zero
        outer_train_fltr[is.na(outer_train_fltr)] <- 0
        outer_test_fltr[is.na(outer_test_fltr)] <- 0
        
        # ElasticNet
        en_model <- glmnet::cv.glmnet(as.matrix(outer_train_fltr), outer_train[, "label"], alpha = 0.5)
        en_pred <- stats::predict(en_model, as.matrix(outer_test_fltr), s = en_model$lambda.min)
        names(en_pred) <- row.names(outer_test)
        en_pred_r2 <- c(en_pred_r2, stats::cor(en_pred, outer_test[, "label"])^2)
        en_obs_pred <- rbind(en_obs_pred, cbind(en_pred, outer_test[, "label"]))
        # Random forest
        rf_model <- randomForest::randomForest(outer_train_fltr, outer_train[, "label"])
        rf_pred <- stats::predict(rf_model, outer_test_fltr)
        rf_pred_r2 <- c(rf_pred_r2, stats::cor(rf_pred, outer_test[, "label"])^2)
        rf_obs_pred <- rbind(rf_obs_pred, cbind(rf_pred, outer_test[, "label"]))
        #SVM
        svm_model <- e1071::svm(outer_train_fltr, outer_train[, "label"])
        svm_pred <- stats::predict(svm_model, outer_test_fltr)
        svm_pred_r2 <- c(svm_pred_r2, stats::cor(svm_pred, outer_test[, "label"])^2)
        svm_obs_pred <- rbind(svm_obs_pred, cbind(svm_pred, outer_test[, "label"]))
        #dt
        dt_model <- rpart::rpart(label ~., data=data.frame(outer_train_fltr, label=outer_train[, "label"]))
        dt_pred <- stats::predict(dt_model, outer_test_fltr)
        dt_pred_r2 <- c(dt_pred_r2, stats::cor(dt_pred, outer_test[, "label"])^2)
        dt_obs_pred <- rbind(dt_obs_pred, cbind(dt_pred, outer_test[, "label"]))
        
      }
      pred_r2 <- c( mean(en_pred_r2), mean(rf_pred_r2), mean(svm_pred_r2), mean(dt_pred_r2))
      obs_pred <- list(en_obs_pred, rf_obs_pred, svm_obs_pred, dt_obs_pred)
      names(pred_r2) <- c("en", "rf", "svm", "dt")
      best_pred <- pred_r2[which.max(pred_r2)]
      best_obs_pred <- obs_pred[[which.max(pred_r2)]]
      
      # consensus genes
      library(dplyr)
      nCV_atts_sig_fltr <- nCV_atts_sig %>%
        group_by(gene) %>%
        filter(pvalue == min(pvalue)) %>%
        ungroup()
      nCV_atts_sig_sorted <- unique(nCV_atts_sig_fltr[order(nCV_atts_sig_fltr$pvalue), ])
      # FDR correction for multiple hypothesis
      nCV_atts_sig_sorted$qvalue <- p.adjust(nCV_atts_sig_sorted$pvalue, method = "BH")
      cncv_res <- list(best_pred, best_obs_pred, nCV_atts_sig_sorted)
      # Save results
      # setwd("...")
      save(cncv_res, file = paste0(ca, "_", drugs_name[dr], "_", dt, "_res.RData"))
    }
  }
}


# Finish
stopCluster(cl)

