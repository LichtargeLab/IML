# 
# Drug response prediction
# This analysis predict IC50 values of PanCancer cell lines based on four single data type: 
# Gene Expression, Mutations, CNV, and Methylations
# 
# Author: Saeid Parvandeh Dec 2020
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
no_cores <- 64

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
drugs_names <- unique(pheno_data$DRUG_NAME)

data_types <- c("rnaseq", "mut", "cnv", "meth")
nfolds=c(10,10)
foreach(drug=1:length(drugs_names), .errorhandling = "pass") %dopar% {
  pheno_data_fltr <- pheno_data[which(pheno_data$DRUG_NAME==drugs_names[drug]), ]
  for (data_type in data_types){
    # setwd("...")
    drug_data <- get(load(paste0(data_type, "_data.RData")))
    common_id <- intersect(row.names(drug_data), pheno_data_fltr$SANGER_MODEL_ID)
    drug_data_fltr <- drug_data[rownames(drug_data)%in%common_id, ]
    # drug_data_fltr[is.na(drug_data_fltr)] <- 0
    pheno_data_fltr <- pheno_data_fltr[pheno_data_fltr$SANGER_MODEL_ID%in%common_id, ]
    pheno_data_fltr <- pheno_data_fltr[match(rownames(drug_data_fltr),  pheno_data_fltr$SANGER_MODEL_ID),]
    
    fold_data <- data.frame(drug_data_fltr, label=as.numeric(pheno_data_fltr$LN_IC50))
    outer_folds <- cvTools::cvFolds(nrow(fold_data), K=nfolds[1])
    
    nCV_atts <- list()
    en_pred_r2 <- NULL
    rf_pred_r2 <- NULL
    svm_pred_r2 <- NULL
    dt_pred_r2 <- NULL
    pred_r2 <- NULL
    folds_list_ <- list()
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
        
        # Statistical Significant
        vars_scores_fltr <- sort(vars_scores[which(vars_scores>0)], T)
        
        # Back-Transformation (square root)
        trans_fun <- function(x){x^(1/2)}
        
        vars_scores_trans <- trans_fun(vars_scores_fltr)
        
        # calculate the p-value for each gene
        Pvals <- NULL
        for (p in vars_scores_trans){
          Pvals <- c(Pvals, pnorm(p, mean(vars_scores_trans), sd(vars_scores_trans), lower.tail = F))
        }
        
        # calculate the q-value for each positive gene
        Qvals <- p.adjust(Pvals, method = "BH")
        
        genelist <- data.frame(gene=names(vars_scores_trans), pvalue=Pvals, qvalue=Qvals)
        
        genelist_sorted <- genelist[order(genelist$qvalue), ]
        atts[[j]] <-  genelist_sorted$gene[genelist_sorted$pvalue<0.05]
        atts_sig <-  rbind(atts_sig, genelist_sorted[genelist_sorted$pvalue<0.05, ])
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
    
    pred_r2 <- c(mean(en_pred_r2), mean(rf_pred_r2), mean(svm_pred_r2), mean(dt_pred_r2))
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
    
    cncv_res <- list(best_pred, best_obs_pred, nCV_atts_sig_sorted)
    # Save results
    setwd("...")
    save(cncv_res, file = paste0(drugs_names[drug], "_", data_type, "_res.RData"))
  }
}

# Finish
stopCluster(cl)
