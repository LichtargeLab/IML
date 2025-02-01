# 
# Integration of heterogeneous data for drug response prediction in Breast cancer
# This analysis predict IC50 based on data combination: Gene Expression, Methylation, CNV, and Mutation
# 

# Author: Saeid Parvandeh Oct 2022
# .libPaths("/storage/lichtarge/home/parvande/.conda/envs/r-env/lib/R/library")
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
no_cores <- 32

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

data_types <- c("rnaseq", "mut", "cnv", "meth")
dt_comb <- combn(data_types, 2)
nfolds=c(10,10)

foreach(dr=1:length(drugs_name), .errorhandling = "pass") %dopar% {
  pheno_data_ca_dr <- pheno_data_ca[pheno_data_ca$DRUG_NAME %in% drugs_name[dr], ]
  drug_vs <- unique(pheno_data_ca_dr$DRUG_ID)
  if (length(drug_vs) >= 2) {# if there is two versions of a drug
    pheno_data_ca_dr$LN_IC50 <- as.numeric(pheno_data_ca_dr$LN_IC50)
    pheno_data_ca_dr <- aggregate(pheno_data_ca_dr$LN_IC50, by=list(pheno_data_ca_dr$SANGER_MODEL_ID), FUN=mean)
    colnames(pheno_data_ca_dr) <- c("SANGER_MODEL_ID", "LN_IC50")
  }
  
  if (nrow(pheno_data_ca_dr)>=30){
    for (m in 1:ncol(dt_comb)){
      # setwd("...")
      dt <- dt_comb[, m]
      dt1 <- get(load(paste0(dt[1],"_data.RData")))
      dt2 <- get(load(paste0(dt[2],"_data.RData")))
      
      common_id <- intersect(intersect(pheno_data_ca_dr$SANGER_MODEL_ID, row.names(dt1)), row.names(dt2))
      dt1_fltr <- dt1[rownames(dt1)%in%common_id, ]
      dt2_fltr <- dt2[rownames(dt2)%in%common_id, ]
      
      pheno_data_ca_dr_fltr <- pheno_data_ca_dr[pheno_data_ca_dr$SANGER_MODEL_ID%in%common_id, ]
      pheno_data_ca_dr_fltr <- pheno_data_ca_dr_fltr[match(rownames(dt1_fltr),  pheno_data_ca_dr_fltr$SANGER_MODEL_ID),]
      
      fold_data1 <- data.frame(dt1_fltr, label=as.numeric(pheno_data_ca_dr_fltr$LN_IC50))
      outer_folds1 <- cvTools::cvFolds(nrow(fold_data1), K=nfolds[1])
      
      fold_data2 <- data.frame(dt2_fltr, label=as.numeric(pheno_data_ca_dr_fltr$LN_IC50))
      outer_folds2 <- cvTools::cvFolds(nrow(fold_data2), K=nfolds[1])
      
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
        outer_train1 <- fold_data1[which(outer_folds1[[5]]!=i), ]
        outer_test1 <- fold_data1[which(outer_folds1[[5]]==i), ]
        
        outer_train2 <- fold_data2[which(outer_folds2[[5]]!=i), ]
        outer_test2 <- fold_data2[which(outer_folds2[[5]]==i), ]
        
        inner_folds1 <- cvTools::cvFolds(nrow(outer_train1), K=nfolds[2])
        inner_folds2 <- cvTools::cvFolds(nrow(outer_train2), K=nfolds[2])
        
        atts1 <- list()
        atts2 <- list()
        atts_sig <- NULL
        for (j in 1:nfolds[1]){
          for (t in 1:length(dt)){
            assign("outer_train", get(noquote(paste0("outer_train", t))))
            assign("inner_folds", get(noquote(paste0("inner_folds", t))))
            inner_train <- outer_train[which(inner_folds[[5]]!=j), ]
            # Percentage of missing values for each gene
            pMiss <- function(x){sum(is.na(x))/length(x)*100}
            # removing columns with missing value
            inner_train_fltr <- inner_train[ ,which(apply(inner_train, 2, pMiss) != 100)]
            inner_train_fltr <- inner_train_fltr[which(apply(t(inner_train_fltr), 2, pMiss) != 100), ]
            
            # Nearest neighbor feature selection
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
            
            # calculate the q-value for each positive gene
            Qvals <- p.adjust(Pvals, method = "BH")
            
            genelist <- data.frame(gene=names(vars_scores_trans), pvalue=Pvals, qvalue=Qvals, type=dt[t])
            
            genelist_sorted <- genelist[order(genelist$pvalue), ]
            if (t == 1) {
              atts1[[j]] <-  genelist_sorted$gene#[genelist_sorted$pvalue<0.05]
              atts_sig <-  rbind(atts_sig, genelist_sorted)#[genelist_sorted$pvalue<0.05, ])
            } else if (t == 2) {
              atts2[[j]] <-  genelist_sorted$gene#[genelist_sorted$pvalue<0.05]
              atts_sig <-  rbind(atts_sig, genelist_sorted)#[genelist_sorted$pvalue<0.05, ])
            }
          }
        }
        nCV_atts1 <- Reduce(intersect, atts1)
        nCV_atts2 <- Reduce(intersect, atts2)      
        
        # Integrate multiple data types 
        # One possibility is that we could integrate all data types before feature selection, but
        # training data becomes large which tends to reduce processing time. Another reason for
        # integrating the data types after feature selection is because, it's possible that certain
        # data types might be preferred over others by a metric in a mixed high dimensional space.
        # And if there is a more of a certain data type, that data type could favored. So, I have
        # thought about an integrative metric for mixed data types. Further, re-training the data 
        # for selecting feature with more predicted power will remove noisy features.
        # Instead of re-training for feature selection, I used forward feature selection strategy.
        
        outer_train_fltr <- cbind(outer_train1[, nCV_atts1, drop=FALSE], outer_train2[, nCV_atts2, drop=FALSE])
        outer_train_fltr$label <- outer_train1[, 'label']
        outer_test_fltr <- cbind(outer_test1[, nCV_atts1, drop=FALSE], outer_test2[, nCV_atts2, drop=FALSE])
        
        # replace NAs with zero
        outer_train_fltr[is.na(outer_train_fltr)] <- 0
        outer_test_fltr[is.na(outer_test_fltr)] <- 0
        
        # Forward feature selection
        integrated_atts <- c(nCV_atts1, nCV_atts2)
        integrated_atts_sig <- atts_sig[which(atts_sig$gene%in%integrated_atts), ]
        integrated_atts_sig_sorted <- as.character(unique(integrated_atts_sig[order(integrated_atts_sig$pvalue), 'gene']))
        best_r2 = 0; best_features = integrated_atts_sig_sorted[1]
        if (length(integrated_atts_sig_sorted)>=2){
          for (f in integrated_atts_sig_sorted[1:length(integrated_atts_sig_sorted)]){
            if (!is.element(f, best_features) & is.element(f, colnames(outer_train_fltr))){
              best_features <- c(best_features, f)
              # Random forest
              rf_model <- randomForest::randomForest(outer_train_fltr[, best_features, drop=FALSE], outer_train_fltr[, "label"])
              rf_pred <- stats::predict(rf_model, outer_test_fltr[, best_features, drop=FALSE])
              rf_r2_cur <- stats::cor(rf_pred, outer_test1[, "label"])^2
              if (rf_r2_cur > best_r2 & !is.na(rf_r2_cur)){
                best_r2 <- rf_r2_cur
              } else {
                best_features <- best_features[! best_features %in% f]
              }
            }
          }
        }
        
        nCV_atts_sig <- rbind(nCV_atts_sig, atts_sig[which(atts_sig$gene%in%best_features), ])
        
        # Trim data to the size of outer atts
        outer_train_trim <- outer_train_fltr[, best_features]
        outer_test_trim <- outer_test_fltr[, best_features]
        
        # Classification with five classifiers
        
        # ElasticNet
        en_model <- glmnet::cv.glmnet(as.matrix(outer_train_trim), outer_train1[, "label"], alpha = 0.5)
        en_pred <- stats::predict(en_model, as.matrix(outer_test_trim), s = en_model$lambda.min)
        names(en_pred) <- row.names(outer_test1)
        en_pred_r2 <- c(en_pred_r2, stats::cor(en_pred, outer_test1[, "label"])^2)
        en_obs_pred <- rbind(en_obs_pred, cbind(en_pred, outer_test1[, "label"]))
        # Random forest
        rf_model <- randomForest::randomForest(outer_train_trim, outer_train1[, "label"])
        rf_pred <- stats::predict(rf_model, outer_test_trim)
        rf_pred_r2 <- c(rf_pred_r2, stats::cor(rf_pred, outer_test1[, "label"])^2)
        rf_obs_pred <- rbind(rf_obs_pred, cbind(rf_pred, outer_test1[, "label"]))
        #SVM
        svm_model <- e1071::svm(outer_train_trim, outer_train1[, "label"])
        svm_pred <- stats::predict(svm_model, outer_test_trim)
        svm_pred_r2 <- c(svm_pred_r2, stats::cor(svm_pred, outer_test1[, "label"])^2)
        svm_obs_pred <- rbind(svm_obs_pred, cbind(svm_pred, outer_test1[, "label"]))
        #dt
        dt_model <- rpart::rpart(label ~., data=data.frame(outer_train_trim, label=outer_train1[, "label"]))
        dt_pred <- stats::predict(dt_model, outer_test_trim)
        dt_pred_r2 <- c(dt_pred_r2, stats::cor(dt_pred, outer_test1[, "label"])^2)
        dt_obs_pred <- rbind(dt_obs_pred, cbind(dt_pred, outer_test1[, "label"]))
        
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
      
      icncv_res <- list(best_pred, best_obs_pred, nCV_atts_sig_sorted)
      # Save results
      # setwd("...")
      save(icncv_res, file = paste0(ca, "_", drugs_name[dr], "_", dt[1], "_", dt[2], "_ffs_res.RData"))
    }
  }
}

# Finish
stopCluster(cl)
