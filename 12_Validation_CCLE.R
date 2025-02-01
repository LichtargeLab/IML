#
# Validate GDSC data on CCLE -- RNA-Seq
#
# Author: Saeid Parvandeh, September 2023

# GDSC
# setwd("...")
gdsc_pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
gdsc_drugs_name <- unique(gdsc_pheno_data$DRUG_NAME)
gdsc_rnaseq_data <- get(load("rnaseq_data.RData"))


# CCLE
setwd("~/Documents/Drug Response/IML/Data/CCLE")
ccle_pheno_data <- get(load('ccle_pheno_data.RData'))
ccle_rnaseq_data <- get(load('ccle_rnaseq_data.RData'))

ccle_drugs_name <- unique(ccle_pheno_data$Compound)

intersect(gdsc_drugs_name, ccle_drugs_name)

ccle_drugs_name <- c("Nilotinib", "PHA-665752", "Lapatinib", "Sorafenib", "Paclitaxel", 
                     "Erlotinib", "Panobinostat", "Irinotecan", "Topotecan", 'Nutlin-3a (-)')

# RNASeq
for (drug in ccle_drugs_name){
  # setwd("...")
  res <- get(load(paste0(drug, '_rnaseq_res.RData')))
  res_genes <- as.character(res[[3]][1]$gene)
  res_genes <- gsub('\\.', '-', res_genes)
  
  # GDSC
  gdsc_pheno_fltr <- gdsc_pheno_data[gdsc_pheno_data$DRUG_NAME == drug, ]
  gdsc_cls <- intersect(rownames(gdsc_rnaseq_data), unique(gdsc_pheno_fltr$SANGER_MODEL_ID))
  res_genes[which(res_genes=="SEPTIN4-1")] <- "SEPTIN4.1"
  gdsc_rnaseq_fltr <- gdsc_rnaseq_data[rownames(gdsc_rnaseq_data) %in% gdsc_cls, res_genes]
  gdsc_pheno_fltr <- gdsc_pheno_fltr[gdsc_pheno_fltr$SANGER_MODEL_ID %in% gdsc_cls, ]
  gdsc_pheno_fltr <- gdsc_pheno_data[match(rownames(gdsc_rnaseq_fltr),  gdsc_pheno_data$SANGER_MODEL_ID),]
  
  # CCLE
  if (drug == 'Nutlin-3a (-)'){
    drug = 'Nutlin-3'
  }
  
  ccle_pheno_fltr <- ccle_pheno_data[ccle_pheno_data$Compound == drug, c('DepMap_ID', 'CCLE_Name', 'IC50..uM.', 
                                                                         'ActArea', 'EC50..uM.', 'Amax')]
  
  ccle_cls <- intersect(colnames(ccle_rnaseq_data), unique(ccle_pheno_fltr$CCLE_Name))
  ccle_rnaseq_fltr <- ccle_rnaseq_data[which(ccle_rnaseq_data$Description%in%res_genes), c('Description', ccle_cls)]
  # removing duplicate genes
  ccle_rnaseq_fltr <- ccle_rnaseq_fltr[order(ccle_rnaseq_fltr$Description, -rowSums(ccle_rnaseq_fltr[, 2:ncol(ccle_rnaseq_fltr)]) ), ] #sort by description and reverse rowsum
  ccle_rnaseq_fltr <- ccle_rnaseq_fltr[ !duplicated(ccle_rnaseq_fltr$Description), ] # take the first row within each description
  # use first column for row names and transpose
  ccle_rnaseq_fltr <- data.frame(ccle_rnaseq_fltr, row.names = 1)
  ccle_rnaseq_fltr <- t(ccle_rnaseq_fltr)
  # remove X beginning of cell lines name
  rownames(ccle_rnaseq_fltr) <- gsub("X","",rownames(ccle_rnaseq_fltr))
  
  ccle_pheno_fltr <- ccle_pheno_fltr[ccle_pheno_fltr$CCLE_Name %in% ccle_cls, ]
  ccle_pheno_fltr <- ccle_pheno_fltr[match(rownames(ccle_rnaseq_fltr),  ccle_pheno_fltr$CCLE_Name),]
  
  fold_data <- data.frame(ccle_rnaseq_fltr, label=as.numeric(log(ccle_pheno_fltr$IC50..uM.)))
  
  
  
  # Train on GDSC
  res_common_genes <- intersect(res_genes, colnames(ccle_rnaseq_fltr))
  
  # Random forest
  rf_model <- randomForest::randomForest(gdsc_rnaseq_fltr[, res_common_genes], as.numeric(gdsc_pheno_fltr$LN_IC50))
  rf_pred <- stats::predict(rf_model, ccle_rnaseq_fltr)
  print(drug)
  print(summary(lm((rf_pred~ as.numeric(-ccle_pheno_fltr$ActArea))))$r.squared^(1/2))
  
  # # ElasticNet
  # en_model <- glmnet::cv.glmnet(as.matrix(gdsc_rnaseq_fltr[, res_common_genes]), as.numeric(gdsc_pheno_fltr$LN_IC50), alpha = 0.5)
  # en_pred <- stats::predict(en_model, as.matrix(ccle_rnaseq_fltr), s = en_model$lambda.min)
  # names(en_pred) <- row.names(ccle_rnaseq_fltr)
  # stats::cor(en_pred, as.numeric(-ccle_pheno_fltr$ActArea))
  # summary(lm((en_pred~ as.numeric(-ccle_pheno_fltr$ActArea))))$r.squared^(1/2)
  # # en_obs_pred <- rbind(en_obs_pred, cbind(en_pred, outer_test[, "label"]))
  # 
  # # SVM
  # svm_model <- e1071::svm(gdsc_rnaseq_fltr[, res_common_genes], as.numeric(gdsc_pheno_fltr$LN_IC50))
  # svm_pred <- stats::predict(svm_model, ccle_rnaseq_fltr)
  # stats::cor(svm_pred, as.numeric(-ccle_pheno_fltr$ActArea))
  # summary(lm((svm_pred~ as.numeric(-ccle_pheno_fltr$ActArea))))$r.squared^(1/2)
  # # svm_obs_pred <- rbind(svm_obs_pred, cbind(svm_pred, as.numeric(ccle_pheno_fltr$LN_IC50)))
  # 
  # #dt
  # dt_model <- rpart::rpart(label ~., data=data.frame(gdsc_rnaseq_fltr[, res_common_genes], label=as.numeric(gdsc_pheno_fltr$LN_IC50)))
  # dt_pred <- stats::predict(dt_model, data.frame(ccle_rnaseq_fltr))
  # stats::cor(dt_pred, as.numeric(-ccle_pheno_fltr$ActArea))
  # summary(lm((dt_pred~ as.numeric(-ccle_pheno_fltr$ActArea))))$r.squared^(1/2)
  # # dt_obs_pred <- rbind(dt_obs_pred, cbind(dt_pred, as.numeric(ccle_pheno_fltr$LN_IC50)))
  
}
