#
# Post cancer specific analysis of selected drugs response
#
# --- drug:gene interactions --- #
# Author: Saeid Parvandeh, Jan 2024

# setwd("...")
# Gene Expression
load("rnaseq_data.RData")


wish_list <- cbind(c("COREAD", "COREAD", "COREAD", "COREAD", "COREAD"),
                   c("5-Fluorouracil", "Cetuximab", "Trametinib", "Oxaliplatin", "Irinotecan"))

cancers <- unique(wish_list[, 1])

SAMPLE_IDs_Class <- NULL
for (ca in cancers){
  for (pheno in c("GDSC1", "GDSC2")){
    # setwd("...")
    pheno_data <- read.csv(paste0(pheno, "_fitted_dose_response_25Feb20.csv"), colClasses = "character")
    pheno_data_ca <- pheno_data[pheno_data$TCGA_DESC %in% ca, ]
    drugs_name <- wish_list[wish_list[, 1] %in% ca, 2]
    
    for (dr in drugs_name) {
      pheno_data_ca_dr <- pheno_data_ca[pheno_data_ca$DRUG_NAME %in% dr, ]
      drug_vs <- unique(pheno_data_ca_dr$DRUG_ID)
      for (vs in drug_vs){
        pheno_data_ca_dr_vs <- pheno_data_ca_dr[which(pheno_data_ca_dr$DRUG_ID == vs &
                                                        pheno_data_ca_dr$SANGER_MODEL_ID %in% rownames(gdsc_rnaseq)), ]

        
        file = paste0("...", ca, "/", pheno, "_", ca, "_", dr, "_", dt, "_res_", vs, ".RData")
        path = file.path(file)
        if(!file.exists(path)) {
          next
        }
        IML_res <- get(load(path))
        
        # extract genes
        IML_genes <-IML_res[[3]]$gene#[which(IML_res[[3]]$pvalue < 0.05)]
        IML_genes <- gsub("\\.","-",IML_genes)
        IML_genes[which(IML_genes=="SEPTIN4-1")] <- "SEPTIN4.1"
        IML_res[[3]]$gene <- IML_genes
        
        common_id <- intersect(rownames(IML_res[[2]]), rownames(gdsc_rnaseq))
        gdsc_rnaseq_fltr <- gdsc_rnaseq[rownames(gdsc_rnaseq)%in%common_id, ]
        gdsc_rnaseq_fltr <- gdsc_rnaseq_fltr[match(rownames(IML_res[[2]]),  rownames(gdsc_rnaseq_fltr)),]
        
        
        # IC50 RNAseq ratio in sensitive vs. resistance cell lines
        if(nrow(pheno_data_ca_dr_vs) >= 18 & nrow(pheno_data_ca_dr_vs) <= 22){
          hi_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = F), ])[1:8]
          lo_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = T), ])[1:8]
        } else if(nrow(pheno_data_ca_dr_vs) > 22 & nrow(pheno_data_ca_dr_vs) <= 30){
          hi_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = F), ])[1:10]
          lo_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = T), ])[1:10]
        } else if(nrow(pheno_data_ca_dr_vs) > 30 & nrow(pheno_data_ca_dr_vs) <= 40){
          hi_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = F), ])[1:12]
          lo_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = T), ])[1:12]
        } else if(nrow(pheno_data_ca_dr_vs) > 40 & nrow(pheno_data_ca_dr_vs) <= 52){
          hi_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = F), ])[1:15]
          lo_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = T), ])[1:15]
        } else if(nrow(pheno_data_ca_dr_vs) > 52){
          hi_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = F), ])[1:16]
          lo_responders <- rownames(IML_res[[2]][order(IML_res[[2]][, 2], decreasing = T), ])[1:16]
        }
        # common_genes <- intersect(genes, colnames(data_fltr))
        rnaseq_sen <- gdsc_rnaseq_fltr[hi_responders, IML_genes]
        rnaseq_res <- gdsc_rnaseq_fltr[lo_responders, IML_genes]
        
        rna_avg_sen <- apply(rnaseq_sen, 2, mean, na.rm=TRUE)
        rna_avg_res <- apply(rnaseq_res, 2, mean, na.rm=TRUE)
        
        rnaseq_sen_res <- as.data.frame(rbind(rnaseq_sen, rnaseq_res))
        if(nrow(pheno_data_ca_dr_vs) >= 18 & nrow(pheno_data_ca_dr_vs) <= 22){
          rnaseq_sen_res$label <- c(rep(0, 8), rep(1, 8))
        } else if(nrow(pheno_data_ca_dr_vs) > 22 & nrow(pheno_data_ca_dr_vs) <= 30){
          rnaseq_sen_res$label <- c(rep(0, 10), rep(1, 10))
        } else if(nrow(pheno_data_ca_dr_vs) > 30 & nrow(pheno_data_ca_dr_vs) <= 40){
          rnaseq_sen_res$label <- c(rep(0, 12), rep(1, 12))
        } else if(nrow(pheno_data_ca_dr_vs) > 40 & nrow(pheno_data_ca_dr_vs) <= 52){
          rnaseq_sen_res$label <- c(rep(0, 15), rep(1, 15))
        } else if(nrow(pheno_data_ca_dr_vs) > 52){
          rnaseq_sen_res$label <- c(rep(0, 16), rep(1, 16))
        }
        rnaseq_ttest_list <- lapply(rnaseq_sen_res[-ncol(rnaseq_sen_res)], function(x) t.test(x ~ rnaseq_sen_res$label))
        rnaseq_ttest_pval <- NULL
        for (x in 1:length(rnaseq_ttest_list)){
          rnaseq_ttest_pval <- c(rnaseq_ttest_pval, rnaseq_ttest_list[[x]]$p.value)
        }
        names(rnaseq_ttest_pval) <- IML_genes
        
        rnaseq_measures <- data.frame(gene = IML_genes, ttest_p_value = rnaseq_ttest_pval, ratio = rna_avg_sen/rna_avg_res)
        
        # Average all cell lines
        rnaseq_avg <- apply(gdsc_rnaseq_fltr, 2, mean, na.rm=TRUE)
        IML_genes_pval <- as.data.frame(IML_res[[3]])
        drug_gene_int <- NULL
        for (g in IML_genes){
          if (all(is.na(gdsc_rnaseq_fltr[, g]))) {
            coeff <- NA; pval <- NA; r_squared <- NA
          } else {
            lm_fit <- lm(gdsc_rnaseq_fltr[, g] ~ IML_res[[2]][, 1])
            if (!is.na(coefficients(summary(lm_fit))[2]) & dim(coefficients(summary(lm_fit)))[1]==2) {
              coeff <- coefficients(summary(lm_fit))[2, 1]
              pval <- coefficients(summary(lm_fit))[2, 4]
              r_squared <- summary(lm_fit)$r.squared
            } else {
              coeff <- NA; pval <- NA; r_squared <- NA
            }
          }
          drug_gene_int <- rbind(drug_gene_int, data.frame(IML_genes_pval[IML_genes_pval$gene==g, ],
                                                           avg_sen = rna_avg_sen[g], avg_res = rna_avg_res[g], rna_ratio = rnaseq_measures[g, 'ratio'],
                                                           ttest_pval = rnaseq_measures[g, 'ttest_p_value'],
                                                           standardized_coefficinet = (sd(gdsc_rnaseq_fltr[, g])/sd(IML_res[[2]][, 1]))*coeff,
                                                           rna_avg = rnaseq_avg[g], lin_reg_pval = pval, lin_reg_r2=r_squared, effect_size=r_squared/(1-r_squared)))
        }
        
        # sum up the sensitive and resistant cell lines RNAs in one column
        drug_gene_int$SUM_AVG_SEN_RES <- rowSums(drug_gene_int[, c('avg_sen', 'avg_res')])
        # remove all genes with average expression < 0.5
        drug_gene_int_fltr <- drug_gene_int[which(drug_gene_int$rna_avg > 0.5), ] # Cutoff changes from 0.5/1 to 1/0.5 to select more/less genes
        # Separate sensitive and resistant genes
        #                                                                                      # Relaxing/strengthen the cutoff to select more/less genes
        drug_gene_int_sen <- drug_gene_int_fltr[which(drug_gene_int_fltr$rna_ratio > 1.25), ]  # 1.5
        drug_gene_int_res <- drug_gene_int_fltr[which(drug_gene_int_fltr$rna_ratio < 0.80), ]  # 0.67
        
        
        drug_gene_int_sen$GROUP <- rep('Sensitive', nrow(drug_gene_int_sen))
        drug_gene_int_res$GROUP <- rep('Resistant', nrow(drug_gene_int_res))
        drug_gene_int_sen_res <- rbind(drug_gene_int_sen, drug_gene_int_res)
        
        # setwd("...")
         write.table(drug_gene_int_sen_res, paste0(pheno, "_", ca, "_", dr, "_", vs, "_", dt, "_drug_gene_interaction_Larry.tsv"), quote = F, row.names = F, sep = "\t")
        
        # Extract sample IDs
        train_ids <- rownames(IML_res[[2]])
        test_ids <- setdiff(unique(pheno_data_ca_dr_vs$SANGER_MODEL_ID), rownames(IML_res[[2]]))
        
        pheno_data_ca_dr_vs$Z_SCORE <- as.numeric(pheno_data_ca_dr_vs$Z_SCORE)
        pheno_data_ca_dr_vs <- aggregate(pheno_data_ca_dr_vs$Z_SCORE, by=list(pheno_data_ca_dr_vs$SANGER_MODEL_ID), FUN=mean)
        colnames(pheno_data_ca_dr_vs) <- c("SANGER_MODEL_ID", "Z_SCORE")
        SAMPLE_IDs_Class <- rbind.data.frame(SAMPLE_IDs_Class, cbind.data.frame(CANCER = ca, DRUG = dr, DRUD_ID = vs,
                                                                                CLASS = c(rep("TRAIN", length(train_ids)), 
                                                                                          rep("TEST", length(test_ids))), 
                                                                                SAMPLE_IDs = c(train_ids, test_ids),
                                                                                Z_SCORES = c(IML_res[[2]][, 2], pheno_data_ca_dr_vs$Z_SCORE[pheno_data_ca_dr_vs$SANGER_MODEL_ID%in%test_ids]),
                                                                                IML_PREDICTION = c(IML_res[[2]][, 1], rep(NA, length(test_ids)))))
        
        # genes
        sen_res_genes <- rbind(drug_gene_int_sen[1:500, c("gene", "GROUP")], drug_gene_int_res[1:500, c("gene", "GROUP")])
        sen_res_genes_fltr <- sen_res_genes[!is.na(sen_res_genes$gene), ]
        sen_res_genes_fltr <- data.frame(sen_res_genes_fltr, row.names = 1)
        # rnaseq
        gdsc_rnaseq_ <- t(gdsc_rnaseq[c(train_ids, test_ids), rownames(sen_res_genes_fltr)])
        gdsc_rnaseq_[which(gdsc_rnaseq_ == 0)] <- 0.01
        gdsc_rnaseq_log2 <- log2(gdsc_rnaseq_)
        gdsc_rnaseq_norm <- t(apply(gdsc_rnaseq_log2, 1, scale))
        colnames(gdsc_rnaseq_norm) <- c(train_ids, test_ids)
        
        sen_res_rnaseq <- cbind.data.frame(sen_res_genes_fltr, gdsc_rnaseq_norm)
        rownames(pheno_data_ca_dr_vs) <- pheno_data_ca_dr_vs[, 1]
        sen_res_pheno <- cbind.data.frame(GROUP = c("", ""), t(pheno_data_ca_dr_vs[c(train_ids, test_ids),]))
        sen_res_rnaseq_pheno <- rbind.data.frame(sen_res_pheno, CLASS = c("", rep("TRAIN", length(train_ids)), 
                                                                          rep("TEST", length(test_ids))),
                                                 sen_res_rnaseq)
        
        # setwd("...")
        write.table(sen_res_rnaseq_pheno, paste0(pheno, "_", ca, "_", dr, "_", vs, "_", dt, "_table.tsv"), quote = F, row.names = T, col.names = F, sep = "\t")
        
        
      }
    }
  }
}

setwd("...")
write.table(SAMPLE_IDs_Class, "STAD 21 DRUGS SAMPLE IDs CLASS.tsv", sep = "\t", row.names = F)


