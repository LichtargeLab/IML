#
# IML figures
#
# Author: Saeid Parvandeh, April, 2024
# -----------------------------------------
rm (list = ls())
library(dplyr)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggbeeswarm)
library(forcats)

#########################################
#### --- boxplot for whole drugs --- ####
setwd("...")
pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
drugs_name <- unique(pheno_data$DRUG_NAME)

incomplete_drugs_name <- c('BMS-536924', 'Crizotinib', 'JW-7-52-1', 'Seliciclib', # error 
                            'Cyclopamine', 'AZ628', 'Sorafenib', 'Tozasertib', 'Saracatinib', # error 
                            'S-Trityl-L-cysteine', 'Bicalutamide', 'Ponatinib',"BX-912", # error
                            "GSK1070916", # cell_11
                            "GSK429286A", "Cabozantinib", "WZ3105", "Olaparib", "KU-55933", "NU7441",
                            "XMD14-99", "Quizartinib", "VNLG/124", "SU11274", "Navitoclax",# incomplete
                            "CGP-082996", "CX-5461", "PI-103", "TPCA-1", "Cytarabine", "Tretinoin",
                            "Tenovin-6", "TWS119", "Torin 2", "Pilaralisib", "GSK1059615", "Voxtalisib",
                            "Brivanib, BMS-540215", "BIBF-1120","AST-1306", "Apitolisib", "LIMK1 inhibitor BMS4",
                            "eEF2K Inhibitor, A-484954", "MetAP2 Inhibitor, A832234", "Venotoclax", "CPI-613",
                            "CAY10566", "Ara-G", "Alisertib", "RO-3306", "EHT-1864", "Cetuximab", "XMD11-85h", 
                            'Rucaparib', 'QL-XII-61', "ICL1100013", "Dyrk1b_0191", "FGFR_0939", "FGFR_3831",
                            "PARP_9482", "PFI-3", "QL-VIII-58")

drugs_name <- unique(drugs_name[!drugs_name %in% incomplete_drugs_name])


######################################
####### ------ One data ------ #######
setwd("...")
data_types <- c('cnv', 'meth', 'mut', 'rnaseq')
iml_cors <- NULL
dr_dt_cors <- data.frame()
for (dr in drugs_name){
  dt_cors <- NULL
  for (dt in data_types){
    pred_res <- get(load(paste0(dr, "_", dt, "_res.RData")))
    iml_cors <- c(iml_cors, pred_res[[1]]^(1/2))
    dt_cors <- c(dt_cors, pred_res[[1]]^(1/2))
  }
  dr_dt_cors <- rbind.data.frame(dr_dt_cors, data.frame(Drug = dr, CNV = dt_cors[1], 
                                                        Methylation = dt_cors[2], 
                                                        Mutation = dt_cors[3], 
                                                        RNASeq = dt_cors[4]))
}


# Save all drugs prediction correlations with four data types.
setwd("...")
write.table(dr_dt_cors, file = "All_Drugs_Correlation_Prediction_PANCAN.tsv", row.names = F, quote = F, sep = "\t")

# Iorio drug cors
setwd("...")
iorio_res <- read.csv("Iorio_results.csv", header = TRUE)
iorio_res_fltr <- iorio_res[which(iorio_res$Cancer.Type == "PANCAN" & iorio_res$Method == "EN"), ]
iorio_cors <-  NULL
drugs_name[drugs_name=='Nutlin-3a (-)'] <- 'Nutlin-3a'
for (dr in drugs_name){
  iorio_res_fltr2 <- iorio_res_fltr[which(grepl(dr, iorio_res_fltr$Drug.name, fixed = TRUE)), 
                                    c("RACS", "gene.exp.", "iCpG", "mutation", "tissue.label", "Pred..power")]
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==1)[c(1, 3, 4, 2)]
  iorio_cors <- c(iorio_cors, iorio_res_fltr2$Pred..power[single_drug_idx])
}
method <- c(rep("Iorio", 1132), rep("IML", 1132))
drug <- c(rep(drugs_name, each=4))
# create a dataset
boxplot_df_single <- data.frame(
  dt=rep(c("CNV", "Methylation", "Mutation", "RNAseq"), 283),
  cor=c(iorio_cors, iml_cors),
  method = fct_rev(method),
  drug = drug
)
boxplot_df_single$dt <- factor(boxplot_df_single$dt, levels = c("Mutation", "CNV", "Methylation", "RNAseq"), ordered = TRUE)
boxplot_df_single$method <- factor(boxplot_df_single$method, levels = c('Iorio', 'IML'), ordered = TRUE)


# Plot
p1 <- boxplot_df_single %>%
  ggplot(aes(x=dt, y=cor,  color= method)) +
  geom_boxplot(outlier.alpha = 0,  
               position = position_dodge(width = 0.78, preserve = "single")) +
  scale_color_manual(name = 'method', values = c('IML' = "red", 'Iorio' = "blue"), labels = c('IML', 'Iorio')) +
  scale_fill_discrete(guide = guide_legend(reverse=FALSE)) +
  geom_point(aes(x = dt, y=cor, group = method, color =factor(method)),
             size = 1.5, shape = 21, alpha=0.8, fill = "gray", position = position_jitterdodge()) +
  stat_compare_means(aes(method = cor), method = "t.test", label = "p.format", hjust = .9, size = 5, show.legend=FALSE) +
  theme_bw()  +
  theme(legend.position = c(0.12, 0.9),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title = element_text(size=15,face="bold")) +
  ggtitle("Individual omics data") +
  labs(x="Omics data type", y= "") +
  coord_flip()


######################################
####### ------ Two data ------ #######
setwd("~/Documents/Drug Response/IML/Data/Phenotype")
pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
drugs_name <- unique(pheno_data$DRUG_NAME)
data_types <- c("rnaseq", "mut", "cnv", "meth")
dt_comb <- combn(data_types, 2)
                             
dr_2dt_cors <- NULL
for (dr in drugs_name){
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/double_data/", dr, "_", dt[1], "_", dt[2], "_ffs_res.RData")
    if (file.exists(path)){
      pred_res <- get(load(path))
      dr_2dt_cors <- rbind(dr_2dt_cors, data.frame(drug=dr, dt=paste0(dt, collapse = '_'), cor=pred_res[[1]]^(1/2)))
    }
  }
}
drugs_name <- names(which(table(dr_2dt_cors$drug)==6))

setwd("~/Documents/Drug Response/IML/Results")
iorio_res <- read.csv("Iorio_results.csv", header = TRUE)
iorio_res_fltr <- iorio_res[which(iorio_res$Cancer.Type == "PANCAN" & iorio_res$Method == "EN"), ]
iorio_cors <-  data.frame()
drugs_name <- drugs_name[!drugs_name %in% c("AICA Ribonucleotide", "Entinostat", "Fedratinib", "GSK319347A", "Ispinesib Mesylate",
                                            "Navitoclax", "Palbociclib", 'PD0325901', "Refametinib", "Saracatinib", "Tanespimycin", 
                                            "Tivozanib", "Tretinoin")]
drugs_name[drugs_name=='Nutlin-3a (-)'] <- 'Nutlin-3a'
for (dr in drugs_name){
  iorio_res_fltr2 <- iorio_res_fltr[which(grepl(dr, iorio_res_fltr$Drug.name, fixed = TRUE)), c("RACS", "gene.exp.", "iCpG", "mutation", "tissue.label", "Pred..power")]
  colnames(iorio_res_fltr2) <- c('cnv', 'rnaseq', 'meth', 'mut', 'tissue.label', 'cor')
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==2)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2])])==2), 'cor']
    iorio_cors <- rbind(iorio_cors, data.frame(drug=dr, dt=paste0(dt, collapse = '_'), cor=rho))
  }
}

# compare two data in iorio and icncv
boxplot_df_double <- rbind(iorio_cors, dr_2dt_cors)
boxplot_df_double$method <- fct_rev(c(rep("Iorio", 528), rep("IML", 640)))
boxplot_df_double$dt <- factor(boxplot_df_double$dt, levels = c("cnv_meth", "mut_cnv", "mut_meth", "rnaseq_mut", "rnaseq_cnv", "rnaseq_meth"), ordered = TRUE)

# Plot
# multiple data - comparing to Iorio
p2 <- boxplot_df_double %>%
  ggplot(aes(x=dt, y=cor,  color= method)) +
  geom_boxplot(outlier.alpha = 0,  
               position = position_dodge(width = 0.78, preserve = "single")) + 
  scale_color_manual(name = 'method', values = c('IML' = "red", 'Iorio' = "blue"), labels = c('IML', 'Iorio')) +
  scale_fill_discrete(guide = guide_legend(reverse=FALSE)) +
  expand_limits(y = c(0.0, 0.85)) +
  geom_point(aes(x = dt, y=cor, group = method, color =factor(method)),
             size = 1.5, shape = 21, alpha=0.8, fill = "gray", position = position_jitterdodge()) +
  stat_compare_means(aes(method = cor), method = "t.test", label = "p.format", hjust = .8, size=5, show.legend=FALSE) +#, angle=90
  theme_bw()  +
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title = element_text(size=15,face="bold")) +
  ggtitle("Two omics data combinations") +
  labs(x="", y= "") +
  coord_flip()


################################################
####### ------ Tree and four data ------ #######
setwd("~/Documents/Drug Response/IML/Data/Phenotype")
pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
drugs_name <- unique(pheno_data$DRUG_NAME)
data_types <- c("rnaseq", "mut", "cnv", "meth")
dt_comb <- combn(data_types, 3)

dr_3_all_dt_cors <- NULL
for (dr in drugs_name){
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/triple_data/", dr, "_", dt[1], "_", dt[2], "_", dt[3], "_ffs_res.RData")
    if (file.exists(path)){
      pred_res <- get(load(path))
      dr_3_all_dt_cors <- rbind(dr_3_all_dt_cors, data.frame(drug=dr, dt=paste0(dt, collapse = '_'), cor=pred_res[[1]]^(1/2)))
    }
  }
  path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/all_data/", dr, "_all_ffs_res.RData")
  if (file.exists(path)){
    pred_res <- get(load(path))
    dr_3_all_dt_cors <- rbind(dr_3_all_dt_cors, data.frame(drug=dr, dt='all', cor=pred_res[[1]]^(1/2)))
    
  }
}
drugs_name <- names(which(table(dr_3_all_dt_cors$drug)==5))

setwd("~/Documents/Drug Response/IML/Results")
iorio_res <- read.csv("Iorio_results.csv", header = TRUE)
iorio_res_fltr <- iorio_res[which(iorio_res$Cancer.Type == "PANCAN" & iorio_res$Method == "EN"), ]
iorio_cors <-  data.frame()
drugs_name <- drugs_name[!drugs_name %in% c("Refametinib", "Navitoclax", "AICA Ribonucleotide", "IC-87114",
                                            "NSC319726")]
drugs_name[drugs_name=='Nutlin-3a (-)'] <- 'Nutlin-3a'
drugs_name[drugs_name == "PD0325901"] <- "PD-0325901"
for (dr in drugs_name){
  iorio_res_fltr2 <- iorio_res_fltr[which(grepl(dr, iorio_res_fltr$Drug.name, fixed = TRUE)), c("RACS", "gene.exp.", "iCpG", "mutation", "tissue.label", "Pred..power")]
  colnames(iorio_res_fltr2) <- c('cnv', 'rnaseq', 'meth', 'mut', 'tissue.label', 'cor')
  # three datasets
  dt_comb <- combn(data_types, 3)
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==3)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2], dt[3])])==3), 'cor']
    iorio_cors <- rbind(iorio_cors, data.frame(drug=dr, dt=paste0(dt, collapse = '_'), cor=rho))
  }
  # four datasets
  dt_comb <- combn(data_types, 4)
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==4)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2], dt[3], dt[4])])==4), 'cor']
    iorio_cors <- rbind(iorio_cors, data.frame(drug=dr, dt='all', cor=rho))
  }
}

# compare three data in iorio and icncv
boxplot_df_triple <- rbind(iorio_cors, dr_3_all_dt_cors)
boxplot_df_triple$method <- fct_rev(c(rep("Iorio", 170), rep("IML", 190)))
boxplot_df_triple$dt <- factor(boxplot_df_triple$dt, levels = c("mut_cnv_meth", "rnaseq_mut_cnv", "rnaseq_cnv_meth", "rnaseq_mut_meth", "all"), ordered = TRUE)


# Plot
# multiple data - comparing to Iorio
p3 <- boxplot_df_triple %>%
  ggplot(aes(x=dt, y=cor,  color= method)) +
  geom_boxplot(outlier.alpha = 0,  
               position = position_dodge(width = 0.78, preserve = "single")) + 
  scale_color_manual(name = 'method', values = c('IML' = "red", 'Iorio' = "blue"), labels = c('IML', 'Iorio')) +
  scale_fill_discrete(guide = guide_legend(reverse=FALSE)) +
  expand_limits(y = c(0.0, 0.8)) +
  geom_point(aes(x = dt, y=cor, group = method, color =factor(method)),
             size = 1.5, shape = 21, alpha=0.8, fill = "gray", position = position_jitterdodge()) +
  stat_compare_means(aes(method = cor), method = "t.test", label = "p.format", hjust = 3.0, size=5, show.legend=FALSE) +#, angle=90
  theme_bw()  +
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title = element_text(size=15,face="bold")) +
  ggtitle("Three and four omics data\nTop 10 IML drugs") +
  labs(x="Omics data type", y= "Pearson correlation") +
  coord_flip()


# top performed Iorio's drug
setwd("~/Documents/Drug Response/IML/Results")
iorio_res <- read.csv("Iorio_results.csv", header = TRUE)
iorio_res_fltr <- iorio_res[which(iorio_res$Cancer.Type == "PANCAN" & iorio_res$Method == "EN"), ]
single_data_idx <- which(rowSums(iorio_res_fltr[, 4:8])==1)
iorio_res_fltr2 <- iorio_res_fltr[single_data_idx, ]
iorio_drugs <- NULL
for (dt in c("RACS", "gene.exp.", "iCpG", "mutation")){
  iorio_res_trim <- iorio_res_fltr2[which(iorio_res_fltr2[, dt]==1), ]
  iorio_res_sort <- iorio_res_trim[order(iorio_res_trim$Pred..power, decreasing = T), 'Drug.name']
  iorio_drugs <- c(iorio_drugs, iorio_res_sort[1:10])
}
iorio_drugs <- unique(iorio_drugs)

iorio_cors <-  data.frame()
for (dr in iorio_drugs){
  iorio_res_fltr2 <- iorio_res_fltr[which(grepl(dr, iorio_res_fltr$Drug.name, fixed = TRUE)), c("RACS", "gene.exp.", "iCpG", "mutation", "tissue.label", "Pred..power")]
  colnames(iorio_res_fltr2) <- c('cnv', 'rnaseq', 'meth', 'mut', 'tissue.label', 'cor')
  # three datasets
  dt_comb <- combn(data_types, 3)
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==3)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2], dt[3])])==3), 'cor']
    iorio_cors <- rbind(iorio_cors, data.frame(drug=dr, dt=paste0(dt, collapse = '_'), cor=rho))
  }
  # four datasets
  dt_comb <- combn(data_types, 4)
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==4)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2], dt[3], dt[4])])==4), 'cor']
    iorio_cors <- rbind(iorio_cors, data.frame(drug=dr, dt='all', cor=rho))
  }
}

setwd("~/Documents/Drug Response/IML/Data/Phenotype")
data_types <- c("rnaseq", "mut", "cnv", "meth")
dt_comb <- combn(data_types, 3)

# iorio_drugs = iorio_drugs[!iorio_drugs%in% c("RDEA119 (rescreen)")]
iorio_drugs[which(iorio_drugs == "Afatinib (rescreen)")] <- "Afatinib"
iorio_drugs[which(iorio_drugs == "Nutlin-3a")] <- "Nutlin-3a (-)"
iorio_drugs[which(iorio_drugs == "PD-0325901")] <- "PD0325901"

dr_3_all_dt_cors <- NULL
for (dr in iorio_drugs){
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/triple_data/", dr, "_", dt[1], "_", dt[2], "_", dt[3], "_ffs_res.RData")
    if (file.exists(path)){
      pred_res <- get(load(path))
      dr_3_all_dt_cors <- rbind(dr_3_all_dt_cors, data.frame(drug=dr, dt=paste0(dt, collapse = '_'), cor=pred_res[[1]]^(1/2)))
    }
  }
  path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/all_data/", dr, "_all_ffs_res.RData")
  if (file.exists(path)){
    pred_res <- get(load(path))
    dr_3_all_dt_cors <- rbind(dr_3_all_dt_cors, data.frame(drug=dr, dt='all', cor=pred_res[[1]]^(1/2)))
    
  }
}

boxplot_df_triple <- rbind(iorio_cors, dr_3_all_dt_cors)
boxplot_df_triple$method <- fct_rev(c(rep("Iorio", 130), rep("IML", 90)))
boxplot_df_triple$dt <- factor(boxplot_df_triple$dt, levels = c("mut_cnv_meth", "rnaseq_mut_cnv", "rnaseq_cnv_meth", "rnaseq_mut_meth", "all"), ordered = TRUE)


# Plot
# multiple data - comparing to Iorio
p4 <- boxplot_df_triple %>%
  ggplot(aes(x=dt, y=cor,  color= method)) +
  geom_boxplot(outlier.alpha = 0,  
               position = position_dodge(width = 0.78, preserve = "single")) + 
  scale_color_manual(name = 'method', values = c('IML' = "red", 'Iorio' = "blue"), labels = c('IML', 'Iorio')) +
  scale_fill_discrete(guide = guide_legend(reverse=FALSE)) +
  expand_limits(y = c(0.0, 0.8)) +
  geom_point(aes(x = dt, y=cor, group = method, color =factor(method)),
             size = 1.5, shape = 21, alpha=0.8, fill = "gray", position = position_jitterdodge()) +
  stat_compare_means(aes(method = cor), method = "t.test", label = "p.format", hjust = 3.0, size=5, show.legend=FALSE) +#, angle=90
  theme_bw()  +
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title = element_text(size=15,face="bold")) +
  ggtitle("Three and four omics data\nTop 10 Iorio's drugs") +
  labs(x="", y= "Pearson correlation") +
  coord_flip()

gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

#########################################################
########---- Multiple data - self comparison ---- #######
setwd("~/Documents/Drug Response/IML/Data/Phenotype")
pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
drugs_name <- unique(pheno_data$DRUG_NAME)
data_types <- c("rnaseq", "mut", "cnv", "meth")
dr_dt_cors <- data.frame()
for (dr in drugs_name){
  for (dt in data_types){
    path <- paste0('~/Documents/Drug Response/IML/Results/PANCAN/single_data/', dr, '_', dt, '_res.RData')
    if (file.exists(path)){
      load(path)
      dr_dt_cors <- rbind(dr_dt_cors, data.frame(drug = dr, data_types=dt, cor=cncv_res[[1]]^(1/2)))
    }
  }
  dt_comb <- combn(data_types, 2)
  for (m in 1:ncol(dt_comb)) {
    path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/double_data/", dr, "_", dt_comb[1, m], "_", dt_comb[2, m], "_ffs_res.RData")
    if (file.exists(path)){
      load(path)
      dr_dt_cors <- rbind(dr_dt_cors, data.frame(drug = dr, data_types=paste0(dt_comb[, m], collapse = "_"), cor = icncv_res[[1]]^(1/2)))
    }
  }
  dt_comb <- combn(data_types, 3)
  for (m in 1:ncol(dt_comb)) {
    path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/triple_data/", dr, "_", dt_comb[1, m], "_", dt_comb[2, m], "_", dt_comb[3, m], "_ffs_res.RData")
    if (file.exists(path)){
      load(path)
      dr_dt_cors <- rbind(dr_dt_cors, data.frame(drug = dr, data_types=paste0(dt_comb[, m], collapse = "_"), cor = icncv_res[[1]]^(1/2)))
    }
  }
  path <- paste0("~/Documents/Drug Response/IML/Results/PANCAN/all_data/", dr, "_all_ffs_res.RData")
  if (file.exists(path)){
    load(path)
    dr_dt_cors <- rbind(dr_dt_cors, data.frame(drug = dr, data_types="all", cor = icncv_res[[1]]^(1/2)))
  }
}
# Find drugs that have all 15 data combinations
unique_drugs <- names(which(table(dr_dt_cors$drug)==15))
drugs_name <- unique_drugs

# plot
dt_res <- dr_dt_cors[dr_dt_cors$drug %in% drugs_name, ]
dt_res$data_types <- factor(dt_res$data_types, levels = c("mut", "cnv", "cnv_meth", "mut_meth", "mut_cnv", "mut_cnv_meth", 
                                                          "meth", "rnaseq", "rnaseq_cnv", "rnaseq_meth", "rnaseq_mut", 
                                                          "rnaseq_mut_cnv", "rnaseq_cnv_meth", "rnaseq_mut_meth", "all"), ordered = TRUE)

my_comparisons <- list( c("rnaseq_mut",  "rnaseq"), c( "rnaseq_meth", "rnaseq"), c("mut_meth", "mut"), 
                        c("all",  "rnaseq"), c("mut_cnv", "mut"),c("rnaseq_meth", "meth"), c("cnv_meth", "cnv"))
dt_res %>%
  ggplot(aes(x=data_types, y=cor, col=data_types))+#, fill=factor(data_types))) +#, levels = c("Iorio", "icncv")))) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  expand_limits(y = c(0.0, 0.75)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes=factor(data_types), size=1.5, alpha=0.5) +
  scale_color_manual(values=c("#DF73D4", "#FF6700", "#046307", "#6495ED", "#0000BB", "#008000", "#808080", "#FF00FF", 
                              "#00FFFF", "#FFA500", "#0000FF", "#9999FF", "#99FFCC", "#FF3300", "#663300")) +
  stat_compare_means(method = "t.test", method.args = list(alternative="greater"), comparisons = my_comparisons, size = 5)+ # Add pairwise comparisons p-value
  guides(colour = guide_legend(reverse=T)) + 
  theme_bw()  +
  theme(
    legend.position = 'none',
    axis.text=element_text(size=15),
    axis.title=element_text(size=15),
    plot.title = element_text(size=15,face="bold")) +
  labs(x="Integrated and Individual Omics Data", y= "Pearson correlation") + #  of predicted vs. observed IC50
  coord_flip()

##################################################
### --- cancer specific comparing to Iorio --- ###
setwd("~/Documents/Drug Response/IML/Data/Phenotype")
pheno_data <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", colClasses = "character")
all_individual_tumors <- unique(pheno_data$TCGA_DESC[pheno_data$TCGA_DESC != ""])
cancers <- NULL
for (ca in all_individual_tumors){
  pheno_data_ca <- pheno_data[pheno_data$TCGA_DESC == ca, ]
  if (length(unique(pheno_data_ca$SANGER_MODEL_ID)) > 30){
    cancers <- c(cancers, ca)
  }
}
cancers <- cancers[which(cancers!="UNCLASSIFIED")]

cancers <- c("COREAD", "SKCM", "HNSC", "OV")#, "BRCA", "GBM", "LUAD", "SCLC"
ca <- cancers[1]
pheno_data_ca <- pheno_data[pheno_data$TCGA_DESC %in% ca, ]

data_types <- c("rnaseq", "mut", "cnv", "meth")
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
IML_drugs_name <- NULL
for (dr in drugs_name){
  for (dt in data_types){
    path <- paste0("~/Documents/Drug Response/IML/Results/Subtypes/", ca, "/", ca, "_", dr, "_", dt, "_res.RData")
    if (file.exists(path)){
      load(path)
      IML_drugs_name <- c(IML_drugs_name, dr)
    }
  }
}

IML_drugs_name <- unique(IML_drugs_name)

# COREAD
IML_drugs_name <- IML_drugs_name[!IML_drugs_name %in% c("AZ628", "Pyrimethamine", "BMS-509744", "BI-2536", "GW843682X")]


# SKCM
IML_drugs_name <- IML_drugs_name[!IML_drugs_name %in% c("AZ628", "MG-132", "Bortezomib")]

# HNSC
IML_drugs_name <- IML_drugs_name[!IML_drugs_name %in% c("TGX221", "Erlotinib", "Lapatinib", "WZ-1-84", "Saracatinib")]


# OV
IML_drugs_name <- IML_drugs_name[!IML_drugs_name %in% c("Bortezomib", "GSK319347A", "NSC-87877", "Pyrimethamine",
                                            "A-770041", "Erlotinib", "PAC-1")]


setwd(paste0("~/Documents/Drug Response/IML/Results/Subtypes/", ca))
IML_res <- NULL
for (dr in IML_drugs_name) {
  for (dt in data_types){
    pred_res <- get(load(paste0(ca, "_", dr, "_", dt, "_res.RData")))
    # pred_res <- get(load(paste0("GDSC1_",ca, "_", dr, "_", dt, "_res.RData")))
    pred_cor <- cor(pred_res[[2]][, 1], pred_res[[2]][, 2])^2
    IML_res <- rbind(IML_res, data.frame(drug=dr, cor = pred_cor, data_types=dt))#^(1/2)
  }
  pred_res <- get(load(paste0(ca, "_", dr, "_all_ffs_res.RData")))
  pred_cor <- cor(pred_res[[2]][, 1], pred_res[[2]][, 2])^2
  IML_res <- rbind(IML_res, data.frame(drug=dr, cor = pred_cor, data_types="all"))#^(1/2)
  dt_comb <- combn(data_types, 2)
  for (m in 1:ncol(dt_comb)) {
    pred_res <- get(load(paste0(ca, "_", dr, "_", dt_comb[1, m], "_", dt_comb[2, m], "_ffs_res.RData")))
    pred_cor <- cor(pred_res[[2]][, 1], pred_res[[2]][, 2])^2
    IML_res <- rbind(IML_res, data.frame(drug=dr, cor = pred_cor, data_types=paste0(dt_comb[, m], collapse = "_")))#^(1/2)
  }
  dt_comb <- combn(data_types, 3)
  for (m in 1:ncol(dt_comb)) {
    pred_res <- get(load(paste0(ca, "_", dr, "_", dt_comb[1, m], "_", dt_comb[2, m], "_", dt_comb[3, m], "_ffs_res.RData")))
    pred_cor <- cor(pred_res[[2]][, 1], pred_res[[2]][, 2])^2
    IML_res <- rbind(IML_res, data.frame(drug=dr, cor = pred_cor, data_types=paste0(dt_comb[, m], collapse = "_")))#^(1/2)
  }
}

# Iorio
setwd("~/Documents/Drug Response/IML/Results")
iorio_res <- read.csv("Iorio_results.csv", header = TRUE)
iorio_res_fltr <- iorio_res[which(iorio_res$Cancer.Type == ca & iorio_res$Method == "EN"), ]
iorio_res_fltr <- iorio_res_fltr[-which(iorio_res_fltr$Drug.name == "Afatinib (rescreen)"), ]

# top sensitive drugs based on average drug response across cell lines
Iorio_drugs_name <- names(sort(avg_dr))[1:50]
# COREAD
Iorio_drugs_name <- Iorio_drugs_name[!Iorio_drugs_name %in% c("CAY10566", "Selumetinib", "Linsitinib", "Refametinib", "Pemetrexed",
                                                              "Enzastaurin", "RAF_9304", "AZD8931", "Pelitinib", "GSK650394",
                                                              "AZD8835", "Tanespimycin", "EHT-1864", "Avagacestat")]

# SKCM
Iorio_drugs_name <- Iorio_drugs_name[!Iorio_drugs_name %in% c("PLX-4720", "LIMK1 inhibitor BMS4", "Selumetinib", "RAF_9304",
                                            "Refametinib", "GSK1059615", "Tanespimycin", "Voxtalisib",
                                            "Pilaralisib", "IAP_7638", "AST-1306", "Selisistat", "LDN-193189",
                                            "WYE-125132", "CAP-232, TT-232, TLN-232", "FGFR_0939", "IAP_5620",
                                            "Motesanib")]
# HNSC
Iorio_drugs_name <- Iorio_drugs_name[!Iorio_drugs_name %in% c("TGX221", "Erlotinib", "Lapatinib", "AZD8931", "WZ-1-84",
                                            "IAP_7638", "IAP_5620", "Bleomycin (10 uM)", "Saracatinib",
                                            "Wee1 Inhibitor", "AZD5582", "AST-1306", "GSK1059615",
                                            "MG-132", "BPTES", "Tanespimycin", "Kobe2602", "MetAP2 Inhibitor, A832234",
                                            "XAV939", "WHI-P97", "AZ628", "Tozasertib", "CPI-613", "PI3Ka_4409",
                                            "EphB4_9721", "Veliparib", "CAP-232, TT-232, TLN-232", "GW441756",
                                            "AZD3514", "Paclitaxel", "PFI-3", "TANK_1366", "C-75", "Bortezomib",
                                            "Cyclopamine")]

# GBM
drugs_name <- drugs_name[!drugs_name %in% c("AZD6094", "SB216763", "JNJ38877605", "LCL161", "Tanespimycin",
                                            "Rucaparib", "Saracatinib", "Tivozanib", "Cabozantinib")]

# OV
Iorio_drugs_name <- Iorio_drugs_name[!Iorio_drugs_name %in% c("IAP_5620", "TANK_1366", "XAV939", "Bortezomib", "GSK319347A",
                                            "LCL161", "Enzastaurin", "AZD5582", "EphB4_9721", "IAP_7638",
                                            "GW441756", "Pyrimethamine", "A-770041", "BPTES", "Erlotinib",
                                            "AZD8931", "Tanespimycin", "KIN001-042", "WHI-P97", "MG-132", 
                                            "Dactolisib", "JW-7-52-1", "Pictilisib", "AZD5438", "AZD1208",
                                            "AZD8835", "PI3Ka_4409", "PARP_9495", "JAK3_7406", "Wee1 Inhibitor")]

Iorio_drugs_name[Iorio_drugs_name=='Nutlin-3a (-)'] <- 'Nutlin-3a'
Iorio_drugs_name[Iorio_drugs_name=='PD0325901'] <- 'PD-0325901'
Iorio_drugs_name[Iorio_drugs_name=='PLX-4720'] <- 'PLX4720'
Iorio_res <- NULL
for (dr in Iorio_drugs_name){
  iorio_res_fltr2 <- iorio_res_fltr[which(grepl(dr, iorio_res_fltr$Drug.name, fixed = TRUE)), c("RACS", "gene.exp.", "iCpG", "mutation", "tissue.label", "Pred..power")]
  colnames(iorio_res_fltr2) <- c('cnv', 'rnaseq', 'meth', 'mut', 'tissue.label', 'cor')
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==1)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  for (dt in data_types){
    rho <- iorio_res_fltr3[which(iorio_res_fltr3[, dt[1]]==1), 'cor']
    Iorio_res <- rbind(Iorio_res, data.frame(cor=(rho)^2, data_types=dt, drug=dr))
  }
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==2)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  dt_comb <- combn(data_types, 2)
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2])])==2), 'cor']
    Iorio_res <- rbind(Iorio_res, data.frame(cor=(rho)^2, data_types=paste0(dt, collapse = '_'), drug=dr))
  }
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==3)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  dt_comb <- combn(data_types, 3)
  for (m in 1:ncol(dt_comb)){
    dt <- dt_comb[, m]
    rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(dt[1], dt[2], dt[3])])==3), 'cor']
    Iorio_res <- rbind(Iorio_res, data.frame(cor=(rho)^2, data_types=paste0(dt, collapse = '_'), drug=dr))
  }
  single_drug_idx <- which(rowSums(iorio_res_fltr2[, 1:5])==4)#[c(1, 3, 4, 2)]
  iorio_res_fltr3 <- iorio_res_fltr2[single_drug_idx, ]
  rho <- iorio_res_fltr3[which(rowSums(iorio_res_fltr3[, c(data_types[1], data_types[2], data_types[3], data_types[4])])==4), 'cor']
  Iorio_res <- rbind(Iorio_res, data.frame(cor=(rho)^2, data_types=paste0("all", collapse = '_'), drug=dr))
}

# compare two data in iorio and icncv
boxplot_df <- rbind(Iorio_res, IML_res)
boxplot_df$method <- fct_rev(c(rep("Iorio", nrow(Iorio_res)), rep("IML", nrow(IML_res))))
boxplot_df$cor=c(Iorio_res$cor, IML_res$cor)
boxplot_df$data_types <- factor(boxplot_df$data_types, levels = c("mut","cnv","meth","rnaseq",    
                                                                  "cnv_meth", "mut_cnv", "mut_meth", "rnaseq_mut", "rnaseq_cnv", "rnaseq_meth",
                                                                  "rnaseq_mut_cnv", "rnaseq_mut_meth", "rnaseq_cnv_meth", "mut_cnv_meth", "all"), ordered = TRUE)
my_comparisons <- list( c("rnaseq_mut",  "rnaseq"), c( "rnaseq_meth", "rnaseq"), c("mut_meth", "mut"), 
                        c("mut_cnv", "mut"),c("rnaseq_meth", "meth"), c("cnv_meth", "cnv"))

p1 <- boxplot_df %>%
  ggplot(aes(x=data_types, y=cor,  color= method)) +
  geom_boxplot(outlier.alpha = 0,  
               position = position_dodge(width = 0.78, preserve = "single")) + 
  scale_color_manual(name = 'method', values = c('IML' = "red", 'Iorio' = "blue"), labels = c('IML', 'Iorio')) +
  scale_fill_discrete(guide = guide_legend(reverse=FALSE)) +
  expand_limits(y = c(0.0, 0.6)) +
  geom_point(aes(x = data_types, y=cor, group = method, color =factor(method)),
             size = 1.5, shape = 21, alpha=0.8, fill = "gray", position = position_jitterdodge()) +
  stat_compare_means(aes(method = cor), method = "t.test", paired = F, label = "p.format", hjust = 0.8, size = 5, show.legend=FALSE) +#, angle=90
  theme_bw()  +
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title = element_text(size=15,face="bold")) +
  ggtitle(paste0(ca)) + #, "-Omics data combinations"
  labs(x="", y= "") + # Omics data type # R-squared
  coord_flip()

gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)


