#
# Preprocessing multi-omics data from GDSC & CCLE
#
#################################

# Reading drugs phenotype data
#
# setwd("...")
drug_pheno <- read.csv("GDSC1_fitted_dose_response_25Feb20.csv", header = TRUE)

# RNA-Seq
# data downloaded from https://cellmodelpassports.sanger.ac.uk/downloads
# Sanger - GDSC
rnaseq_data <- read.csv("rnaseq_fpkm_20210329.csv", header = T)
rnaseq_data$gene <- make.unique(rnaseq_data$gene)
rownames(rnaseq_data) <- rnaseq_data$gene
rnaseq_data <- rnaseq_data[-c(1,2,3), -c(1, 2)]
rnaseq_data <- t(rnaseq_data)
# convert character to numeric in matrix
mode(rnaseq_data) <- "numeric"
save(rnaseq_data, file = "rnaseq_data.RData")
# Broad - CCLE
# data downloaded from https://data.broadinstitute.org/ccle/
gene_expr <- read.table("CCLE_RNAseq_genes_rpkm_20180929.gct", skip = 2, row.names = 1, header = T)
# colnames(gene_expr) <- gsub("\\..*","",colnames(gene_expr))
colnames(gene_expr) <- gsub("X","",colnames(gene_expr))
ccle_rnaseq_data <- gene_expr[, c('Description', intersect(colnames(gene_expr), unique(ccle_pheno_data$CCLE_Name)))]
save(ccle_rnaseq_data, file = "ccle_rnaseq_data.RData")
write.table(ccle_rnaseq_data, file= "CCLE_RNAseq_genes_rpkm_20180929.tsv", row.names = F, quote = F, sep = "\t")

# Methylation
# Pre-Processed beta values for all CpG islands across all the cell-lines
# Data downloaded from https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html
meth_data <- read.table("F2_METH_CELL_DATA.txt", header = TRUE)
colnames(meth_data) <- gsub("X","",colnames(meth_data))
sample_id <- read.table("sample_id.txt", header = F)
sample_id <- unique(sample_id)
sample_id_file <- read.csv("sample.csv", header = TRUE)
sample_id$cellLine_name <- sample_id_file$Title[match(sample_id$V1, sample_id_file$Accession)]
sample_id$cellLine_name <- gsub("CellLine_","",sample_id$cellLine_name)
sample_id$SANGER_MODEL_ID <- drug_pheno$SANGER_MODEL_ID[match(sample_id$cellLine_name, drug_pheno$CELL_LINE_NAME)]

colnames(meth_data) <- sample_id$SANGER_MODEL_ID[match(colnames(meth_data), sample_id$V2)]

meth_data <- t(meth_data)
meth_data <- meth_data[!is.na(rownames(meth_data)), ]
save(meth_data, file="meth_data.RData")

# GDSC Cell Lines Whole Exome Sequencing
# Caveman and Pindel processed and annotated VCFs for whole exome sequencing of cell lines.
# All cavemans merged into one VCF and annotated with EA (VEP tool and hg19 reference genome) 
# data downloaded from https://cellmodelpassports.sanger.ac.uk/downloads
# dots removed before reading using bcftools
vcf_file <- read.delim("normalized.merge.info.all.samples.SMID.rm.EA.fltr_vep_dbNSFP4.2a.nonCanonical.vcf",
                       fill = T, header = T, colClasses = "character", sep = "\t", skip = 152)
variants_info <- data.frame(matrix(NA, nrow = nrow(vcf_file), ncol = 8))

for (r in 1:nrow(vcf_file)){
  info <- unlist(strsplit(vcf_file[r, ]$INFO, ';'))
  # print(info)
  info_ENSP <- grep("ENSP=", info, value = TRUE)
  info_ENSP_trim <- gsub('.*=', '', info_ENSP)
  info_ens_pid <- grep("Ensembl_proteinid", info, value = TRUE)
  info_ens_pid_trim <- gsub('.*=', '', info_ens_pid)
  ens_pid <- unlist(strsplit(info_ens_pid_trim, ','))
  canonical_trans_idx <- which(info_ENSP_trim==ens_pid)
  if (!is.element(info_ENSP_trim, ens_pid)){
    canonical_trans_idx <- 1
  }
  info_cons <- grep("Consequence", info, value = TRUE)
  info_cons_trim <- gsub('.*=', '', info_cons)
  cons_vec <- unlist(strsplit(info_cons_trim, ','))
  vairant_class <- ifelse(length(cons_vec)==length(ens_pid), cons_vec[canonical_trans_idx], cons_vec)
  if (is.element(vairant_class,c('frameshift_variant',
                                 'missense_variant',
                                 'inframe_insertion',
                                 'inframe_deletion',
                                 'splice_acceptor_variant',
                                 'splice_donor_variant',
                                 'start_lost',
                                 'stop_gained',
                                 'stop_lost',
                                 'missense_variant&splice_region_variant&NMD_transcript_variant',
                                 'missense_variant&splice_region_variant',
                                 'splice_acceptor_variant&NMD_transcript_variant',
                                 'stop_lost&NMD_transcript_variant'))){
    info_gene <- grep("SYMBOL", info, value = TRUE)
    info_gene_trim <- gsub('.*=', '', info_gene)
    gene_vec <- unlist(strsplit(info_gene_trim, '_'))
    gene <- ifelse(length(gene_vec)==length(ens_pid), gene_vec[canonical_trans_idx], gene_vec)
    info_EA <- grep("EA", info, value = TRUE)
    info_EA_trim <- gsub('.*=', '', info_EA)
    iso_ea <- unlist(strsplit(info_EA_trim, ','))
    ea <- ifelse(length(iso_ea)==length(ens_pid), iso_ea[canonical_trans_idx], iso_ea)
    ea <- ifelse(vairant_class%in%c('frameshift_variant',
                                    'splice_acceptor_variant',
                                    'splice_donor_variant',
                                    'start_lost', 'stop_gained', 'stop_lost',
                                    'splice_acceptor_variant&NMD_transcript_variant',
                                    'stop_lost&NMD_transcript_variant'), "100", ea)
    variants_info[r, ] <- c(vcf_file[r, c(1, 2, 4, 5)], info_ENSP_trim, gene, vairant_class, ea)
  }
}

colnames(variants_info) <- c("CHROM", "POS", "REF", "ALT", "ENSP", "GENE", "Variant_Classification", "EA")


write.csv(variants_info, "GDSC_maf_EA.csv", row.names = F, quote = F)

genotype <- sapply(vcf_file[, 10:ncol(vcf_file)], substring, 0, 3)


vcf_file_ <- cbind(variants_info, genotype)
vcf_file_fltr <- vcf_file_[!is.na(vcf_file_$EA), ]

# Rcpp
vcf_file_fltr$EA <- as.numeric(vcf_file_fltr$EA)

all_genes <- unique(vcf_file_fltr$GENE)
all_samples <- colnames(genotype)

Rcpp::sourceCpp("VCF2DM.cpp")

GDSC_Mutation_Data <- VCFtoDM(vcf_file_fltr, all_genes, all_samples)

rownames(GDSC_Mutation_Data) <- all_samples
colnames(GDSC_Mutation_Data) <- all_genes

write.table(GDSC_Mutation_Data, "GDSC_Mutation_Data.tsv", sep = "\t")
rownames(mut_data) <- gsub(" ","",rownames(mut_data))
save(mut_data, file = "mut_data.RData")

# CNV
# data downloaded from (https://cellmodelpassports.sanger.ac.uk/downloads)
cnv_data <- read.csv("cnv_gistic_20191101.csv", header = T, row.names = 2)
cnv_data <- t(cnv_data[-c(1,2), -1])
# convert character to numeric in matrix
mode(cnv_data) <- "numeric"
save(cnv_data, file = "cnv_data.RData")
