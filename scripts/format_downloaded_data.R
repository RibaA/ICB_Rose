## load libraries

library(GEOquery)
library(data.table)
library(readxl)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
#work_dir <- args[1]
work_dir <- "data"

# CLIN.txt
gse <- getGEO("GSE176307", GSEMatrix = TRUE, AnnotGPL = TRUE)
dat <- gse[[1]]

# Extract clinical/phenotype (metadata / sample info)
clin <- pData(dat)
clin <- as.data.frame(clin)
clin$patientid <- substr(clin$title, 16, nchar(clin$title))
clin <- clin[order(clin$patientid), ]

# EXP_TPM.tsv 
expr <- read_tsv( file = file.path(work_dir, "GSE176307_salmon_tpm_gene.matrix.tsv") ) 
expr <- as.data.frame(expr)
rownames(expr) <- expr[, 1]
expr <- expr[, -1]
expr <- expr[order(rownames(expr)), ]

# samples with both clinical and expression data
sample_id <- read.csv(file.path(work_dir, 'GSE176307_BACI_Omniseq_Sample_Name_Key_submitted_GEO_v2.csv'))
sample_id <- sample_id[order(sample_id$'Omniseq_RS_ID..RNAseq.'), ]
int <- intersect(sample_id$'Omniseq_RS_ID..RNAseq.', colnames(expr))
expr <- expr[, colnames(expr) %in% int]
expr <- expr[, order(colnames(expr))]
sample_id <- sample_id[sample_id$'Omniseq_RS_ID..RNAseq.' %in% int, ]
colnames(expr) <- sample_id$Sample.ID
expr <- expr[, order(colnames(expr))]

int <- intersect(clin$patientid, colnames(expr))
clin <- clin[clin$patientid %in% int, ]
rownames(clin) <- clin$patientid

write.table(clin, file=file.path(work_dir, 'CLIN.txt'), sep = "\t" , quote = FALSE , row.names = FALSE)

expr <- expr[, colnames(expr) %in% int]
write.table(expr, file=file.path(work_dir, 'EXP_TPM.tsv'), sep = "\t" , quote = FALSE , row.names = TRUE, col.names=TRUE)
