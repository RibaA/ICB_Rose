## load libraries
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]
#annot_dir <- args[3]

input_dir <- "data/input"
output_dir <- "data/output"
annot_dir <- "data/annot"

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" , header=TRUE )
rownames(clin) <- clin$patientid
cols <- c('patientid', 'gender.ch1', 'age.ch1', 'histology.ch1', 'n.stage.at.diagnosis.ch1',
          'io.response.ch1',
          'alive.ch1', 'overall.survival.ch1', 'progressed.ch1', 'pfs.ch1')

clin <- clin[, c(cols, colnames(clin)[!colnames(clin) %in% cols])]
colnames(clin)[colnames(clin) %in% cols] <- c('patient', 'sex', 'age', 'histo', "stage", "recist",
                                              'os', 't.os', 'pfs', 't.pfs')

new_cols <- c( "primary" , "response.other.info" , "response" , "drug_type" , "dna" , "rna", "rna.info")
clin[new_cols] <- NA

clin$drug_type <- 'PD-1/PD-L1'
clin$rna <- 'rnaseq'
clin$rna.info <- 'tpm'
clin$primary <- ifelse(clin$primary.tumor.location.ch1 != 'Bladder', 'Urothelial', 'Bladder')
clin$t.os <- clin$t.os / 30.44
clin$t.pfs <- clin$t.pfs /30.44

clin$response = Get_Response( data=clin )

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Rose', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, treatmentid=clin$drug_type, .after='tissueid')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

