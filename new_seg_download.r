rm(list = ls())
library("TCGAbiolinks")
library(SummarizedExperiment)
setwd("desktop")
makeNset_dir <- function(dir_name)
{if(dir.exists(dir_name)  == FALSE)
{dir.create(dir_name)}
  setwd(dir_name)  
}

makeNset_dir("UCEC_CNV")

#geting CNV data
query.cnv.UCEC <- GDCquery(project = "TCGA-UCEC",
                           data.category = "Copy Number Variation",
                           platform="Affymetrix SNP 6.0",
                           data.type = "Masked Copy Number Segment",
                           sample.type = "Primary Tumor",
                           #barcode = c("TCGA-C8-A26W-01A", "TCGA-D8-A1JA-01A", "TCGA-D8-A27N-01A")
                           )
GDCdownload(query.cnv.UCEC, files.per.chunk = 150)

query.cnv.UCEC <- GDCprepare(query = query.cnv.UCEC,
                             save = TRUE,
                             save.filename = "query.cnv.UCEC",
                             summarizedExperiment = TRUE)

save(query.cnv.UCEC,file="query.cnv.UCEC.RData")



seg_file <- data.frame(query.cnv.UCEC$Sample, query.cnv.UCEC$Chromosome,query.cnv.UCEC$Start, query.cnv.UCEC$End, 
                       query.cnv.UCEC$Num_Probes, query.cnv.UCEC$Segment_Mean)

colnames(seg_file) <- c("Sample","Chromosome","Start", "End","Num_Probes","Segment_Mean")
write.table(seg_file, "UCEC_total.nocnv_grch38.seg", quote = F, row.names = F, sep = "\t")
