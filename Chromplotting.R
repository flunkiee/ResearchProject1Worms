# read the gtf file, this probably explodes if you dont put a gtf in here 
gtf_file <- "A.australiensis.gtf"
gtf_data <- read.table(gtf_file, sep = "\t", header = FALSE)

#useful ones start end chrom
gtf_annot <- gtf_data[, c(1, 4, 5)] 
colnames(gtf_annot) <- c("Chrom", "Start", "End")


#TRIM1 -- really good comment, past me
gtf_annot$Chrom <- gsub(".*(_chromosome(?:[:_])+[xy0-9]+|_scaffold_[0-9]+).*", "\\1", 
                        gtf_annot$Chrom, ignore.case = TRUE)
gtf_annot$Chrom <- ifelse(gtf_annot$Chrom == "OZ176061.1_Paragordius_varius_genome_assembly,_organelle:_mitochondrion", 
                   "mitochondrion", gtf_annot$Chrom)

# remove rows where Chrom starts with "unn"
gtf_annot <- gtf_annot[!grepl("^un", gtf_annot$Chrom, ignore.case = TRUE), ]


head(gtf_annot$Chrom)

print(df)

head(gtf_annot)

#call chromplot
library(chromPlot)

chromPlot(
  annot1 = gtf_annot,
  plotRndchr = T,
  bin = 500000,
  chr = NULL #i dont think i need to do this actually but im scared to remove it since it works now
)

#this code is a mess, do not use if you have too many scaffolds as it will explode 
