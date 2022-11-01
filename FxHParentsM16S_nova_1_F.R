#list.files() # make sure what we think is here is actually here

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("dada2")

library(dada2)
library(tidyverse)

## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("ParentsM_16S_ALL_novaseq_samplesF.txt", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_R1_trimmed.fq.gz")
# and one with the reverse
reverse_reads <- paste0(samples, "_R2_trimmed.fq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq.gz")
names(filtered_forward_reads) <- samples
names(filtered_reverse_reads) <- samples

plotQualityProfile(forward_reads[155:155])

plotQualityProfile(reverse_reads[150:155])

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=170, truncLen=c(210,190), multithread=10)
filtered_out


set.seed(100)
err_forward_reads <- learnErrors(filtered_forward_reads, nbases = 1e8, multithread=10, randomize=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, nbases = 1e8, multithread=10, randomize=TRUE)

plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)


ddsF <- vector("list", length(samples))
names(ddsF) <- samples
ddsR <- vector("list", length(samples))
names(ddsR) <- samples
merge_seqs <- vector("list", length(samples))
names(merge_seqs) <- samples

for(sam in samples) {
  cat("Processing:", sam, "\n")
  derep_forward <- derepFastq(filtered_forward_reads[[sam]])
  derep_reverse <- derepFastq(filtered_reverse_reads[[sam]])
  ddsF[[sam]] <- dada(derep_forward, err=err_forward_reads, multithread=10)
  ddsR[[sam]] <- dada(derep_reverse, err=err_reverse_reads, multithread=10)
  merge_seqs[[sam]] <- mergePairs(ddsF[[sam]], derep_forward, ddsR[[sam]], derep_reverse, trimOverhang = T)
}

class(ddsF)
sapply(ddsF, class)
head(names(filtered_forward_reads))
head(names(filtered_reverse_reads))
head(names( merge_seqs))
seqtab <- makeSequenceTable(merge_seqs)
class(seqtab) # matrix
dim(seqtab)




table(nchar(colnames(seqtab)))

seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10)
write_rds(seqtab_nochim, "asv_table16ParentsMnova.rds")

#seqtab_nochim <- read_rds("asv_table16ParentsMnova.rds")
  
table(nchar(colnames(seqtab_nochim)))
seqtab_nochim2 <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% 250:256]
write_rds(seqtab_nochim2, "asv_table16ParentsM2novaF.rds")

#tax <- assignTaxonomy(seqtab_nochim, "silva_nr_v132_train_set.fa.gz", multithread=1)
#write_rds(tax, "tax16ParentsMnova.rds")

tax <- assignTaxonomy(seqtab_nochim2, "silva_nr_v132_train_set.fa.gz", multithread=10)
write_rds(tax, "tax16ParentsMnova2F.rds")
