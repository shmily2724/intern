
# ape, stringr were already installed

install.packages("adegenet")
install.packages("devtools")
install.packages("tidyverse")
library(adegenet)
library(devtools)
library(tidyverse)
library(ape)
library(stringr)
install_github("helixcn/phylotools", build_vignettes = TRUE)
library(phylotools)

Exclude8regions = function(infile, outfile) {
  
  tabular = data.frame(read.dna(infile, format = "fasta", as.character = T))
  
  tabular_new <- tabular[1,!grepl("-", tabular[1,])]
  
  start_pattern <- c("aatttccaccaaa", # 303-315
                     "ccatcctacccag", # 514-522
                     "aaccccaaagaca", # 568-573
                     "gtgttttagatca", # 956-965
                     "ccattttacctca", # 5895-5899
                     "tttaccctatagc", # 8271-8289
                     "caatccacatcaa"  # 16182-16193
  )
  
  
  end_pattern <- c("gcttctggccaca", # 303-315
                   "nncgctgctaacc", # 514-522
                   "acagtttatgtag", # 568-573
                   "aataaagctaaaa", # 956-965
                   "actgatgttcgcc", # 5895-5899
                   "gagcccactgtaa", # 8271-8289
                   "atgcttacaagca"  # 16182-16193
  )
  
  for (i in 1:length(start_pattern)) {
    start <- 0
    for (j in 1:(length(tabular_new)-nchar(start_pattern[1])+1)) {
      if (apply(tabular_new[,j:(j+nchar(start_pattern[i])-1)], 1 , paste , collapse = "" ) == start_pattern[i]) {
        start <- colnames(tabular_new[(j+nchar(start_pattern[i]))])
      }
    }
    
    end <- 0
    for (j in 1:(length(tabular_new)-nchar(end_pattern[1])+1)) {
      if (apply(tabular_new[,j:(j+nchar(end_pattern[i])-1)], 1 , paste , collapse = "" ) == end_pattern[i]) {
        end <- colnames(tabular_new[(j-1)])
      }
    }
    
    tabular[,gsub("\\D", "", start):gsub("\\D", "", end)] <- "?"
  }
  
  tabularSeq = data.frame(paste(">", rownames(tabular), sep = ""))
  tabularSeq[,2] <- tabular %>% unite(tabular, colnames(tabular), sep = "")
  
  enter = "\n"
  
  #delete if any existing file 
  
  unlink(outfile, force=TRUE)
  
  for (i in 2:nrow(tabularSeq)) {
    cat(tabularSeq[i,1], file = outfile, append = TRUE, enter)
    cat(tabularSeq[i,2], file = outfile, append = TRUE, enter)
  }
  
}

Exclude8regions(infile = "/Population_genetics_Evol_Mol/mtDNA_2/Processed_sequence_dataset/VN_CN_2522_RSRS_Thao_aln.fasta",
                outfile = "/Population_genetics_Evol_Mol/mtDNA_2/Processed_sequence_dataset/VN_CN_2522_RSRS_Thao_aln_8r.fasta")


Exclude8regions(infile = "Processed_sequence_dataset/CN_VN_mtDNA_2531_reorder_RSRS_Dinh_aln.fasta",
                outfile = "Processed_sequence_dataset/CN_VN_mtDNA_2531_reorder_RSRS_Dinh_aln_8r.fasta")

Exclude8regions(infile = "Processed_sequence_dataset/VN_CN_mtDNA_2522_RSRS_reorder_aln.fasta",
                outfile = "Processed_sequence_dataset/VN_CN_mtDNA_2522_RSRS_reorder_aln_8r.fasta")
