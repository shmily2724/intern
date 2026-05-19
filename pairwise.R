library(ape)
library(dplyr)
library(openxlsx)

# =========================
# 1. Đọc alignment
# =========================
fasta_file <- "/media/shmily/writable/working/data/Kinh_all_8reg.fasta"
dna <- read.dna(fasta_file, format = "fasta")

# Metadata mẫu
meta <- data.frame(
  SampleID = rownames(dna),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Population = case_when(
      grepl("^Kinh", SampleID) ~ "North",
      grepl("^KICN", SampleID) ~ "Central",
      grepl("^HG", SampleID) ~ "South",
      TRUE ~ "Unknown"
    )
  )

print(table(meta$Population))
print(meta %>% filter(Population == "Unknown"))

# =========================
# 2. Hàm tính pi và k
# =========================
calc_pi <- function(dna_sub) {
  if (nrow(dna_sub) < 2) return(NA_real_)
  d <- as.matrix(dist.dna(
    dna_sub,
    model = "raw",
    pairwise.deletion = TRUE
  ))
  mean(d[upper.tri(d)], na.rm = TRUE)
}

calc_k <- function(dna_sub) {
  if (nrow(dna_sub) < 2) return(NA_real_)
  d <- as.matrix(dist.dna(
    dna_sub,
    model = "N",
    pairwise.deletion = TRUE
  ))
  mean(d[upper.tri(d)], na.rm = TRUE)
}

# số site alignment
alignment_length <- ncol(dna)

# =========================
# 3. Tính toàn bộ dataset
# =========================
overall_summary <- data.frame(
  Group = "All",
  Sample_size = nrow(dna),
  Alignment_length = alignment_length,
  Average_pairwise_differences_k = calc_k(dna),
  Nucleotide_diversity_pi = calc_pi(dna)
)

print(overall_summary)

# =========================
# 4. Tính theo population
# =========================
pop_levels <- unique(meta$Population)
pop_levels <- pop_levels[pop_levels != "Unknown"]

pop_summary <- lapply(pop_levels, function(pop) {
  ids <- meta$SampleID[meta$Population == pop]
  dna_sub <- dna[ids, , drop = FALSE]
  
  data.frame(
    Population = pop,
    Sample_size = nrow(dna_sub),
    Alignment_length = ncol(dna_sub),
    Average_pairwise_differences_k = calc_k(dna_sub),
    Nucleotide_diversity_pi = calc_pi(dna_sub)
  )
}) %>%
  bind_rows()

print(pop_summary)

# =========================
# 5. Pairwise distance matrix (optional)
# =========================
pairwise_pi_matrix <- as.matrix(
  dist.dna(dna, model = "raw", pairwise.deletion = TRUE)
)

pairwise_k_matrix <- as.matrix(
  dist.dna(dna, model = "N", pairwise.deletion = TRUE)
)

# đổi thành data.frame để ghi Excel
pairwise_pi_df <- data.frame(SampleID = rownames(pairwise_pi_matrix), pairwise_pi_matrix, check.names = FALSE)
pairwise_k_df  <- data.frame(SampleID = rownames(pairwise_k_matrix),  pairwise_k_matrix,  check.names = FALSE)

# =========================
# 6. Xuất Excel
# =========================
wb <- loadWorkbook("Kinh_all_haplotype_result.xlsx")

addWorksheet(wb, "Overall_pi")
writeData(wb, "Overall_pi", overall_summary)

addWorksheet(wb, "Population_pi")
writeData(wb, "Population_pi", pop_summary)

addWorksheet(wb, "Pairwise_pi_matrix")
writeData(wb, "Pairwise_pi_matrix", pairwise_pi_df)

addWorksheet(wb, "Pairwise_k_matrix")
writeData(wb, "Pairwise_k_matrix", pairwise_k_df)

saveWorkbook(wb, "Kinh_all_haplotype_result.xlsx", overwrite = TRUE)
