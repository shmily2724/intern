library(ape)
library(dplyr)
library(tidyr)
library(openxlsx)

# =========================
# 1. Đọc alignment FASTA
# =========================
fasta_file <- "/media/shmily/writable/working/data/Kinh_all_8reg.fasta"

dna <- read.dna(fasta_file, format = "fasta")

# Kiểm tra
print(dim(dna))       # số mẫu x số site
print(rownames(dna)[1:10])

# =========================
# 2. Chuyển mỗi sequence thành chuỗi ký tự
# =========================
seq_mat <- as.character(dna)
seq_string <- apply(seq_mat, 1, paste0, collapse = "")

hap_df <- data.frame(
  SampleID = rownames(dna),
  Sequence = seq_string,
  stringsAsFactors = FALSE
)

# =========================
# 3. Gán population
#    (tạm dùng rule cũ của bạn)
# =========================
hap_df <- hap_df %>%
  mutate(
    Population = case_when(
      grepl("^Kinh", SampleID) ~ "North",
      grepl("^KICN", SampleID) ~ "Central",
      grepl("^HG", SampleID) ~ "South",
      TRUE ~ "Unknown"
    )
  )

print(table(hap_df$Population))
print(hap_df %>% filter(Population == "Unknown"))

# =========================
# 4. Gán nhãn haplotype
#    Cùng sequence = cùng haplotype
# =========================
hap_lookup <- hap_df %>%
  distinct(Sequence) %>%
  mutate(Haplotype = paste0("H", row_number()))

hap_df <- hap_df %>%
  left_join(hap_lookup, by = "Sequence")

# Xem vài dòng
print(head(hap_df[, c("SampleID", "Population", "Haplotype")]))

# =========================
# 5. Hàm tính haplotype diversity
# =========================
calc_hd <- function(hap_vector) {
  n <- length(hap_vector)
  if (n <= 1) return(NA_real_)

  freq <- table(hap_vector) / n
  hd <- (n / (n - 1)) * (1 - sum(freq^2))
  return(hd)
}

# =========================
# 6. Tính toàn bộ dataset
# =========================
overall_n <- nrow(hap_df)
overall_h <- n_distinct(hap_df$Haplotype)
overall_hd <- calc_hd(hap_df$Haplotype)

overall_summary <- data.frame(
  Group = "All",
  Sample_size = overall_n,
  Number_of_haplotypes = overall_h,
  Haplotype_diversity = overall_hd
)

print(overall_summary)

# =========================
# 7. Tính theo population
# =========================
pop_summary <- hap_df %>%
  filter(Population != "Unknown") %>%
  group_by(Population) %>%
  summarise(
    Sample_size = n(),
    Number_of_haplotypes = n_distinct(Haplotype),
    Haplotype_diversity = calc_hd(Haplotype),
    .groups = "drop"
  )

print(pop_summary)

# =========================
# 8. Frequency của từng haplotype theo vùng
# =========================
hap_freq_by_pop <- hap_df %>%
  filter(Population != "Unknown") %>%
  count(Population, Haplotype, name = "Count") %>%
  group_by(Population) %>%
  mutate(Frequency = Count / sum(Count)) %>%
  ungroup() %>%
  arrange(Population, desc(Count))

print(hap_freq_by_pop)

# =========================
# 9. Danh sách mẫu thuộc từng haplotype
# =========================
hap_members <- hap_df %>%
  group_by(Haplotype, Population) %>%
  summarise(
    N = n(),
    Samples = paste(SampleID, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Haplotype, Population)

print(hap_members)

# =========================
# 10. Xuất Excel
# =========================
wb <- createWorkbook()

addWorksheet(wb, "Overall_summary")
writeData(wb, "Overall_summary", overall_summary)

addWorksheet(wb, "Population_summary")
writeData(wb, "Population_summary", pop_summary)

addWorksheet(wb, "Haplotype_table")
writeData(wb, "Haplotype_table", hap_df[, c("SampleID", "Population", "Haplotype")])

addWorksheet(wb, "Haplotype_freq_by_pop")
writeData(wb, "Haplotype_freq_by_pop", hap_freq_by_pop)

addWorksheet(wb, "Haplotype_members")
writeData(wb, "Haplotype_members", hap_members)

saveWorkbook(wb, "Kinh_haplotype_summary.xlsx", overwrite = TRUE)

#===============================================================================
library(ape)
library(pegas)
library(readxl)
library(openxlsx)

# =========================
# 1. Khai bao file
# =========================
fasta_file    <- "Kinh_all_mafft.fasta"
metadata_file <- "Kinh_all_metadata.xlsx"
output_file   <- "/media/shmily/writable/working/result/Kinh_diversity_result.xlsx"

# =========================
# 2. Doc FASTA
# =========================
dna <- read.dna(fasta_file, format = "fasta")

cat("So sample:", nrow(dna), "\n")
cat("So vi tri:", ncol(dna), "\n")
cat("Ten sample dau tien:\n")
print(head(rownames(dna)))

# =========================
# 3. Doc metadata
# =========================
meta <- read_excel(metadata_file)

# chuan hoa ten cot
names(meta) <- tolower(trimws(names(meta)))

# doi ten cot neu can
# yeu cau metadata co 2 cot: sample, population
meta$id <- trimws(as.character(meta$id))
meta$region <- trimws(as.character(meta$region))

names(meta)[names(meta) == "id"] <- "sample"
names(meta)[names(meta) == "region"] <- "population"

# =========================
# 4. Kiem tra sample co khop khong
# =========================
dna_samples <- rownames(dna)

cat("So sample trong FASTA:", length(dna_samples), "\n")
cat("So sample trong metadata:", nrow(meta), "\n")

missing_in_meta <- setdiff(dna_samples, meta$sample)
missing_in_dna  <- setdiff(meta$sample, dna_samples)

if (length(missing_in_meta) > 0) {
  cat("Sample co trong FASTA nhung khong co trong metadata:\n")
  print(missing_in_meta)
}

if (length(missing_in_dna) > 0) {
  cat("Sample co trong metadata nhung khong co trong FASTA:\n")
  print(missing_in_dna)
}

# chi giu lai sample co mat o ca 2 ben
common_samples <- intersect(dna_samples, meta$sample)

dna  <- dna[common_samples, ]
meta <- meta[match(common_samples, meta$sample), ]

# =========================
# 5. Ham tinh thong ke cho 1 population
# =========================
calc_pop_stats <- function(dna_pop, pop_name) {
  
  n_samples <- nrow(dna_pop)
  
  # so unique haplotypes
  haps <- haplotype(dna_pop)
  n_haplotypes <- nrow(haps)
  
  # haplotype diversity
  H <- hap.div(dna_pop)
  
  # nucleotide diversity
  pi <- nuc.div(dna_pop)
  
  # MPD = mean pairwise differences
  # dist.dna model = "N" la so khac biet nucleotide
  d <- dist.dna(dna_pop, model = "N", pairwise.deletion = TRUE)
  MPD <- mean(as.numeric(d))
  
  # variance cua pi (neu co)
  pi_var <- NA
  if (!is.null(attr(pi, "variance"))) {
    pi_var <- attr(pi, "variance")
  }
  
  # variance cua H (neu co)
  H_var <- NA
  if (!is.null(attr(H, "variance"))) {
    H_var <- attr(H, "variance")
  }
  
  data.frame(
    population = pop_name,
    sample_size = n_samples,
    unique_haplotypes = n_haplotypes,
    haplotype_diversity_H = as.numeric(H),
    H_variance = as.numeric(H_var),
    nucleotide_diversity_pi = as.numeric(pi),
    pi_variance = as.numeric(pi_var),
    MPD = as.numeric(MPD)
  )
}

# =========================
# 6. Tinh cho tung population
# =========================
pop_list <- unique(meta$population)

result_list <- list()

for (pop in pop_list) {
  samples_pop <- meta$sample[meta$population == pop]
  dna_pop <- dna[samples_pop, ]
  
  cat("Dang tinh population:", pop, "- n =", nrow(dna_pop), "\n")
  
  result_list[[pop]] <- calc_pop_stats(dna_pop, pop)
}

result_df <- do.call(rbind, result_list)

# =========================
# 7. Tinh % difference from mean cho H va pi
# =========================
mean_H  <- mean(result_df$haplotype_diversity_H, na.rm = TRUE)
mean_pi <- mean(result_df$nucleotide_diversity_pi, na.rm = TRUE)

result_df$H_percent_diff_from_mean  <- ((result_df$haplotype_diversity_H - mean_H) / mean_H) * 100
result_df$pi_percent_diff_from_mean <- ((result_df$nucleotide_diversity_pi - mean_pi) / mean_pi) * 100

# =========================
# 8. Xuat Excel
# =========================
wb <- createWorkbook()

addWorksheet(wb, "diversity_stats")
writeData(wb, "diversity_stats", result_df)

saveWorkbook(wb, output_file, overwrite = TRUE)

cat("Da xuat file:", output_file, "\n")
