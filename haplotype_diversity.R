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
