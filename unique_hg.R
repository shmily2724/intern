library(dplyr)
library(tidyr)
library(openxlsx)

hg <- read.delim(
  "/media/shmily/writable/working/haplogrep2/Kinh_all_rcrs.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

hg2 <- hg %>%
  select(SampleID, Haplogroup) %>%
  mutate(
    Haplogroup = trimws(Haplogroup),
    Population = case_when(
      grepl("^Kinh", SampleID) ~ "North",
      grepl("^KICN", SampleID) ~ "Central",
      grepl("^HG", SampleID) ~ "South",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(Population != "Unknown")

# số mẫu mỗi vùng
sample_size_tbl <- hg2 %>%
  count(Population, name = "Sample size")

# haplogroup xuất hiện ở bao nhiêu vùng
hg_presence <- hg2 %>%
  distinct(Population, Haplogroup) %>%
  add_count(Haplogroup, name = "n_region")

# haplogroup private = chỉ xuất hiện ở 1 vùng
private_hg_tbl <- hg_presence %>%
  filter(n_region == 1) %>%
  count(Population, name = "Private haplogroup")

# bảng tóm tắt
unique_hg_table <- sample_size_tbl %>%
  left_join(private_hg_tbl, by = "Population") %>%
  mutate(`Private haplogroup` = replace_na(`Private haplogroup`, 0)) %>%
  pivot_longer(
    cols = c(`Sample size`, `Private haplogroup`),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Population,
    values_from = Value
  )

print(unique_hg_table)

wb <- loadWorkbook("Kinh_all_haplotype_result.xlsx")
addWorksheet(wb, "Unique_HG_summary")
writeData(wb, "Unique_HG_summary", unique_hg_table)
saveWorkbook(wb, "Kinh_all_haplotype_result.xlsx", overwrite = TRUE)
