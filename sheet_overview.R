library(dplyr)
library(tidyr)
library(openxlsx)

# =========================
# 1. Đọc file haplogroup
# =========================
hg <- read.delim("Kinh_all_rcrs.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

hg2 <- hg %>%
  select(SampleID, Haplogroup)

# =========================
# 2. Gán population + reference
# =========================
hg2 <- hg2 %>%
  mutate(
    Population = case_when(
      grepl("^Kinh", SampleID) ~ "North",
      grepl("^KICN", SampleID) ~ "Central",
      grepl("^HG", SampleID) ~ "South",
      TRUE ~ "Unknown"
    ),
    Reference = case_when(
      Population == "North" ~ "Study 2018",
      Population == "Central" ~ "Present study",
      Population == "South" ~ "1kG",
      TRUE ~ NA
    ),
    Country = "Vietnam"
  )

# kiểm tra số mẫu theo vùng
print(table(hg2$Population))

# nếu có sample không khớp
print(hg2 %>% filter(Population == "Unknown"))

# =========================
# 3. Đếm haplogroup theo vùng
# =========================
hg_count <- hg2 %>%
  count(Country, Reference, Population, Haplogroup, name = "Count")

# =========================
# 4. Pivot wide
# =========================
sheet1_hg <- hg_count %>%
  pivot_wider(
    names_from = Haplogroup,
    values_from = Count,
    values_fill = 0
  )

# =========================
# 5. Thêm sample size
# =========================
sample_size <- hg2 %>%
  count(Country, Reference, Population, name = "Sample size")

sheet1 <- sample_size %>%
  left_join(sheet1_hg, by = c("Country", "Reference", "Population"))

# =========================
# 6. Sắp xếp thứ tự dòng
# =========================
sheet1$Population <- factor(sheet1$Population, levels = c("North", "Central", "South"))

sheet1 <- sheet1 %>%
  arrange(Population)

# =========================
# 7. Xuất Excel
# =========================
wb <- createWorkbook()
addWorksheet(wb, "Overview")
writeData(wb, "Overview", sheet1)

saveWorkbook(wb, "Kinh_overview.xlsx", overwrite = TRUE)
