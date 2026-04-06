library(dplyr)
library(tidyr)
library(openxlsx)
round_to_100 <- function(freq, digits = 2) {
  freq_round <- round(freq, digits)
  diff <- round(100 - sum(freq_round), digits)
  i <- which.max(freq_round)
  freq_round[i] <- round(freq_round[i] + diff, digits)
  return(freq_round)
}
hg <- read.delim("/media/shmily/writable/working/haplogrep2/Kinh_all_rcrs.txt")
hg2 <- hg %>% 
  select(SampleID, Haplogroup)
sheet2<- hg2 %>% 
  count(Haplogroup, name = "Count") %>% 
  mutate(
    `Frequency (%)`= round_to_100(Count / sum(Count) * 100, 2)
    ) %>% 
  arrange(desc(Count), Haplogroup)
View(sheet2)

wb <- loadWorkbook("Kinh_overview.xlsx")
addWorksheet(wb, "HG_frequency")
writeData(wb, "HG_frequency", sheet2)
saveWorkbook(wb, "Kinh_overview.xlsx", overwrite = TRUE)
  
 # Tạo bảng có Population từ SampleID
hg2 <- hg %>%
  select(SampleID, Haplogroup) %>%
  mutate(
    Population = case_when(
      grepl("^Kinh", SampleID) ~ "North",
      grepl("^KICN", SampleID) ~ "Central",
      grepl("^HG", SampleID) ~ "South",
      TRUE ~ "Unknown"
    )
  )

# Kiểm tra xem có sample nào chưa được gán không
hg2 %>% filter(Population == "Unknown")

# Tạo bảng frequency theo region
sheet3 <- hg2 %>%
  filter(Population != "Unknown") %>%
  count(Population, Haplogroup, name = "Count") %>%
  group_by(Population) %>%
  mutate(Frequency = round_to_100(Count / sum(Count) * 100, 2)) %>%
  ungroup() %>%
  select(Haplogroup, Population, Frequency) %>%
  pivot_wider(
    names_from = Population,
    values_from = Frequency,
    values_fill = 0
  ) %>%
  select(Haplogroup, North, Central, South)

colnames(sheet3) <- c("Haplogroup", "Kinh North", "Kinh Central", "Kinh South")

sheet3 <- sheet3 %>%
  arrange(Haplogroup)
# Add vào workbook có sẵn
wb <- loadWorkbook("Kinh_overview.xlsx")
addWorksheet(wb, "HG_frequency_by_region")
writeData(wb, "HG_frequency_by_region", sheet3)
saveWorkbook(wb, "Kinh_overview.xlsx", overwrite = TRUE)
