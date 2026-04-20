library(readxl)
library(dplyr)
library(stringr)
library(plotly)

file_path <- "Kinh_all_haplotype_result.xlsx"

# 1) Đọc bảng frequency tổng
hg <- read_excel(file_path, sheet = "HG_frequency") %>%
  transmute(
    haplogroup = str_squish(Haplogroup),
    n = as.numeric(Count)
  ) %>%
  filter(!is.na(haplogroup), !is.na(n), n > 0)

# 2) Gán subgroup và macro group
hg <- hg %>%
  mutate(
    subgroup = case_when(
      str_detect(haplogroup, "^B")  ~ "B",
      str_detect(haplogroup, "^F")  ~ "F",
      str_detect(haplogroup, "^R")  ~ "R",
      
      str_detect(haplogroup, "^A")  ~ "A",
      str_detect(haplogroup, "^N9") ~ "N9",
      str_detect(haplogroup, "^Y")  ~ "Y",
      str_detect(haplogroup, "^N")  ~ "N*",
      
      str_detect(haplogroup, "^C")  ~ "C",
      str_detect(haplogroup, "^D")  ~ "D",
      str_detect(haplogroup, "^G")  ~ "G",
      str_detect(haplogroup, "^Z")  ~ "Z",
      str_detect(haplogroup, "^M7") ~ "M7",
      str_detect(haplogroup, "^M8") ~ "M8",
      str_detect(haplogroup, "^M9") ~ "M9",
      str_detect(haplogroup, "^M")  ~ "M*",
      
      TRUE ~ "Other"
    ),
    macro = case_when(
      subgroup %in% c("C", "D", "G", "Z", "M7", "M8", "M9", "M*") ~ "M",
      subgroup %in% c("A", "N9", "Y", "N*") ~ "N",
      subgroup %in% c("B", "F", "R") ~ "R",
      TRUE ~ "Other"
    )
  )

# 3) Tùy chọn: gộp các haplogroup quá nhỏ để hình đỡ rối
# bỏ comment nếu muốn
 hg <- hg %>%
  mutate(
   haplogroup_plot = if_else(n / sum(n) < 0.01,
                              paste0(subgroup, "_other"),
                              haplogroup)
  )
# nếu không gộp thì dùng haplogroup gốc
#hg <- hg %>%
#  mutate(haplogroup_plot = haplogroup)

# 4) Tạo bảng node cho sunburst
nodes_macro <- hg %>%
  count(macro, wt = n, name = "value") %>%
  transmute(
    ids = macro,
    labels = macro,
    parents = "",
    value = value
  )

nodes_sub <- hg %>%
  count(macro, subgroup, wt = n, name = "value") %>%
  transmute(
    ids = paste(macro, subgroup, sep = "|"),
    labels = subgroup,
    parents = macro,
    value = value
  )

nodes_hg <- hg %>%
  count(macro, subgroup, haplogroup_plot, wt = n, name = "value") %>%
  transmute(
    ids = paste(macro, subgroup, haplogroup_plot, sep = "|"),
    labels = haplogroup_plot,
    parents = paste(macro, subgroup, sep = "|"),
    value = value
  )

nodes <- bind_rows(nodes_macro, nodes_sub, nodes_hg)

# 5) Vẽ
p <- plot_ly(
  data = nodes,
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~value,
  type = "sunburst",
  branchvalues = "total",
  textinfo = "label+percent entry",
  insidetextorientation = "radial",
  hovertemplate = paste0(
    "<b>%{label}</b><br>",
    "Count: %{value}<br>",
    "Percent: %{percentEntry:.2%}<extra></extra>"
  )
) %>%
  layout(
    margin = list(l = 20, r = 20, b = 20, t = 20)
  )

p
