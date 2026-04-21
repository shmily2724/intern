setwd("/media/shmily/writable/working/result")

library(readxl)
library(dplyr)
library(stringr)
library(plotly)

file_path <- "Kinh_all_haplotype_result.xlsx"

# =========================================================
# 1) Doc du lieu
# =========================================================
hg <- read_excel(file_path, sheet = "HG_frequency") %>%
  transmute(
    haplogroup = str_squish(Haplogroup),
    n = as.numeric(Count)
  ) %>%
  filter(!is.na(haplogroup), !is.na(n), n > 0)
# =========================================================
# 1) Doc du lieu cho kinh trung
# =========================================================
#tmp <- read_excel(file_path, sheet = "HG_frequency_by_region") %>%
#  rename_with(str_squish)
#
#central_col <- names(tmp)[str_detect(names(tmp), regex("central", ignore_case = TRUE))][1]

#hg <- tmp %>%
#  transmute(
#    haplogroup = str_squish(Haplogroup),
#    n = as.numeric(.data[[central_col]])
#  ) %>%
#  filter(!is.na(haplogroup), !is.na(n), n > 0)

# =========================================================
# 2) Map moi haplogroup vao backbone path
#    R nam trong N
#    F nam trong R9
#    M10/11/12/... gom vao M*
# =========================================================
hg <- hg %>%
  mutate(
    path = case_when(
      # ----- N > R > B -----
      str_detect(haplogroup, "^B4")  ~ "N|R|B|B4",
      str_detect(haplogroup, "^B5")  ~ "N|R|B|B5",
      str_detect(haplogroup, "^B6")  ~ "N|R|B|B6",
      
      # ----- N > R > R9 > F -----
      str_detect(haplogroup, "^F1")  ~ "N|R|R9|F|F1",
      str_detect(haplogroup, "^F2")  ~ "N|R|R9|F|F2",
      str_detect(haplogroup, "^F3")  ~ "N|R|R9|F|F3",
      
      # direct R9 lineages
      str_detect(haplogroup, "^R9")  ~ "N|R|R9",
      
      # ----- N > R > R* -----
      str_detect(haplogroup, "^R11") ~ "N|R|R*|R11",
      str_detect(haplogroup, "^R22") ~ "N|R|R*|R22",
      str_detect(haplogroup, "^R(\\*|\\+|$)") ~ "N|R|R*",
      
      # ----- N direct branches -----
      str_detect(haplogroup, "^A")   ~ "N|A",
      str_detect(haplogroup, "^N9")  ~ "N|N9",
      str_detect(haplogroup, "^Y1")  ~ "N|Y|Y1",
      
      # ----- M branches -----
      str_detect(haplogroup, "^C4")  ~ "M|C|C4",
      str_detect(haplogroup, "^C7")  ~ "M|C|C7",
      
      str_detect(haplogroup, "^D4")  ~ "M|D|D4",
      str_detect(haplogroup, "^D5")  ~ "M|D|D5",
      
      str_detect(haplogroup, "^G2")  ~ "M|G|G2",
      
      # direct M branches
      str_detect(haplogroup, "^M7")  ~ "M|M7",
      str_detect(haplogroup, "^M9")  ~ "M|M9",
      
      # M* bucket for rare M lineages
      str_detect(haplogroup, "^M10") ~ "M|M*|M10",
      str_detect(haplogroup, "^M11") ~ "M|M*|M11",
      str_detect(haplogroup, "^M12") ~ "M|M*|M12",
      str_detect(haplogroup, "^M23") ~ "M|M*|M23",
      str_detect(haplogroup, "^M51") ~ "M|M*|M51",
      str_detect(haplogroup, "^M59") ~ "M|M*|M59",
      str_detect(haplogroup, "^M71") ~ "M|M*|M71",
      str_detect(haplogroup, "^M74") ~ "M|M*|M74",
      
      TRUE ~ NA_character_
    )
  )

# Neu con haplogroup nao chua map thi dung lai de sua them
if (any(is.na(hg$path))) {
  print(hg %>% filter(is.na(path)) %>% count(haplogroup, sort = TRUE))
  stop("Co haplogroup chua duoc map vao backbone.")
}

# =========================================================
# 3) Thu tu backbone muon hien thi
# =========================================================
path_order <- c(
  "N|R|B|B4",
  "N|R|B|B5",
  "N|R|B|B6",
  
  "N|R|R9",
  "N|R|R9|F|F1",
  "N|R|R9|F|F2",
  "N|R|R9|F|F3",
  
  "N|R|R*",
  "N|R|R*|R11",
  "N|R|R*|R22",
  
  "N|A",
  "N|N9",
  "N|Y|Y1",
  
  "M|C|C4",
  "M|C|C7",
  "M|D|D4",
  "M|D|D5",
  "M|G|G2",
  "M|M7",
  "M|M9",
  "M|M*|M10",
  "M|M*|M11",
  "M|M*|M12",
  "M|M*|M23",
  "M|M*|M51",
  "M|M*|M59",
  "M|M*|M71",
  "M|M*|M74"
)

paths <- hg %>%
  group_by(path) %>%
  summarise(direct = sum(n), .groups = "drop") %>%
  mutate(path = factor(path, levels = path_order)) %>%
  arrange(path)

total_n <- sum(paths$direct)

# =========================================================
# 4) Tao node tree cho plotly
#    dung branchvalues = 'remainder'
#    direct = gia tri tai node do
#    total  = tong ca nhanh con
# =========================================================
prefixes_of <- function(path_string) {
  x <- strsplit(path_string, "\\|")[[1]]
  vapply(seq_along(x), function(i) paste(x[1:i], collapse = "|"), character(1))
}

node_ids <- unique(unlist(lapply(as.character(paths$path), prefixes_of)))

nodes <- tibble(id = node_ids) %>%
  mutate(
    parent = if_else(str_detect(id, "\\|"), sub("\\|[^|]+$", "", id), ""),
    label  = sub("^.*\\|", "", id),
    depth  = str_count(id, "\\|") + 1L
  )

direct_map <- setNames(paths$direct, as.character(paths$path))

nodes <- nodes %>%
  mutate(
    direct = unname(ifelse(id %in% names(direct_map), direct_map[id], 0)),
    total = vapply(
      id,
      function(node) {
        sum(paths$direct[
          as.character(paths$path) == node |
            startsWith(as.character(paths$path), paste0(node, "|"))
        ])
      },
      numeric(1)
    ),
    percent = total / total_n * 100
  )

# =========================================================
# 5) Mau
#    N la nhom lon, R la nhanh con trong N
# =========================================================
col_map <- c(
  "M"   = "#C86A50",
  "C"   = "#D69A89",
  "C4"  = "#E4B8AD",
  "C7"  = "#F0D4CC",
  "D"   = "#CF7E69",
  "D4"  = "#DA8F7B",
  "D5"  = "#E5A79A",
  "G"   = "#E7C9C0",
  "G2"  = "#F3E0DA",
  "M7"  = "#C25B47",
  "M9"  = "#D98F7D",
  "M*"  = "#B84F3B",
  "M10" = "#C96D59",
  "M11" = "#CD7663",
  "M12" = "#D27F6C",
  "M23" = "#D78977",
  "M51" = "#DD9383",
  "M59" = "#E29E8F",
  "M71" = "#E8AA9C",
  "M74" = "#EEB8AB",
  
  "N"   = "#148DB5",
  "A"   = "#A8D5EE",
  "N9"  = "#7FC3E8",
  "Y"   = "#C7E4F2",
  "Y1"  = "#DDEEF8",
  
  "R"   = "#0F79A9",
  "B"   = "#59BE6B",
  "B4"  = "#8ECF98",
  "B5"  = "#70C17E",
  "B6"  = "#B8E1BE",
  
  "R9"  = "#33A6CF",
  "F"   = "#E09B8E",
  "F1"  = "#E7AAA0",
  "F2"  = "#EDC2BB",
  "F3"  = "#F3D8D3",
  
  "R*"  = "#89CFE8",
  "R11" = "#BCE2F2",
  "R22" = "#DCEFF8"
)

nodes <- nodes %>%
  mutate(color = if_else(label %in% names(col_map), col_map[label], "#D9D9D9"))

# =========================================================
# 6) Xac dinh node nao la internal node / leaf
# =========================================================
nodes <- nodes %>%
  mutate(
    has_child = id %in% parent[parent != ""]
  )

# label het tat ca
nodes <- nodes %>%
  mutate(
    nodes <- nodes %>%
      mutate(
        text_show = paste0(label, "<br>", sprintf("%.1f%%", percent))
      ),
    hover_txt = paste0(
      "<b>", label, "</b>",
      "<br>Direct count: ", direct,
      "<br>Subtree total: ", total,
      "<br>Percent of all samples: ", sprintf("%.2f%%", percent),
      "<br>Path: ", id
    )
  )

# =========================================================
# 7) Ve
# =========================================================
p <- plot_ly(
  data = nodes,
  ids = ~id,
  labels = ~label,
  parents = ~parent,
  values = ~direct,
  type = "sunburst",
  branchvalues = "remainder",
  sort = FALSE,
  text = ~text_show,
  textinfo = "text",
  insidetextorientation = "radial",
  hovertext = ~hover_txt,
  hoverinfo = "text",
  marker = list(
    colors = nodes$color,
    line = list(color = "white", width = 2)
  ),
  maxdepth = 6
) %>%
  layout(
    margin = list(l = 20, r = 20, b = 20, t = 20),
    font = list(family = "Arial", size = 16)
  )

p
