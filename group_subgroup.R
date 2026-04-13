library(readxl)
library(openxlsx)

# =========================
# 1. Đọc file
# =========================
df <- read.delim("/media/shmily/writable/working/haplogrep2/Kinh_all_rcrs.txt", sep = "\t", stringsAsFactors = FALSE)

print(names(df))

# =========================
# 2. Lấy cột haplogroup
# =========================
haplo_df <- data.frame(
  haplogroup = df[["Haplogroup"]],
  stringsAsFactors = FALSE
)

haplo_df <- haplo_df[!is.na(haplo_df$haplogroup) & haplo_df$haplogroup != "", , drop = FALSE]
haplo_df$haplogroup <- trimws(haplo_df$haplogroup)
haplo_df$haplogroup <- gsub("[^A-Za-z0-9]", "", haplo_df$haplogroup)

# =========================
# 3. Hàm tách level
# =========================
split_haplogroup_levels <- function(hg) {
  sapply(1:nchar(hg), function(i) substr(hg, 1, i))
}

# =========================
# 4. Tách tất cả level
# =========================
all_levels_list <- lapply(haplo_df$haplogroup, split_haplogroup_levels)

# =========================
# 5. Đếm số sample mỗi level
# =========================
all_levels <- unlist(all_levels_list)

level_count <- as.data.frame(table(all_levels), stringsAsFactors = FALSE)
names(level_count) <- c("level", "count")

# =========================
# 6. Tạo sheet1
# =========================
max_depth <- max(sapply(all_levels_list, length))

sheet1_list <- list()

for (i in seq_along(all_levels_list)) {
  levs <- all_levels_list[[i]]
  row_list <- list()
  
  for (j in seq_along(levs)) {
    row_list[[paste0("level_", j)]] <- levs[j]
    row_list[[paste0("n_", j)]] <- level_count$count[level_count$level == levs[j]]
  }
  
  if (length(levs) < max_depth) {
    for (j in (length(levs) + 1):max_depth) {
      row_list[[paste0("level_", j)]] <- NA
      row_list[[paste0("n_", j)]] <- NA
    }
  }
  
  sheet1_list[[i]] <- as.data.frame(row_list, stringsAsFactors = FALSE)
}

sheet1 <- do.call(rbind, sheet1_list)
sheet1$full_haplogroup <- haplo_df$haplogroup
sheet1 <- sheet1[!duplicated(sheet1$full_haplogroup), ]

ordered_cols <- c("full_haplogroup")
for (j in 1:max_depth) {
  ordered_cols <- c(ordered_cols, paste0("level_", j), paste0("n_", j))
}
sheet1 <- sheet1[, ordered_cols]
# =========================
# 7. Xuất Excel
# =========================
wb <- createWorkbook()
addWorksheet(wb, "sheet1")
writeData(wb, "sheet1", sheet1)

saveWorkbook(wb, "sheet1_macro_to_smallest_subgroup.xlsx", overwrite = TRUE)
