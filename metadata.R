library(readxl)
library(openxlsx)

# =========================
# 1. Khai bao duong dan file
# =========================
haplogroup_file <- "/media/shmily/writable/working/haplogrep2/Kinh_all_rcrs.txt"
haplotype_file  <- "/media/shmily/writable/working/data/Kinh_all_haplotype_result.xlsx"
output_file     <- "/media/shmily/writable/working/data/Kinh_all_metadata.meta"

clean_colnames <- function(df) {
  names(df) <- tolower(trimws(names(df)))
  names(df) <- gsub("\\s+", "_", names(df))
  names(df) <- gsub("\\.+", "_", names(df))
  names(df) <- gsub("[^a-z0-9_]", "", names(df))
  df
}

find_col <- function(df, candidates) {
  nms <- names(df)
  hit <- nms[nms %in% candidates]
  if (length(hit) > 0) return(hit[1])
  
  for (cand in candidates) {
    hit2 <- nms[grepl(cand, nms, ignore.case = TRUE)]
    if (length(hit2) > 0) return(hit2[1])
  }
  
  NA
}

strip_quotes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('^"|"$', "", x)
  x
}

# doc haplogroup txt
haplogroup_df <- read.table(
  haplogroup_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

haplogroup_df <- clean_colnames(haplogroup_df)

id_col_hg <- find_col(haplogroup_df, c("id", "sample", "sampleid", "sample_id", "samples"))
hg_col    <- find_col(haplogroup_df, c("haplogroup", "hg", "group"))

haplogroup_df2 <- haplogroup_df[, c(id_col_hg, hg_col)]
names(haplogroup_df2) <- c("id", "haplogroup")

haplogroup_df2$id <- strip_quotes(haplogroup_df2$id)
haplogroup_df2$haplogroup <- strip_quotes(haplogroup_df2$haplogroup)

# doc haplotype xlsx
haplotype_df <- read.xlsx(haplotype_file, sheet = 1)
haplotype_df <- clean_colnames(haplotype_df)

id_col_ht <- find_col(haplotype_df, c("id", "sample", "sampleid", "sample_id", "samples"))
ht_col    <- find_col(haplotype_df, c("haplotype", "ht"))

haplotype_df2 <- haplotype_df[, c(id_col_ht, ht_col)]
names(haplotype_df2) <- c("id", "haplotype")

haplotype_df2$id <- strip_quotes(haplotype_df2$id)
haplotype_df2$haplotype <- strip_quotes(haplotype_df2$haplotype)

# clean id de merge
haplogroup_df2$id_clean <- tolower(trimws(haplogroup_df2$id))
haplotype_df2$id_clean  <- tolower(trimws(haplotype_df2$id))

meta <- merge(
  haplotype_df2,
  haplogroup_df2[, c("id_clean", "haplogroup")],
  by = "id_clean",
  all = TRUE
)

if ("id.x" %in% names(meta) && "id.y" %in% names(meta)) {
  meta$id <- ifelse(!is.na(meta$id.x) & meta$id.x != "", meta$id.x, meta$id.y)
}

id_lower <- tolower(meta$id)
meta$region <- ifelse(grepl("^kinh", id_lower), "north",
                      ifelse(grepl("^kicn", id_lower), "central",
                             ifelse(grepl("^hg", id_lower), "south", NA)))

meta_final <- meta[, c("id", "haplogroup")]
names(meta_final) <- c("ID", "Haplogroup")

meta_final[] <- lapply(meta_final, as.character)
meta_final[] <- lapply(meta_final, as.character)

write.table(
  meta_final,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  na = ""
)

cat("Da luu file:", output_file, "\n")
