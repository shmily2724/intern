
# 1. Doc file XML
xml_file <- "/media/shmily/writable/working/data/Kinh_all_8reg.res/Kinh_all_8reg.xml"
txt <- paste(readLines(xml_file, warn = FALSE, encoding = "latin1"), collapse = "\n")

# 2. Ham lay block giua 2 tag
get_block <- function(text, tag) {
  pattern <- paste0("<", tag, "[^>]*>([\\s\\S]*?)</", tag, ">")
  m <- regexec(pattern, text, perl = TRUE)
  x <- regmatches(text, m)[[1]]
  
  if (length(x) < 2) {
    stop(paste("Khong tim thay block:", tag))
  }
  
  x[2]
}

# 3. Lay labels
lab_block <- get_block(txt, "pairDistPopLabels")
lab_lines <- unlist(strsplit(lab_block, "\n"))
lab_lines <- trimws(lab_lines)
lab_lines <- lab_lines[lab_lines != ""]
lab_lines <- lab_lines[grepl("^[0-9]+:", lab_lines)]
labs <- sub("^[0-9]+:\\s*", "", lab_lines)

print(labs)

# 4. Ham doc ma tran tam giac duoi
parse_lower_matrix <- function(block_text, n) {
  lines <- unlist(strsplit(block_text, "\n"))
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  mat <- matrix("", nrow = n, ncol = n)
  
  for (ln in lines) {
    parts <- strsplit(ln, "\\s+")[[1]]
    
    if (length(parts) >= 2 && grepl("^[0-9]+$", parts[1])) {
      row_i <- as.integer(parts[1])
      vals <- parts[-1]
      
      for (j in seq_along(vals)) {
        mat[row_i, j] <- vals[j]
      }
    }
  }
  
  mat
}

# 5. Lay FST matrix
fst_block <- get_block(txt, "PairFstMat")
fst_mat <- parse_lower_matrix(fst_block, length(labs))

# 6. Lay P-value matrix
pval_block <- get_block(txt, "PairFstPvalMat")
pval_mat <- parse_lower_matrix(pval_block, length(labs))

print(fst_mat)
print(pval_mat)

# 7. Ghep matrix cuoi
final_mat <- matrix("", nrow = length(labs), ncol = length(labs))

for (i in seq_along(labs)) {
  for (j in seq_along(labs)) {
    if (i == j) {
      final_mat[i, j] <- "0.00000"
    } else if (i > j) {
      fst_val <- fst_mat[i, j]
      pval_val <- gsub("\\+-", " ± ", pval_mat[i, j])
      final_mat[i, j] <- paste0(fst_val, " (", pval_val, ")")
    } else {
      final_mat[i, j] <- ""
    }
  }
}

rownames(final_mat) <- labs
colnames(final_mat) <- labs

print(final_mat, quote = FALSE)

# 8. Ghi file
out_file <- "/media/shmily/writable/working/result/fst_matrix.tsv"
write.table(final_mat, file = out_file, sep = "\t", quote = FALSE, col.names = NA)

cat("\nDa ghi file:", out_file, "\n")
