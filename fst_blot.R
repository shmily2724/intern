library(XML)
library(corrplot)

xml_file <- "/media/shmily/writable/working/data/Kinh_all_8reg.res/Kinh_all_8reg.xml"
stopifnot(file.exists(xml_file))

arqxml <- xmlParse(xml_file)

# labels
lab_nodes <- xpathSApply(arqxml, "//pairDistPopLabels", xmlValue)
stopifnot(length(lab_nodes) >= 1)

labs <- unlist(strsplit(lab_nodes[1], "\n"))
labs <- trimws(labs)
labs <- labs[labs != ""]
labs <- labs[!grepl("^Label", labs, ignore.case = TRUE)]
labs <- labs[!grepl("^-+", labs)]
labs <- gsub("^\\d+:\\s*", "", labs)

print(labs)

# fst block
fst_nodes <- xpathSApply(arqxml, "//PairFstMat", xmlValue)
stopifnot(length(fst_nodes) >= 1)

fst_lines <- unlist(strsplit(fst_nodes[1], "\n"))
fst_lines <- trimws(fst_lines)
fst_lines <- fst_lines[fst_lines != ""]

# bỏ dòng "Distance method..."
fst_lines <- fst_lines[!grepl("^Distance method", fst_lines)]

# chỉ giữ các dòng dữ liệu thật, không lấy dòng header "1 2 3"
fst_lines <- fst_lines[grepl("^\\d+\\s+-?\\d*\\.\\d+", fst_lines)]

print(fst_lines)

vals_list <- lapply(strsplit(fst_lines, "\\s+"), function(x) as.numeric(x[-1]))

print(vals_list)

stopifnot(length(labs) == length(vals_list))

n <- length(labs)
mat <- matrix(NA_real_, nrow = n, ncol = n)

for (i in seq_len(n)) {
  mat[i, 1:length(vals_list[[i]])] <- vals_list[[i]]
}

mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
diag(mat) <- 0

rownames(mat) <- labs
colnames(mat) <- labs

mat_plot <- mat
mat_plot[mat_plot < 0] <- 0
mat_plot[!is.finite(mat_plot)] <- 0

print(mat)
print(mat_plot)

corrplot(mat_plot,
         method = "color",
         type = "lower",
         tl.col = "black",
         tl.srt = 45,
         order = "original",
         is.corr = FALSE,
         na.label = " ")
