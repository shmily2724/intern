# ==============================================================================
# SCRIPT XỬ LÝ FST & P-VALUE TỪ ARLEQUIN XML SANG GGCORRPLOT
# ==============================================================================

library(XML)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggcorrplot) # Thêm thư viện ggcorrplot

# --- 1. ĐỌC DỮ LIỆU & LẤY NHÃN POPULATION THÔNG MINH ---
arqxml <- xmlParse("/media/shmily/writable/working/data/Kinh_all_8reg.res/Kinh_all_8reg.xml")

raw_labels <- xpathSApply(arqxml, "//pairDistPopLabels", xmlValue)[[1]]
label_lines <- strsplit(raw_labels, "\n")[[1]]

valid_label_lines <- grep("^\\s*\\d+:\\s+", label_lines, value = TRUE)
arq_rowname <- sub("^\\s*\\d+:\\s+", "", valid_label_lines)
arq_rowname <- trimws(arq_rowname)

N_pop <- length(arq_rowname)

# --- 2. XỬ LÝ MA TRẬN FST ---
raw_fst <- xpathSApply(arqxml, "//PairFstMat", xmlValue)[[1]]

fst_df <- read.table(text = raw_fst, fill = TRUE, stringsAsFactors = FALSE)
fst_mat <- as.matrix(fst_df)

arq_fstmat_mt <- fst_mat[(nrow(fst_mat) - N_pop + 1):nrow(fst_mat), 
                         (ncol(fst_mat) - N_pop + 1):ncol(fst_mat)]

mode(arq_fstmat_mt) <- "numeric"

arq_fstmat_mt[upper.tri(arq_fstmat_mt)] <- t(arq_fstmat_mt)[upper.tri(arq_fstmat_mt)]
diag(arq_fstmat_mt) <- 0
arq_fstmat_mt[arq_fstmat_mt < 0] <- 0 

rownames(arq_fstmat_mt) <- arq_rowname
colnames(arq_fstmat_mt) <- arq_rowname

# --- 3. XỬ LÝ MA TRẬN P-VALUE ---
raw_pval <- xpathSApply(arqxml, "//PairFstPvalMat", xmlValue)[[1]]

raw_pval_clean <- gsub("\\s*\\+-\\s*[0-9.]+", "", raw_pval)

pval_df <- read.table(text = raw_pval_clean, fill = TRUE, stringsAsFactors = FALSE)
pval_mat <- as.matrix(pval_df)

arq_fstp <- pval_mat[(nrow(pval_mat) - N_pop + 1):nrow(pval_mat), 
                     (ncol(pval_mat) - N_pop + 1):ncol(pval_mat)]

mode(arq_fstp) <- "numeric"

arq_fstp[upper.tri(arq_fstp)] <- t(arq_fstp)[upper.tri(arq_fstp)]
diag(arq_fstp) <- 0

rownames(arq_fstp) <- arq_rowname
colnames(arq_fstp) <- arq_rowname

# ==============================================================================
# --- 4. XUẤT ĐỒ THỊ BẰNG GGCORRPLOT (DẠNG TAM GIÁC CŨ MỊN MÀNG) ---
# ==============================================================================

p <- ggcorrplot(
  corr = arq_fstmat_mt,
  type = "lower",               # Dạng tam giác dưới
  outline.color = "black",      # Kẻ viền đen sắc nét
  p.mat = arq_fstp,             # Add thẳng ma trận p-value để đóng dấu X
  sig.level = 0.05,             
  insig = "pch",                
  pch.cex = 4,                  
  digits = 5                    
) + 
  # Nạp thẳng dải mã màu chuẩn: Trắng tinh -> Xanh nhạt -> Xanh lam -> Xanh Navy đậm
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#C6DBEF", "#4292C6", "#084594"), 
    limits = c(0, 0.01)
  ) +
  theme(panel.grid = element_blank())
# Hiển thị đồ thị
print(p)
