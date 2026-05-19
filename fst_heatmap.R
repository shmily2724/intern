# ==============================================================================
# SCRIPT XỬ LÝ FST & P-VALUE TỪ ARLEQUIN XML SANG CORRPLOT (BULLETPROOF VERSION)
# ==============================================================================

library(XML)
library(corrplot)
library(magrittr)
library(dplyr)

# --- 1. ĐỌC DỮ LIỆU & LẤY NHÃN POPULATION THÔNG MINH ---
arqxml <- xmlParse("/media/shmily/writable/working/data/Kinh_all_8reg.res/Kinh_all_8reg.xml")

# Lấy text thô của nhãn
raw_labels <- xpathSApply(arqxml, "//pairDistPopLabels", xmlValue)[[1]]
label_lines <- strsplit(raw_labels, "\n")[[1]]

# CHỐNG LỖI 1: Tự động lọc đúng các dòng chứa tên Pop (định dạng " 1: PopA")
valid_label_lines <- grep("^\\s*\\d+:\\s+", label_lines, value = TRUE)
arq_rowname <- sub("^\\s*\\d+:\\s+", "", valid_label_lines) # Cắt bỏ phần "1: "
arq_rowname <- trimws(arq_rowname)

N_pop <- length(arq_rowname) # Lưu lại số lượng Pop để cắt ma trận chuẩn xác


# --- 2. XỬ LÝ MA TRẬN FST ---
raw_fst <- xpathSApply(arqxml, "//PairFstMat", xmlValue)[[1]]

# Đọc thành dạng bảng tự bù đắp khoảng trống
fst_df <- read.table(text = raw_fst, fill = TRUE, stringsAsFactors = FALSE)
fst_mat <- as.matrix(fst_df)

# CHỐNG LỖI 2: Cắt chính xác khối ma trận N x N chứa toàn số ở góc dưới cùng phải
arq_fstmat_mt <- fst_mat[(nrow(fst_mat) - N_pop + 1):nrow(fst_mat), 
                         (ncol(fst_mat) - N_pop + 1):ncol(fst_mat)]

# Ép kiểu an toàn (Vì bây giờ ma trận 100% chỉ còn là các con số text)
mode(arq_fstmat_mt) <- "numeric"

# Xây dựng Full Matrix (Gắn tam giác dưới lên tam giác trên)
arq_fstmat_mt[upper.tri(arq_fstmat_mt)] <- t(arq_fstmat_mt)[upper.tri(arq_fstmat_mt)]
diag(arq_fstmat_mt) <- 0
arq_fstmat_mt[arq_fstmat_mt < 0] <- 0 # Xóa Fst âm

rownames(arq_fstmat_mt) <- arq_rowname
colnames(arq_fstmat_mt) <- arq_rowname


# --- 3. XỬ LÝ MA TRẬN P-VALUE ---
raw_pval <- xpathSApply(arqxml, "//PairFstPvalMat", xmlValue)[[1]]

# CHỐNG LỖI 3: Dùng Regex dọn sạch sẽ dấu +- và độ lệch chuẩn SD (Ví dụ: " +- 0.000")
raw_pval_clean <- gsub("\\s*\\+-\\s*[0-9.]+", "", raw_pval)

pval_df <- read.table(text = raw_pval_clean, fill = TRUE, stringsAsFactors = FALSE)
pval_mat <- as.matrix(pval_df)

# Cắt khối P-value N x N
arq_fstp <- pval_mat[(nrow(pval_mat) - N_pop + 1):nrow(pval_mat), 
                     (ncol(pval_mat) - N_pop + 1):ncol(pval_mat)]

# Ép kiểu an toàn
mode(arq_fstp) <- "numeric"

# Xây dựng Full Matrix
arq_fstp[upper.tri(arq_fstp)] <- t(arq_fstp)[upper.tri(arq_fstp)]
diag(arq_fstp) <- 0

rownames(arq_fstp) <- arq_rowname
colnames(arq_fstp) <- arq_rowname


# --- 4. XUẤT ĐỒ THỊ CORRPLOT ---
# Vẽ biểu đồ với matrix FST và mark điểm P-value không có ý nghĩa thống kê
corrplot(arq_fstmat_mt, method = "color", type = "lower", na.label = " ",
         tl.col = "black", tl.srt = 45, order = "FPC", cl.lim = c(0, 1),
         p.mat = arq_fstp, sig.level = 0.05, insig = "pch", 
         pch.cex = 1.2, pch.col = "#6b717a")
