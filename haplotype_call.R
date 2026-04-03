library(ape)
library(pegas)
library(openxlsx)

# tên bộ dữ liệu
sample_name <- "sample name"

# đọc fasta
dna <- read.dna(paste0(sample_name, ".fasta"), format = "fasta")

# gọi haplotype
haps <- haplotype(dna)

# lấy tên sample thật từ header fasta
sample_names <- rownames(dna)
hap_index <- attr(haps, "index")

# tạo bảng sample - haplotype
hap_assign <- do.call(
  rbind,
  lapply(seq_along(hap_index), function(i) {
    data.frame(
      sample_name = sample_names[hap_index[[i]]],
      haplotype = rownames(haps)[i],
      stringsAsFactors = FALSE
    )
  })
)

rownames(hap_assign) <- NULL
hap_assign <- hap_assign[order(hap_assign$haplotype, hap_assign$sample_name), ]

# in ra để kiểm tra
print(hap_assign)

# xuất excel
wb <- createWorkbook()
addWorksheet(wb, "sample_haplotype")
writeData(wb, "sample_haplotype", hap_assign)

saveWorkbook(
  wb,
  paste0(sample_name, "_haplotype_result.xlsx"),
  overwrite = TRUE
)
