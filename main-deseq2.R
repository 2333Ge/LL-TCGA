# 通过[DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 得到logfc和p值

if(!require(stringr))
  install.packages("stringr")  # 常规R包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require(DESeq2))
  BiocManager::install("DESeq2")  # 差异分析R包
if (!require(apeglm))
  BiocManager::install("apeglm")

library('stringr')
library('DESeq2')

DIR = './workspace/gdc-temp/'

# 读取JSON文件
json_data <-
  fromJSON(file = "./workspace/metadata.cart.2023-09-21.json")

#=======  样本分组信息  =======


groups = data.frame(entity_submitter_id = character(), group = character())
simples = data.frame()

for (item in json_data) {
  entity_submitter_id <-
    item$associated_entities[[1]]$entity_submitter_id
  
  # 判断患病组和对照组（注意：这是从[文档](https://www.jianshu.com/p/f74601e8b379)中判断的分组条件，真实的分组是否这样需要再次确认下）
  group <-
    ifelse(substr(entity_submitter_id, 14, 15) < 10, "case", "control")
  
  full_file_path <- paste0(DIR, item$file_id, '/', item$file_name)
  
  
  if (file.exists(full_file_path)) {
    groups = rbind(groups, data.frame(entity_submitter_id, group))
    temp_table <- read.table(full_file_path,
                             header = TRUE,
                             sep = "\t",
                             fill = TRUE)
    if (length(simples) == 0) {
      simples = data.frame(gene_name = temp_table[['gene_name']])
      simples[[entity_submitter_id]] = temp_table[['unstranded']]
      # 如何动态指定列名
      # simples = data.frame(`entity_submitter_id` = temp_table[['gene_name']])
      
    } else{
      simples[[entity_submitter_id]] = temp_table[['unstranded']]
    }
    
  }
}
#=================
# print(ncol(simples))
# print(nrow(groups))

dds <- DESeqDataSetFromMatrix(countData = simples[, 2:ncol(simples)],
                              colData = groups,
                              design = ~ group)

dds <- DESeq(dds)
# res <- results(dds)
# summary(res)
resLFC <- lfcShrink(dds, coef = "group_control_vs_case", type = "apeglm")
resLFC
# resSig <- subset(res, padj < 0.05)
