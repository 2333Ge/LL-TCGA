# 通过[DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 得到logfc和p值

if (!require(stringr))
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
    rowCount = nrow(temp_table)
    if (ncol(simples) == 0) {
      simples <-
        data.frame(matrix(nrow = rowCount , ncol = 0))
      # 一个样本中gen_name有重复的，所以用id，但是看官方没有示例用的gene_name,可能提取过程有误
      gene_ids <- temp_table[['gene_id']]
      rownames(simples) = temp_table[['gene_id']]
    }
    simples[[entity_submitter_id]] = temp_table[['unstranded']]
  }
}
#=================

dds <-
  DESeqDataSetFromMatrix(countData = simples,
                         colData = groups,
                         # 这个字段作用不是很懂
                         
                         design = ~ group)
dds <- DESeq(dds)
# 得到dds后想看其他的数据可以在控制台执行，会快一点
res = results(dds)
resLFC <-
  lfcShrink(dds, coef = "group_control_vs_case", type = "apeglm")
