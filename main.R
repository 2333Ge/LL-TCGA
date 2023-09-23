# 安装和加载所需的包
if (!require('rjson'))
  install.packages("rjson")
library('rjson')

# 样本文件存放目录，（注意：可以先放少部分样本跑跑验证一下，太多了需要跑很久）
DIR = './gdc-temp/'

# 读取JSON文件
json_data <- fromJSON(file = "./metadata.cart.2023-09-21.json")

#=======得到样本文件分别属于哪个分组

# 患病组文件列表
disease_files = list()
# 正常组文件列表
control_files = list()

# 循环处理
for (item in json_data) {
  entity_submitter_id <-
    item$associated_entities[[1]]$entity_submitter_id
  
  # 判断患病组和对照组（注意：这是从[文档](https://www.jianshu.com/p/f74601e8b379)中判断的分组条件，真实的分组是否这样需要再次确认下）
  group <-
    ifelse(substr(entity_submitter_id, 14, 15) < 10, "Disease", "Control")
  
  full_file_path <- paste0(DIR, item$file_id, '/', item$file_name)
  
  if (file.exists(full_file_path)) {
    if (group == 'Disease') {
      disease_files = append(disease_files, full_file_path)
    } else{
      control_files = append(control_files, full_file_path)
    }
  }
}
#=======

get_mean = function(file_list, key) {
  temp_sum = data.frame()
  for (file in file_list) {
    temp_table <-
      read.table(file,
                 header = TRUE,
                 sep = "\t",
                 fill = TRUE)
    # 提取指定列
    temp_column <- temp_table[[key]]
    
    # 将值为0的元素替换为NA
    temp_column[temp_column == 0] <- NA
    
    # 将 "temp_column" 列添加到结果数据框中，按行分组相加
    if (length(temp_sum) != 0) {
      temp_sum <- cbind(temp_sum, temp_column)
    } else {
      temp_sum <- temp_column
    }
  }
  result_mean <-
    rowMeans(temp_sum, na.rm = TRUE)
  return(result_mean)
}

# # 患病组处理，得到understand均值
# disease_result = data.frame()
# for (file in disease_files) {
#   temp_df <-
#     read.table(file,
#                header = TRUE,
#                sep = "\t",
#                fill = TRUE)
#   # 提取  "unstranded" 列
#   unstranded <- temp_df[["unstranded"]]
# 
#   # 将值为0的元素替换为NA
#   unstranded[unstranded == 0] <- NA
# 
#   # 将 "unstranded" 列添加到结果数据框中，按行分组相加
#   if (length(disease_result) != 0) {
#     disease_result <- cbind(disease_result, unstranded)
#   } else {
#     disease_result <- unstranded
#   }
# }
# 
# disease_result_mean <-
#   rowMeans(disease_result, na.rm = TRUE)  # 计算每一行的平均值
# 
# # 同上，正常组处理，得到understand均值
# control_result = data.frame()
# for (file in control_files) {
#   temp_df <-
#     read.table(file,
#                header = TRUE,
#                sep = "\t",
#                fill = TRUE)
#   # 提取  "unstranded" 列
#   unstranded <- temp_df[["unstranded"]]
# 
#   # 将值为0的元素替换为NA
#   unstranded[unstranded == 0] <- NA
# 
#   # 将 "unstranded" 列添加到结果数据框中，按行分组相加
#   if (length(control_result) != 0) {
#     control_result <- cbind(control_result, unstranded)
#   } else {
#     control_result <- unstranded
#   }
# }
# control_result_mean <-
#   rowMeans(control_result, na.rm = TRUE)  # 计算每一行的平均值
# 

disease_result_mean = get_mean(disease_files, 'understanded')
control_result_mean = get_mean(control_files, 'understanded')
# 计算logfc
log_frame1 <- log2(as.matrix(disease_result_mean))
log_frame2 <- log2(as.matrix(control_result_mean))
logfc <- log_frame1 - log_frame2


# 其他需要保留的行信息
temp_file =  read.table(disease_files[[1]],
                        header = TRUE,
                        sep = "\t",
                        fill = TRUE)
gene_name = temp_file[['gene_name']]
gene_id = temp_file[['gene_id']]



# 注意：从chatgpt代码review还提到，需要统计检验或标准化处理，暂时不理解这个概念，数据可能要进一步处理

# 合并成新的数据框得到目标结果
result_final <- data.frame(gene_id,
                           gene_name,
                           disease_result_mean,
                           control_result_mean,
                           logfc)
