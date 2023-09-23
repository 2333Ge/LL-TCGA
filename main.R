# 安装和加载所需的包
if (!require('rjson'))
  install.packages("rjson")
library('rjson')

# 样本文件存放目录，（注意：可以先放少部分样本跑跑验证一下，太多了需要跑很久）
DIR = './workspace/gdc-temp/'

# 读取JSON文件
json_data <-
  fromJSON(file = "./workspace/metadata.cart.2023-09-21.json")

#=======  得到样本文件分别属于哪个分组  =======

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
#=================

#=======  计算logfc  =======

# 定义一个方法，求指定文件指定列的平均值
get_mean = function(file_list, key) {
  temp_sum = data.frame()
  for (file in file_list) {
    temp_table <- read.table(file,
                             header = TRUE,
                             sep = "\t",
                             fill = TRUE)
    print(key)
    print(temp_table[[key]])
    print(temp_table[['unstranded']])
    # 提取指定列
    temp_column <- temp_table[[key]]
    
    # 将值为0的元素替换为NA
    temp_column[temp_column == 0] <- NA
    
    # 按列求和
    if (length(temp_sum) != 0) {
      temp_sum <- cbind(temp_sum, temp_column)
    } else {
      temp_sum <- temp_column
    }
  }
  # 得到平局值返回
  return(rowMeans(temp_sum, na.rm = TRUE))
}

# 得到患病组均值
disease_result_mean = get_mean(disease_files, 'unstranded')
# 得到对照组均值
control_result_mean = get_mean(control_files, 'unstranded')
# 计算logfc
logfc <- log2(as.matrix(disease_result_mean)) - log2(as.matrix(control_result_mean))

#=================

# 其他需要保留的行信息，每个样本一样的，随便读一个文件
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
