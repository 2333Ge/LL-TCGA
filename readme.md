# git 说明

用git管理代码，方便保存、回退、共享

环境配置：https://www.runoob.com/git/git-install-setup.html

拉取代码：控制台执行

```
git clone [仓库上方的clone地址]
```

存代码
```sh
git add .
git commit -m "保存信息"
```


# R 语言环境配置

- 环境：https://www.runoob.com/r/r-setup.html
- 开发工具：https://posit.co/download/rstudio-desktop/

# 代码说明

## 目录说明

```
├── main.R // 核心代码，暂时弃用
├── main-deseq2.R // v2版，使用deseq2进行分析
├── main.ipynb // 暂时用不到
├── readme.md
└── workspace // 放样本数据
    ├── gdc // 压缩包解压后放到了这
    ├── gdc-temp // 取少部分压缩包内的数据，调试用，写代码时看结果快一点
    ├── gdc_manifest_20230921_011847.txt
    └── metadata.cart.2023-09-21.json
```
## 执行代码

RStudio 导入这个目录后选择 `main.R` 文件，点击文件上方的 `run`

控制台执行`View(result_final)`可视化查看最终结果

# 其他资料

其他似乎有助于数据分析的资料

- https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html#TCGAanalyze_Preprocessing:_Preprocessing_of_Gene_Expression_data_(IlluminaHiSeq_RNASeqV2)

- https://rpubs.com/tiagochst/TCGAbiolinks_to_DESEq2

- https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

- https://www.math.pku.edu.cn/teachers/lidf/docs/Rbook/html/_Rbook/intro.html

- https://www.huber.embl.de/msmb/index.html

- https://mform.ecology.wang/