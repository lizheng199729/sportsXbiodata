library(shiny)
library(shinydashboard)
library(DT)
library(xtable)
library(ggplot2)
library(ggpubr)
library(tidyverse)      # 已包含 dplyr, tidyr, stringr, tibble, readr, purrr
library(latex2exp)
library(rstatix)
library(ggsignif)
library(RColorBrewer)
library(Hmisc)
library(ggrepel)
library(shinyjs)
library(GSEABase)
library(tidyr)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ActivePathways)
library(clusterProfiler)
library(enrichplot)
select <- dplyr::select
#导入数据
tissueexpdata <-read.csv(tissueexpdata1,row.names = 1, header = T)
pyGENEdata <- read.table(pyGENEdata1, header = T,row.names = 1, check.names = F, sep = "\t") 
py <- read.csv(py1, header = T,row.names = 1, check.names = F) 
genedata <- read.csv(genedata1, header = T, check.names = F,row.names = 1)
immunedata <- read.csv(immunedata1, header = T, check.names = F,row.names = 1)
Proteindata <- read.csv(Proteindata1,row.names = 1, header = T)
cortex <-read.csv(cortex1,header = T,row.names = 1)
gastrocnemius <-read.csv(gastrocnemius1,header = T,row.names = 1)
heart <-read.csv(heart1,header = T,row.names = 1)
kidney <-read.csv(kidney1,header = T,row.names = 1)
liver <-read.csv(liver1,header = T,row.names = 1)
lung <-read.csv(lung1,header = T,row.names = 1)
white <-read.csv(white1,header = T,row.names = 1)
GEOmouse <-read.csv(GEOmouse1,header = T,row.names = 1)
genedata <-read.csv(genedata1,header = T,row.names = 1)
difexpdata <- read.csv(difexpdata1,header = T,row.names = 1)
HumanDifferentialGenes <- read.csv(HumanDifferentialGenes1, header = T, check.names = F,row.names = 1)
mouseGEOdata <- read.csv(mouseGEOdata)
df_split <- read.csv(df_split)
keywords <- read.csv(keywords)
terma <- read.csv(terma)
RNA_vena <- read.csv(RNA_vena1, header = T, row.names = 1,check.names = F) 
RNA_White <- read.csv(RNA_White1, header = T, row.names = 1,check.names = F) 
RNA_Brown <- read.csv(RNA_Brown1, header = T, row.names = 1,check.names = F) 
RNA_Liver <- read.csv(RNA_Liver1, header = T, row.names = 1,check.names = F) 
RNA_Small <- read.csv(RNA_Small1, header = T, row.names = 1,check.names = F)
RNA_Lung <- read.csv(RNA_Lung1, header = T, row.names = 1,check.names = F)
RNA_Ovaries <- read.csv(RNA_Ovaries1, header = T, row.names = 1,check.names = F)
RNA_Testes <- read.csv(RNA_Testes1, header = T, row.names = 1,check.names = F)
RNA_Spleen <- read.csv(RNA_Spleen1, header = T, row.names = 1,check.names = F)
RNA_Colon <- read.csv(RNA_Colon1, header = T, row.names = 1,check.names = F)
RNA_Adrenals <- read.csv(RNA_Adrenals1, header = T, row.names = 1,check.names = F)
RNA_Kidney <- read.csv(RNA_Kidney1, header = T, row.names = 1,check.names = F)
RNA_Heart <- read.csv(RNA_Heart1, header = T, row.names = 1,check.names = F)
RNA_Vastus <- read.csv(RNA_Vastus1, header = T, row.names = 1,check.names = F)
RNA_Gastrocnemius <- read.csv(RNA_Gastrocnemius1, header = T, row.names = 1,check.names = F)
proteomics_white <- read.csv(proteomics_white1, header = T, row.names = 1,check.names = F)
proteomics_lung <- read.csv(proteomics_lung1, header = T, row.names = 1,check.names = F)
proteomics_liver <- read.csv(proteomics_liver1, header = T, row.names = 1,check.names = F)
proteomics_heart <- read.csv(proteomics_heart1, header = T, row.names = 1,check.names = F)
proteomics_gastrocnemius <- read.csv(proteomics_gastrocnemius1, header = T, row.names = 1,check.names = F)
proteomics_cortex <- read.csv(proteomics_cortex1, header = T, row.names = 1,check.names = F)
proteomics_kidney <- read.csv(proteomics_kidney1, header = T, row.names = 1,check.names = F)
group <- read.csv(group1, header = T,check.names = F)
#导入富集分析
human <- getGmt(human)
mouse <- getGmt(mouse)
human <- geneIds(human)
human <- stack(human)
colnames(human) <- c("gene", "term")
human <- human[, c("term", "gene")]  # 转换为 clusterProfiler 需要的顺序
mouse <- geneIds(mouse)
mouse <- stack(mouse)
colnames(mouse) <- c("gene", "term")
mouse <- mouse[, c("term", "gene")]  # 转换为 clusterProfiler 需要的顺序
#读取测试数据

#导入功能
######
HumanGEO<- function(genename){
  data_matrix <- HumanDifferentialGenes[genename,]
  m<-ncol(data_matrix)
  #data_matrix <- t(data_matrix)
  #提取FC值与显著性
  padj_with_total <- grep("padj", names(data_matrix), value = TRUE)
  selected_padj <- data_matrix[, padj_with_total]
  FC_with_total <- grep("log2FoldChange", names(data_matrix), value = TRUE)
  selected_FC <- data_matrix[, FC_with_total]
  #转换长型矩阵
  data_padj <- gather(as.data.frame(t(selected_padj)),key = "DataSet",value = "padj")
  padj_underscore <- sub("\\..*", "", (colnames(selected_padj)))
  data_padj$dataSet <- rep(paste(padj_underscore))
  data_FC <- gather(as.data.frame(t(selected_FC)),key = "DataSet",value = "log2FoldChange")
  data_result <- cbind(data_padj,data_FC$log2FoldChange)
  colnames(data_result)[4] <- "log2FoldChange"
  #修改pos表示方法
  rows_count <- nrow(data_result)
  Neg_Pos <- rep("Neg",rows_count)
  Neg_Pos[which(data_result$log2FoldChange > 0)] <- "Pos"
  data_result$Neg_Pos <- Neg_Pos
  # 修改p值表示方法：
  pvalue <- rep(NA,rows_count)
  pvalue[which(data_result$padj > 0.05)] <- ">0.05"
  pvalue[which(data_result$padj < 0.05)] <- "<0.05"
  pvalue[which(data_result$padj < 0.01)] <- "<0.01"
  pvalue[which(data_result$padj < 0.001)] <- "<0.001"
  pvalue[which(data_result$padj < 0.0001)] <- "<0.0001"
  data_result$significance <- pvalue
  data_result$DataSet <- sub("-.*", "", data_result$DataSet)
  p <- ggplot(data_result, aes(dataSet, DataSet)) +
    geom_point(aes(fill = significance, size = abs(log2FoldChange)),
               color = "#999999", shape = 21) +
    scale_fill_manual(values = c("#212c5f", "#3366b1", "#42b0e4", "#7bc6ed", "#dfe1e0")) +
    geom_point(data = data_result[which(data_result$Neg_Pos == "Pos"), ],
               aes(color = significance, size = abs(log2FoldChange)),
               shape = 16) +
    scale_color_manual(values = c("#f26666", "#f49699", "#facccc", "#facccc", "#d9dbd9")) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x =  element_text(angle = 60, hjust = 1, vjust = 1),
      legend.margin = margin(20, unit = 'pt'),
      legend.position = "top",  # <<< 图例放在上方
      plot.margin = unit(c(1, 1, 0.1, 1), "cm")
    ) +
    xlab("") +
    ylab("") +
    guides(
      size = guide_legend(title = "log2FoldChange"),
      fill = guide_legend(title = expression("Adjusted P-values (padj)\nfor Downregulated Genes")),
      col = guide_legend(title = expression("Adjusted P-values (padj)\nfor Upregulated Genes"))
    )
  selected_columns <- data_result[, c(1, 2, 3, 4)]
  selected_columns <- selected_columns %>%
    select("DataSet","dataSet","padj","log2FoldChange",everything())
  colnames(selected_columns)[colnames(selected_columns) == "DataSet"] <- "gene_name"
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}
#######
mouseGEO<- function(genename){
  data_matrix <- mouseGEOdata[genename,]
  m<-ncol(data_matrix)
  #data_matrix <- t(data_matrix)
  #提取FC值与显著性
  padj_with_total <- grep("padj", names(data_matrix), value = TRUE)
  selected_padj <- data_matrix[, padj_with_total]
  FC_with_total <- grep("log2FoldChange", names(data_matrix), value = TRUE)
  selected_FC <- data_matrix[, FC_with_total]
  #转换长型矩阵
  data_padj <- gather(as.data.frame(t(selected_padj)),key = "DataSet",value = "padj")
  padj_underscore <- sub("\\..*", "", (colnames(selected_padj)))
  data_padj$dataSet <- rep(paste(padj_underscore))
  data_FC <- gather(as.data.frame(t(selected_FC)),key = "DataSet",value = "log2FoldChange")
  data_result <- cbind(data_padj,data_FC$log2FoldChange)
  colnames(data_result)[4] <- "log2FoldChange"
  #修改pos表示方法
  rows_count <- nrow(data_result)
  Neg_Pos <- rep("Neg",rows_count)
  Neg_Pos[which(data_result$log2FoldChange > 0)] <- "Pos"
  data_result$Neg_Pos <- Neg_Pos
  # 修改p值表示方法：
  pvalue <- rep(NA,rows_count)
  pvalue[which(data_result$padj > 0.05)] <- ">0.05"
  pvalue[which(data_result$padj < 0.05)] <- "<0.05"
  pvalue[which(data_result$padj < 0.01)] <- "<0.01"
  pvalue[which(data_result$padj < 0.001)] <- "<0.001"
  pvalue[which(data_result$padj < 0.0001)] <- "<0.0001"
  data_result$significance <- pvalue
  data_result$DataSet <- sub("-.*", "", data_result$DataSet)
  p <- ggplot(data_result, aes(dataSet, DataSet)) +
    geom_point(aes(fill = significance, size = abs(log2FoldChange)),
               color = "#999999", shape = 21) +
    scale_fill_manual(values = c("#212c5f", "#3366b1", "#42b0e4", "#7bc6ed", "#dfe1e0")) +
    geom_point(data = data_result[which(data_result$Neg_Pos == "Pos"), ],
               aes(color = significance, size = abs(log2FoldChange)),
               shape = 16) +
    scale_color_manual(values = c("#f26666", "#f49699", "#facccc", "#facccc", "#d9dbd9")) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x =  element_text(angle = 80, hjust = 1, vjust = 1),
      legend.margin = margin(20, unit = 'pt'),
      legend.position = "top",  # <<< 图例放在上方
      plot.margin = unit(c(1, 1, 0.1, 1), "cm")
    ) +
    xlab("") +
    ylab("") +
    guides(
      size = guide_legend(title = "log2FoldChange"),
      fill = guide_legend(title = expression("Adjusted P-values (padj)\nfor Downregulated Genes")),
      col = guide_legend(title = expression("Adjusted P-values (padj)\nfor Upregulated Genes"))
    )
  selected_columns <- data_result[, c(1, 2, 3, 4)]
  selected_columns <- selected_columns %>%
    select("DataSet","dataSet","padj","log2FoldChange",everything())
  colnames(selected_columns)[colnames(selected_columns) == "DataSet"] <- "gene_name"
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}
##########
sportgenerich <- function(data,species,Pathway){
  gene_df <- data %>%
    arrange(desc(log2FC)) %>%
    distinct(gene, .keep_all = TRUE)
  
  # 构建 GSEA 所需的 geneList 向量
  geneList <- gene_df$log2FC
  names(geneList) <- gene_df$gene
  pathway_id <- Pathway
  species_single <- species %>% filter(term == pathway_id)
  #使用 list 格式的自定义基因集
  gsea_result <- GSEA(
    geneList = geneList,
    TERM2GENE =species_single,
    verbose = FALSE,
    pvalueCutoff = 1
  )
  p <- gseaplot2(gsea_result, geneSetID = pathway_id, base_size = 14,
                 title = paste("Enrichment for", pathway_id))
  result <- list(plotData = p, tableData = gsea_result@result)
  return(result)
}
#########
graph <- function(Keywords){
  result1 <- df_split[df_split$Entry == Keywords, ]
  result2 <- result1[, c(1, 3)]
  # 1️ 创建边数据
  edges <- result2 %>% select(from = Entry,to = Molecule)
  # 2️创建图对象
  graph <- tbl_graph(edges = edges, directed = FALSE)
  
  # 3️ 增加节点度（连接数），用于大小
  # 自动标记节点类型
  graph <- graph %>%
    activate(nodes) %>%
    mutate(type = case_when(
      name %in% result2$Entry ~ "Entry",
      name %in% result2$Molecule ~ "Molecule"
    ))
  
  # 绘图
  p <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(alpha = 0.8), color = "grey") +
    geom_node_point(aes(size = 5, color = type)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void() +
    scale_color_manual(values = c("Entry" = "red", "Molecule" = "blue", "Other" = "green")) +
    guides(color = guide_legend(title = "Node Type"))
  result <- list(plotData = p, tableData = result1)
  return(result)
}
#########
user_rich <- function(Control, Experimental, RNA_upload = NULL, RNA, protein) {
  group <- paste0(Control ,"_vs_", Experimental)
  # 提取RNA列
  RNA_columns <- grep(group, colnames(RNA))
  RNA_mat <- RNA[, RNA_columns, drop = FALSE]
  
  # 提取蛋白质列
  protein_columns <- grep(group, colnames(protein))
  protein_mat <- protein[, protein_columns, drop = FALSE]
  # 合并
  common_rows <- intersect(rownames(protein_mat), rownames(RNA_mat))
  RNA_mat1 <- RNA_mat[common_rows, , drop = FALSE]
  protein_mat1 <- protein_mat[common_rows, , drop = FALSE]
  merged_mat <- cbind(RNA_mat1, protein_mat1)
  
  # 如果上传了数据
  if (!is.null(RNA_upload)) {
    names(RNA_upload)[names(RNA_upload) == "log2FoldChange"] <- "log2fc"
    names(RNA_upload)[names(RNA_upload) == "padj"] <- "FDR"
    common_rows_upload <- intersect(rownames(merged_mat), rownames(RNA_upload))
    # 对上传的RNA进行处理
    RNA_upload_mat1 <- RNA_upload[common_rows_upload, , drop = FALSE]
    merged_mat <- merged_mat[common_rows_upload, , drop = FALSE]
    
    # 上传数据也合并成一个上传版 merged
    merged_mat <- cbind(RNA_upload_mat1, merged_mat)
    
  }
  
  # 形成P值矩阵
  p_matrix <- data.frame(
    row.names = rownames(merged_mat),
    rna = merged_mat[, grep("padj", colnames(merged_mat))],
    protein = merged_mat[, grep("adj.P.Val", colnames(merged_mat))],
    user = merged_mat[,grep("FDR", colnames(merged_mat))]
  )
  p_matrix[is.na(p_matrix)] <- 1
  
  # 形成Fold Change矩阵
  FC_matrix <- data.frame(
    row.names = rownames(merged_mat),
    rna = merged_mat[, grep("log2FoldChange", colnames(merged_mat))],
    protein = merged_mat[, grep("logFC", colnames(merged_mat))],
    user = merged_mat[,grep("log2fc", colnames(merged_mat))]
  )
  dir_matrix <- sign(FC_matrix)
  dir_matrix[is.na(dir_matrix)] <- 0
  
  if (!is.null(RNA_upload)) {
    constraints_vector <- c(1, 1, 1)
  } else {
    constraints_vector <- c(1, 1)
  }
  
  p_matrix <- as.matrix(p_matrix)
  dir_matrix <- as.matrix(dir_matrix)
  
  # 合并P值
  directional_merged_pvals <- merge_p_values(p_matrix, method = "DPM", dir_matrix, constraints_vector)
  merged_pvals <- merge_p_values(p_matrix, method = "Brown")
  
  # 阈值设定
  threshold <- -log10(0.05)
  
  lineplot_df <- data.frame(
    original = -log10(merged_pvals),
    modified = -log10(directional_merged_pvals)
  )
  
  # 分类
  lineplot_df <- lineplot_df %>%
    mutate(category = ifelse(original < threshold, "No change",
                             ifelse(modified > threshold, "Joint expression", "Opposite expression")))
  
  # 画图
  p <- ggplot(lineplot_df) +
    geom_point(size = 2.4, shape = 19,
               aes(original, modified,
                   color = ifelse(original < threshold, "gray",
                                  ifelse(modified > threshold, "#1F449C", "#F05039")))) +
    labs(title = "",
         x = "Merged -log10(P)",
         y = "Directional Merged -log10(P)") +
    geom_hline(yintercept = threshold, linetype = "dashed",
               col = 'black', size = 0.5) +
    geom_vline(xintercept = threshold, linetype = "dashed",
               col = "black", size = 0.5) +
    geom_abline(size = 0.5, slope = 1, intercept = 0) +
    scale_color_identity()
  # 返回结果
  result <- list(plotData = p, table = lineplot_df)
  return(result)
}
#########
Multi_omics <- function(Control, Experimental,Threshold, RNA, protein){
  group <- paste0(Control ,"_vs_", Experimental)
  # 提取RNA列
  RNA_columns <- grep(group, colnames(RNA))
  RNA_mat <- RNA[, RNA_columns, drop = FALSE]
  colnames(RNA_mat) <- paste0("mRNA_", sub(".*\\.", "", colnames(RNA_mat)))
  # 提取蛋白质列
  protein_columns <- grep(group, colnames(protein))
  protein_mat <- protein[, protein_columns, drop = FALSE]
  colnames(protein_mat) <- paste0("protein_", sub(".*\\.", "", colnames(protein_mat)))
  colnames(protein_mat) <- sub("protein_Val", "protein_P.Val", colnames(protein_mat))
  # 合并
  common_rows <- intersect(rownames(protein_mat), rownames(RNA_mat))
  RNA_mat1 <- RNA_mat[common_rows, , drop = FALSE]
  protein_mat1 <- protein_mat[common_rows, , drop = FALSE]
  merged_mat <- cbind(RNA_mat1, protein_mat1)
  mRNA_fc <- merged_mat[["mRNA_log2FoldChange"]]
  protein_fc <- merged_mat[["protein_logFC"]]
  
  merged_mat$group <- case_when(
    mRNA_fc > Threshold & protein_fc > Threshold ~ "mRNA+Protein_both",
    mRNA_fc < -Threshold & protein_fc < -Threshold ~ "mRNA+Protein_both",
    
    abs(mRNA_fc) > Threshold & abs(protein_fc) <= Threshold ~ "mRNA_only",
    abs(protein_fc) > Threshold & abs(mRNA_fc) <= Threshold ~ "Protein_only",
    
    abs(mRNA_fc) > Threshold & abs(protein_fc) > Threshold & sign(mRNA_fc) != sign(protein_fc) ~
      ifelse(abs(mRNA_fc) >= abs(protein_fc), "mRNA_only", "Protein_only"),
    
    TRUE ~ "nonsig"
  )
  
  # 定义绘制坐标轴函数：
  p <- draw_axis_line(4, 4)
  p1 <- p + geom_point(data= merged_mat, aes(mRNA_log2FoldChange, protein_logFC, color = group))+
    scale_color_manual(values = c("mRNA+Protein_both" = "#dd8653", 
                                  "mRNA_only" = "#59a5d7", 
                                  "Protein_only" = "#aa65a4", 
                                  "#878787"),
                       breaks = c("mRNA+Protein_both","mRNA_only","Protein_only"))+
    xlab("mRNA:FC")+
    ylab("Protein:FC")+
    theme(legend.position = "bottom")+
    annotate("text", label = "", parse = TRUE, 
             x = -2, y = 2, size = 4, colour = "black")+
    guides(color = guide_legend(title = "", ncol = 1, byrow = TRUE))+
    coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4))
  result <- list(plotData = p1, table = merged_mat)
  return(result)
}

draw_axis_line <- function(length_x, length_y, 
                           tick_step = NULL, lab_step = NULL){
  axis_x_begin <- -1*length_x
  axis_x_end <- length_x
  
  axis_y_begin  <- -1*length_y
  axis_y_end    <- length_y
  
  if (missing(tick_step))
    tick_step <- 1
  
  if (is.null(lab_step))
    lab_step <- 2
  
  # axis ticks data
  tick_x_frame <- data.frame(ticks = seq(axis_x_begin, axis_x_end, 
                                         by = tick_step))
  
  tick_y_frame <-  data.frame(ticks = seq(axis_y_begin, axis_y_end, 
                                          by = tick_step))
  
  # axis labels data
  lab_x_frame <- subset(data.frame(lab = seq(axis_x_begin, axis_x_end, 
                                             by = lab_step), zero = 0), 
                        lab != 0)
  
  lab_y_frame <- subset(data.frame(lab = seq(axis_y_begin, axis_y_end,
                                             by = lab_step),zero = 0), 
                        lab != 0)
  
  # set tick length
  tick_x_length = 0.05
  tick_y_length = 0.05
  
  # set zero point
  
  data <- data.frame(x = 0, y = 0)
  p <- ggplot(data = data) +
    
    # draw axis line
    geom_segment(y = 0, yend = 0, 
                 x = axis_x_begin, 
                 xend = axis_x_end,
                 size = 0.5) + 
    geom_segment(x = 0, xend = 0, 
                 y = axis_y_begin, 
                 yend = axis_y_end,
                 size = 0.5) +
    # x ticks
    geom_segment(data = tick_x_frame, 
                 aes(x = ticks, xend = ticks, 
                     y = 0, yend = 0 - tick_x_length)) +
    # y ticks
    geom_segment(data = tick_y_frame, 
                 aes(x = 0, xend = 0 - tick_y_length, 
                     y = ticks, yend = ticks)) + 
    
    # labels
    geom_text(data=lab_x_frame, aes(x=lab, y=zero, label=lab), vjust = 1.5) +
    geom_text(data=lab_y_frame, aes(x=zero, y=lab, label=lab), hjust= 1.5) +
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.text = element_blank())
  return(p)
}
#########
pansport <- function(gene_name) {
  if (!gene_name %in% rownames(difexpdata)) {
    stop("Gene not found in the dataset.")
  }
  
  # 提取基因表达数据
  gene_data <- difexpdata[gene_name, , drop = FALSE]
  
  # 提取控制组和有氧组数据
  pre_data <- gene_data[, grep("control", colnames(gene_data)), drop = FALSE]
  post_data <- gene_data[, grep("aerobics", colnames(gene_data)), drop = FALSE]
  
  # 提取组织名称
  tissues_pre <- gsub(".*_", "", colnames(pre_data))
  tissues_post <- gsub(".*_", "", colnames(post_data))
  tissues_pre <- sub("\\..*", "", tissues_pre)
  tissues_post <- sub("\\..*", "", tissues_post)
  # 构建结果数据框
  data <- data.frame(
    sample = c(colnames(pre_data), colnames(post_data)),
    tissue = c(tissues_pre, tissues_post),
    group = rep(c("Control", "aerobics"), times = c(length(pre_data), length(post_data))),
    expr = c(t(pre_data), t(post_data))
  )
  # 按中位数由高到低排列：
  data_new <- data %>% 
    group_by(tissue) %>% 
    mutate(median = median(expr), group_max = max(expr)) %>% 
    arrange(desc(median))
  
  # 调整因子顺序：
  data_new$tissue <- factor(data_new$tissue, levels = unique(data_new$tissue))
  
  data_new1 <- data_new %>% 
    group_by(tissue) %>% 
    mutate(group_number = length(unique(group)))
  
  # 此数据用于添加显著性检验的label：
  stat.test <- data_new1[which(data_new1$group_number == 2),] %>%
    group_by(tissue) %>%
    pairwise_t_test(
      expr ~ group, paired = FALSE, 
      p.adjust.method = "fdr"
    ) %>% 
    add_xy_position(x = "tissue")
  data_new1$group <- factor(data_new1$group, levels = c("Control", "aerobics"))
  # 重新画图：
  p <- ggplot()+
    # 基础图形：
    geom_boxplot(data = data_new1, aes(x = tissue, y =expr, fill = group),outlier.shape = 21, outlier.fill = "white")+
    # 设置颜色：
    scale_fill_manual(name = NULL, values = c("Control" = "#5488ef","aerobics" = "#d6503a" ))+
    # 设置主题：
    theme_bw()+
    labs(y = bquote(.(gene_name) ~ " mRNA Level (" * .("Count") * ")"))+
    #labs(y = TeX("mRNA Level ($log_{2}$CPM)", bold = T))+
    theme(axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.99, 0.99), 
          legend.justification = c(1,1))+
    # 显著性检验：
    stat_pvalue_manual(
      stat.test, label = "p.adj.signif", 
      bracket.size = 0.3, # 粗细
      tip.length = 0.01  # 两边竖线的长度
    )
  selected_columns <- stat.test[, c(1, 3, 4, 7, 8, 9, 10)]
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}
########
tissueexp <- function(genename, tissuename,Statistical){
  data_matrix <- tissueexpdata[genename,]
  data_matrix <- data_matrix %>% dplyr::select(contains(tissuename))
  # 提取第一个单词并分组
  tissue_names <- str_extract(colnames(data_matrix), "(?<=_).*")
  tissue_names <- sub("\\..*$", "", tissue_names)
  # 将数据框转换为长格式
  data_long <- data_matrix %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "value") %>%
    mutate(x = tissue_names) %>%
    select(-sample)
  # 添加IQR和中位数信息计算：
  data_long <- data_long %>%
    group_by(x) %>%
    mutate(
      IQR = IQR(value, na.rm = TRUE),   # 在计算IQR时忽略缺失值
      Med = median(value, na.rm = TRUE) # 在计算中位数时忽略缺失值
    ) %>%
    #mutate(IQR = IQR(value), Med = median(value)) %>%
    ungroup()
  # 绘图：
  p<-ggplot(data_long, aes(x, value, fill = x))+
    # 误差棒：
    geom_errorbar(aes(ymin = Med - IQR, ymax = Med + IQR),
                  width = 0.5, linetype = "dashed")+
    # 箱线图：
    geom_boxplot(color = NA, size = 1,
                 notch = TRUE)+
    # 中心白色矩形：
    geom_rect(aes(xmin = rep(seq(0.8, 4.8, 1), length.out = nrow(data_long)),
                  xmax = rep(seq(1.2, 5.2, 1), length.out = nrow(data_long)),
                  ymin = Med-0.0001, ymax = Med+0.0001), fill = "white")+
    # 颜色模式：
    scale_fill_manual(values = c("#4087c7", "#ec0f80",
                                 "#a24b9c", "#56bc85", "#8bc53f"))+
    # 显著性检验：
    geom_signif(aes(x, value),
                # 指定你要比较的组，以list形式传入：
                comparisons = list(c("Control", "w1"),
                                   c("Control", "w2"),
                                   c("Control", "w4"),
                                   c("Control", "w8")),
                map_signif_level=T, # 展示P值或显著性
                textsize=3, # label大小
                test=Statistical, # 检验类型
                step_increase=0.1 # 位置变化
    )+
    xlab("Aerobic duration")+
    #ylab("Aerobic duration")+
    ylab(bquote(.(genename) ~ " at " ~ .(tissuename) ~ " mRNA Level (Count)"))+
    # 主题调整：
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          legend.position = "none")
  result <- list(plotData = p)
  return(result)
} 

########
Phenotypic_correlation<- function(genename){
  data_matrix <- pyGENEdata[genename,c(1:899)]
  # 确保列名唯一
  colnames(data_matrix) <- make.unique(colnames(data_matrix))
  # 提取第一个单词并分组
  tissue_names <- sub("_.*", "", colnames(data_matrix))
  tissue_names <- sub("\\..*$", "", tissue_names)
  # 将数据框转换为长格式
  long_df <- data_matrix %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    mutate(tissue = tissue_names) %>%
    select(-sample)
  
  # 将长格式数据转换为宽格式，每个组织作为一列
  wide_df <- long_df %>%
    pivot_wider(names_from = tissue, values_from = expression, values_fn = list(expression = list))
  
  # 展开列表列为多个列
  final_data <- as.data.frame(do.call(cbind, lapply(wide_df[-1], unlist)))
  cor_list <- rcorr(as.matrix(final_data),as.matrix(py))
  m<-ncol(final_data)
  n<-ncol(py)
  #提取相关性与显著性
  cor_value <- cor_list$r[1:m, (m+1):(m+n)]
  p_value <- cor_list$P[1:m, (m+1):(m+n)]
  #转换长型矩阵
  data_cor <- gather(as.data.frame(t(cor_value)),key = "tissue",value = "cor")
  data_cor$Phenotype <- rep(paste(colnames(cor_value)),m)
  data_p <- gather(as.data.frame(t(p_value)),key = "tissue",value = "pvalue")
  data_p$Phenotype <- rep(paste(colnames(p_value)),m)
  data_result <- cbind(data_cor,data_p$pvalue)
  colnames(data_result)[4] <- "pvalue"
  #修改pos表示方法
  rows_count <- nrow(data_result)
  Neg_Pos <- rep("Neg",rows_count)
  Neg_Pos[which(data_result$cor>0)] <- "Pos"
  data_result$Neg_Pos <- Neg_Pos
  # 修改p值表示方法：
  pvalue <- rep(NA,rows_count)
  pvalue[which(data_result$pvalue > 0.05)] <- ">0.05"
  pvalue[which(data_result$pvalue < 0.05)] <- "<0.05"
  pvalue[which(data_result$pvalue < 0.01)] <- "<0.01"
  pvalue[which(data_result$pvalue < 0.001)] <- "<0.001"
  pvalue[which(data_result$pvalue < 0.0001)] <- "<0.0001"
  data_result$pvalue2 <- pvalue
  p <- ggplot(data_result,aes(tissue,Phenotype))+
    # 蓝色气泡图：
    geom_point(aes(fill=pvalue2, size=abs(cor)),color = "#999999",shape=21)+
    scale_fill_manual(values = c("#212c5f","#3366b1","#42b0e4","#7bc6ed","#dfe1e0"))+
    # 红色气泡图：
    geom_point(data = data_result[which(data_result$Neg_Pos == "Pos"),],
               aes(color=pvalue2, size=abs(cor)),shape=16)+
    scale_color_manual(values = c("#f26666","#f49699","#facccc","#facccc","#d9dbd9"))+
    # 主题：
    theme_bw()+
    theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
          # 坐标轴label方向：
          axis.text.x = element_text(angle = 45, hjust = 1),
          # 图例间距：
          legend.margin = margin(20,unit = 'pt'))+
    xlab("")+
    ylab("")+
    # 图例：
    guides(size = guide_legend(title = "Spearman's "),
           fill = guide_legend(title = expression("Negtive \ncorrelation \nFDR q-value")),
           col = guide_legend(title = expression("Positive \ncorrelation \nFDR q-value")))
  
  # 返回结果
  selected_columns <- data_result[, c(1, 2, 3, 4)]
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}
###########

immgene <- function(genenames){
  data_matrix <- genedata[genenames,c(1:899)]
  data_matrix <- t(data_matrix)
  cor_list <- rcorr(as.matrix(data_matrix),as.matrix(immunedata))
  m<-ncol(data_matrix)
  n<-ncol(immunedata)
  #提取相关性与显著性
  cor_value <- cor_list$r[1:m, (m+1):(m+n)]
  p_value <- cor_list$P[1:m, (m+1):(m+n)]
  #转换长型矩阵
  data_cor <- gather(as.data.frame(t(cor_value)),key = "genename",value = "cor")
  data_cor$immune_cell <- rep(paste(colnames(cor_value)),m)
  data_p <- gather(as.data.frame(t(p_value)),key = "genename",value = "pvalue")
  data_p$immune_cell <- rep(paste(colnames(p_value)),m)
  data_result <- cbind(data_cor,data_p$pvalue)
  colnames(data_result)[4] <- "pvalue"
  #修改pos表示方法
  rows_count <- nrow(data_result)
  Neg_Pos <- rep("Neg",rows_count)
  Neg_Pos[which(data_result$cor>0)] <- "Pos"
  data_result$Neg_Pos <- Neg_Pos
  # 修改p值表示方法：
  pvalue <- rep(NA,rows_count)
  pvalue[which(data_result$pvalue > 0.05)] <- ">0.05"
  pvalue[which(data_result$pvalue < 0.05)] <- "<0.05"
  pvalue[which(data_result$pvalue < 0.01)] <- "<0.01"
  pvalue[which(data_result$pvalue < 0.001)] <- "<0.001"
  pvalue[which(data_result$pvalue < 0.0001)] <- "<0.0001"
  data_result$pvalue2 <- pvalue
  p <- ggplot(data_result,aes(immune_cell,genename))+
    # 蓝色气泡图：
    geom_point(aes(fill=pvalue2, size=abs(cor)),color = "#999999",shape=21)+
    scale_fill_manual(values = c("#212c5f","#3366b1","#42b0e4","#7bc6ed","#dfe1e0"))+
    # 红色气泡图：
    geom_point(data = data_result[which(data_result$Neg_Pos == "Pos"),],
               aes(color=pvalue2, size=abs(cor)),shape=16)+
    scale_color_manual(values = c("#f26666","#f49699","#facccc","#facccc","#d9dbd9"))+
    # 主题：
    theme_bw()+
    theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
          # 坐标轴label方向：
          axis.text.x = element_text(angle = 45, hjust = 1),
          # 图例间距：
          legend.margin = margin(20,unit = 'pt'))+
    xlab("")+
    ylab("")+
    # 图例：
    guides(size = guide_legend(title = "Spearman's"),
           fill = guide_legend(title = expression("Negtive \ncorrelation \nFDR q-value")),
           col = guide_legend(title = expression("Positive \ncorrelation \nFDR q-value")))
  selected_columns <- data_result[, c(1, 2, 3, 4)]
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}

#######
panprotein <- function(gene_name) {
  if (!gene_name %in% rownames(Proteindata)) {
    stop("Gene not found in the dataset.")
  }
  
  # 提取基因表达数据
  gene_data <- Proteindata[gene_name, , drop = FALSE]
  
  # 提取控制组和有氧组数据
  pre_data <- gene_data[, grep(".Control", colnames(gene_data)), drop = FALSE]
  post_data <- gene_data %>%
    select(-matches("control"))
  
  # 提取组织名称
  tissues_pre <- str_extract(colnames(pre_data), "(?<=_)[^_\\s]+")
  tissues_post <- str_extract(colnames(post_data), "(?<=_)[^_\\s]+")
  tissues_pre <- sub("\\..*", "", tissues_pre)
  tissues_post <- sub("\\..*", "", tissues_post)
  # 构建结果数据框
  data <- data.frame(
    sample = c(colnames(pre_data), colnames(post_data)),
    tissue = c(tissues_pre, tissues_post),
    group = rep(c("Control", "aerobics"), times = c(length(pre_data), length(post_data))),
    expr = c(t(pre_data), t(post_data))
  )
  # 按中位数由高到低排列：
  data_new <- data %>% 
    group_by(tissue) %>% 
    mutate(median = median(expr), group_max = max(expr)) %>% 
    arrange(desc(median))
  
  # 调整因子顺序：
  data_new$tissue <- factor(data_new$tissue, levels = unique(data_new$tissue))
  
  data_new1 <- data_new %>% 
    group_by(tissue) %>% 
    mutate(group_number = length(unique(group)))
  
  # 此数据用于添加显著性检验的label：
  stat.test <- data_new1[which(data_new1$group_number == 2),] %>%
    group_by(tissue) %>%
    pairwise_t_test(
      expr ~ group, paired = FALSE, 
      p.adjust.method = "fdr"
    ) %>% 
    add_xy_position(x = "tissue")
  data_new1$group <- factor(data_new1$group, levels = c("Control", "aerobics"))
  # 重新画图：
  p <- ggplot()+
    # 基础图形：
    geom_boxplot(data = data_new1, aes(x = tissue, y =expr, fill = group),outlier.shape = 21, outlier.fill = "white")+
    # 设置颜色：
    scale_fill_manual(name = NULL, values = c("Control" = "#5488ef","aerobics" = "#d6503a" ))+
    # 设置主题：
    theme_bw()+
    labs(y = bquote(.(gene_name) ~ " Log-transformed Protein Abundance"))+
    #labs(y = TeX("mRNA Level ($log_{2}$CPM)", bold = T))+
    theme(axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.99, 0.99), 
          legend.justification = c(1,1))+
    # 显著性检验：
    stat_pvalue_manual(
      stat.test, label = "p.adj.signif", 
      bracket.size = 0.3, # 粗细
      tip.length = 0.01  # 两边竖线的长度
    )
  selected_columns <- stat.test[, c(1, 3, 4, 7, 8, 9, 10)]
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}
##########
panphosphoproteomic <- function(exp,gene_name) {
  # 提取基因表达数据
  gene_data <- exp %>% 
    filter(str_detect(row.names(exp), gene_name))
  # 提取控制组和有氧组数据
  Control <- gene_data[, grep(".Control", colnames(gene_data)), drop = FALSE]
  post_data <- gene_data %>%
    select(-matches("control"))
  # 构建结果数据框
  data <- data.frame(
    sample = c(colnames(Control), colnames(post_data)),
    phosphorylation = c(row.names(Control), row.names(post_data)),
    group = rep(c("Control", "aerobics"), times = c(length(Control), length(post_data))),
    expr = c(t(Control), t(post_data))
  )
  # 按中位数由高到低排列：
  data_new <- data %>% 
    group_by(phosphorylation) %>% 
    mutate(median = median(expr), group_max = max(expr)) %>% 
    arrange(desc(median))
  
  # 调整因子顺序：
  data_new$phosphorylation <- factor(data_new$phosphorylation, levels = unique(data_new$phosphorylation))
  
  data_new1 <- data_new %>% 
    group_by(phosphorylation) %>% 
    mutate(group_number = length(unique(group)))
  
  # 此数据用于添加显著性检验的label：
  stat.test <- data_new1[which(data_new1$group_number == 2),] %>%
    group_by(phosphorylation) %>%
    pairwise_t_test(
      expr ~ group, paired = FALSE, 
      p.adjust.method = "fdr"
    ) %>% 
    add_xy_position(x = "phosphorylation")
  data_new1$group <- factor(data_new1$group, levels = c("Control", "aerobics"))
  # 重新画图：
  p <- ggplot()+
    # 基础图形：
    geom_boxplot(data = data_new1, aes(x = phosphorylation, y =expr, fill = group),outlier.shape = 21, outlier.fill = "white")+
    # 设置颜色：
    scale_fill_manual(name = NULL, values = c("Control" = "#354898","aerobics" = "#972D36"))+
    # 设置主题：
    theme_bw()+
    #labs(y = bquote(.(gene_name)  " Log-transformed Protein Abundance"))
    labs(y =   bquote(.~" Log-Ratio of Normalized Phosphorylation in "~ .(gene_name)))+
    #labs(y = TeX("mRNA Level ($log_{2}$CPM)", bold = T))+
    theme(axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.99, 0.99), 
          legend.justification = c(1,1))+
    # 显著性检验：
    stat_pvalue_manual(
      stat.test, label = "p.adj.signif", 
      bracket.size = 0.3, # 粗细
      tip.length = 0.01  # 两边竖线的长度
    )
  selected_columns <- stat.test[, c(1, 3, 4, 7, 8, 9, 10)]
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
} 
####
GEOmouseFC<- function(genenames){
  data_matrix <- GEOFCmouse[genenames,c(1:8)]
  data_matrix <- na.omit(data_matrix)
  m<-ncol(data_matrix)
  #data_matrix <- t(data_matrix)
  #提取FC值与显著性
  padj_with_total <- grep("padj", names(data_matrix), value = TRUE)
  selected_padj <- data_matrix[, padj_with_total]
  FC_with_total <- grep("log2FoldChange", names(data_matrix), value = TRUE)
  selected_FC <- data_matrix[, FC_with_total]
  #转换长型矩阵
  data_padj <- gather(as.data.frame(t(selected_padj)),key = "DataSet",value = "padj")
  padj_underscore <- sub("_.*$", "", (colnames(selected_padj)))
  data_padj$dataSet <- rep(paste(padj_underscore))
  data_FC <- gather(as.data.frame(t(selected_FC)),key = "DataSet",value = "log2FoldChange")
  data_result <- cbind(data_padj,data_FC$log2FoldChange)
  colnames(data_result)[4] <- "log2FoldChange"
  #修改pos表示方法
  rows_count <- nrow(data_result)
  Neg_Pos <- rep("Neg",rows_count)
  Neg_Pos[which(data_result$log2FoldChange > 0)] <- "Pos"
  data_result$Neg_Pos <- Neg_Pos
  # 修改p值表示方法：
  pvalue <- rep(NA,rows_count)
  pvalue[which(data_result$padj > 0.05)] <- ">0.05"
  pvalue[which(data_result$padj < 0.05)] <- "<0.05"
  pvalue[which(data_result$padj < 0.01)] <- "<0.01"
  pvalue[which(data_result$padj < 0.001)] <- "<0.001"
  pvalue[which(data_result$padj < 0.0001)] <- "<0.0001"
  data_result$significance <- pvalue
  p<-ggplot(data_result,aes(dataSet,DataSet))+
    # 蓝色气泡图：
    geom_point(aes(fill=significance, size=abs(log2FoldChange)),color = "#999999",shape=21)+
    scale_fill_manual(values = c("#212c5f","#3366b1","#42b0e4","#7bc6ed","#dfe1e0"))+
    # 红色气泡图：
    geom_point(data = data_result[which(data_result$Neg_Pos == "Pos"),],
               aes(color=significance, size=abs(log2FoldChange)),shape=16)+
    scale_color_manual(values = c("#f26666","#f49699","#facccc","#facccc","#d9dbd9"))+
    # 主题：
    theme_bw()+
    theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
          # 坐标轴label方向：
          axis.text.x = element_text(angle = 45, hjust = 1),
          # 图例间距：
          legend.margin = margin(20,unit = 'pt'))+
    xlab("")+
    ylab("")+
    # 图例：
    guides(size = guide_legend(title = "log2FoldChange"),
           fill = guide_legend(title = expression("Adjusted P-values (padj)\nfor Downregulated Genes")),
           col = guide_legend(title = expression("Adjusted P-values (padj)\nfor upregulated Genes")))
  selected_columns <- data_result[, c(1, 2, 3, 4)]
  result <- list(plotData = p, tableData = selected_columns)
  return(result)
}
###
GEOdatamouse<- function(genename,Controlname,Experimentalname,Statistics){
  selected_gene <- GEOmouse[genename, ]
  #提取对照组和实验组的数据
  Controlcol <- c(Controlname)
  Control <-  selected_gene[, Controlcol]
  Experimentalcol <- c(Experimentalname)
  Experimental <- selected_gene[, Experimentalcol]
  Control$group <- "Control"
  Experimental$group <- "Experimental"
  df <- bind_rows(Control, Experimental)
  # 将数据框转换为长格式
  data_long <- df %>%
    pivot_longer(cols = -group, 
                 names_to = "variable", 
                 values_to = "value")
  data_long <- na.omit(data_long)
  # 添加IQR和中位数信息计算：
  data_long <- data_long %>%
    group_by(group) %>%
    mutate(
      IQR = IQR(value, na.rm = TRUE),   # 在计算IQR时忽略缺失值
      Med = median(value, na.rm = TRUE) # 在计算中位数时忽略缺失值
    ) %>%
    #mutate(IQR = IQR(value), Med = median(value)) %>%
    ungroup()
  data_long$group <- factor(data_long$group, levels = c( "Control","Experimental"))
  # 绘图：
  #####
  p <- ggplot(data_long, aes(group, value, fill = group)) +
    # 误差棒：
    geom_errorbar(aes(ymin = Med - IQR, ymax = Med + IQR),
                  width = 0.5, linetype = "dashed", position = position_dodge(width = 0.75)) +
    # 箱线图：
    geom_boxplot(color = NA, size = 1, notch = TRUE, position = position_dodge(width = 0.75)) +
    # 加入散点图：
    geom_jitter(aes(fill = group), 
                width = 0.2, # 控制散点的水平抖动
                size = 2.5, # 调整散点的大小以更好地看到边缘
                shape = 21, # 形状21（带边缘的点形状）
                color = "black", # 设置边缘颜色为黑色
                alpha = 0.7) + # 控制透明度
    # 颜色模式：
    scale_fill_manual(values = c("#4087c7", "#ec0f80")) +
    # 显著性检验：
    geom_signif(aes(group, value),
                # 指定你要比较的组，以list形式传入：
                comparisons = list(c("Control", "Experimental")),
                map_signif_level = TRUE, # 展示P值或显著性
                textsize = 3, # label大小
                test = Statistics, # 检验类型
                step_increase = 0.1 # 位置变化
    ) +
    xlab("") +
    ylab(bquote(.(genename)  ~ " mRNA Level (" * count * ")")) +
    # 主题调整：
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic"),
      legend.position = "none",
      axis.text.x = element_text(size = 12),    # 修改x轴刻度文字大小
      axis.text.y = element_text(size = 12),    # 修改y轴刻度文字大小
      axis.title.x = element_text(size = 14),   # 修改x轴标题文字大小
      axis.title.y = element_text(size = 14)    # 修改y轴标题文字大小
    )
  result <- list(plotData = p,table = data_long)
}

