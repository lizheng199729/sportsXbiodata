server <- function(input, output, session) { 
  observe({
    selected_tab <- input$tabs  # 获取当前选择的选项卡
    
    # 如果没有选择选项卡，默认为"01"
    if (is.null(selected_tab) || selected_tab == "") {
      selected_tab <- "01"
    }
    
    # 根据选中的tab更新控件
    if (selected_tab == "01") {
        updateSelectizeInput(session, "selectGene01", 
                             choices = row.names(HumanDifferentialGenes), 
                             server = TRUE)
      
    } else if (selected_tab == "02") {
      updateSelectizeInput(session, "selectGene02", 
                           choices = row.names(difexpdata), 
                           server = TRUE)
    } else if (selected_tab == "03") {
      updateSelectizeInput(session, "selectkeywords03", 
                           choices = keywords$Entry, 
                           server = TRUE)
      
    }  else if (selected_tab == "two") {
        updateSelectizeInput(session, "selectGeneTwo", 
                           choices = row.names(tissueexpdata), 
                           server = TRUE)
      
    } else if (selected_tab == "four") {
      updateSelectizeInput(session, "selectGenefour", 
                           choices = row.names(pyGENEdata), 
                           server = TRUE)
      
    } else if (selected_tab == "six") {
      updateSelectizeInput(session, "selectGene6", 
                           choices = row.names(genedata), 
                           server = TRUE)
      
    } else if (selected_tab == "eight") {
      updateSelectizeInput(session, "selectGeneeight", 
                           choices = row.names(Proteindata), 
                           server = TRUE)
      
    } else if (selected_tab == "nine") {
      updateSelectizeInput(session, "selectGenenine", 
                           choices = row.names(Proteindata), 
                           server = TRUE)
    }
    else if (selected_tab == "c") {
      updateSelectizeInput(session, "selectDataset_c", 
                           choices = c("White.Adipose", "Lung", "Liver", "Kidney", 
                                       "Heart", "Gastrocnemius", "cortex"), 
                           server = TRUE)
     }
    })
 
  observeEvent(input$btnSubmit01, {
    withProgress(message = 'Processing...', value = 0, {
      req(input$selectGene01)
      Humangeo_result <- HumanGEO(input$selectGene01)
      output$plot01_1 <- renderPlot({
        Humangeo_result$plotData
      })
      output$geneDataTable01_1 <- renderDataTable({ Humangeo_result$tableData })
      output$downloadResult01_1 <- downloadHandler(
        filename = function() { paste0(input$selectGene01, "_GEO_plot.png") },
        content = function(file) { 
          ggsave(file, Humangeo_result$plotData, width = 14, height = 6, dpi = 300)  # 宽10英寸，高8英寸，分辨率300
        }
      )
      output$downloadTable01_1 <- downloadHandler(
        filename = function() { paste0(input$selectGene01, "_GEO_table.csv") },
        content = function(file) { write.csv(geo_result$tableData, file) }
      )
      incProgress(1, detail = "Completed.")
    })
  })
  observeEvent(input$btnSubmit02, {
    withProgress(message = 'Processing...', value = 0, {
      req(input$selectGene02)
      motrpac_result <- pansport(input$selectGene02)
      output$plot02_1 <- renderPlot({ motrpac_result$plotData })
      output$geneDataTable02_1 <- renderDataTable({ motrpac_result$tableData })
      output$downloadResult02_1 <- downloadHandler(
        filename = function() { paste0(input$selectGene02, "_MoTrPAC_plot.png") },
        content = function(file) { ggsave(file, motrpac_result$plotData) }
      )
      output$downloadTable02_1 <- downloadHandler(
        filename = function() { paste0(input$selectGene02, "_MoTrPAC_table.csv") },
        content = function(file) { write.csv(motrpac_result$tableData, file) }
      )
      
      geo_result <- mouseGEO(input$selectGene02)
      output$plot02_2 <- renderPlot({
        geo_result$plotData
      })
      output$geneDataTable02_2 <- renderDataTable({ geo_result$tableData })
      output$downloadResult02_2 <- downloadHandler(
        filename = function() { paste0(input$selectGene02, "_GEO_plot.png") },
        content = function(file) { 
          ggsave(file, geo_result$plotData, width = 10, height = 10, dpi = 300)  # 宽10英寸，高8英寸，分辨率300
        }
      )
      output$downloadTable02_2 <- downloadHandler(
        filename = function() { paste0(input$selectGene02, "_GEO_table.csv") },
        content = function(file) { write.csv(geo_result$tableData, file) }
      )
      incProgress(1, detail = "Completed.")
    })
  })
  observeEvent(input$btnSubmit03, {
    withProgress(message = 'Processing...', value = 0, {
      req(input$selectkeywords03)
      result03 <- graph(input$selectkeywords03)
      output$plot03_1 <- renderPlot({ result03$plotData })
      output$geneDataTable03_1 <- renderDataTable({ result03$tableData })
      output$downloadResult03_1 <- downloadHandler(
        filename = function() { paste0(input$selectkeywords03, "_graph_plot.png") },
        content = function(file) { ggsave(file, result03$plotData) }
      )
      output$downloadTable03_1 <- downloadHandler(
        filename = function() { paste0(input$selectGene03, "_graph_table.csv") },
        content = function(file) { write.csv(motrpac_result$tableData, file) }
      )
      incProgress(1, detail = "Completed.")
    })
  })
  # 基于关键词从 df_split 中筛选
  filtered_data <- reactive({
    keyword <- trimws(input$keyword)
    if (keyword == "") {
      df_split
    } else {
      pattern <- paste0("\\b", keyword, "\\b")  # 精确匹配完整基因名
      df_split[grepl(pattern, df_split$Gene_Link, ignore.case = TRUE), ]
    }
  })
  
  # 展示筛选后的数据表格
  output$gene_table <- renderDT({
    datatable(
      filtered_data()[, c("Gene_Link", "Verification",  
                          "Forms of exercise", "Species", "Sex", "PMID_Link")],
      escape = FALSE,
      rownames = FALSE,
      colnames = c("Gene (GeneCards)", "Verification", 
                   "Forms of Exercise", "Species", "Sex", "PMID (PubMed)")
    )
  })
  observeEvent(input$btnSubmit2, {
    withProgress(message = 'Processing...', value = 0, {
      req(input$selectGeneTwo)  # 确保基因名已选择
      req(input$selectTissue)
      selected_data2 <- tissueexp(input$selectGeneTwo, input$selectTissue,input$Statistical2)
      
      # 假设您有一个 plotOutput 和一个 dataTableOutput
      output$plot2 <- renderPlot({
        plot(selected_data2$plotData)  # 生成并渲染图表
      },width = input$textWidthtwo, height = input$textHeighttwo)
      output$downloadResultTwo <- downloadHandler(
        filename = function() {
          paste("data-plot", Sys.Date(), ".png", sep="")
        },
        content = function(file) {
          # 假设图像已经存储在一个名为 plotVar 的 ggplot 对象中
          width_in_inches <- as.numeric(input$textWidthtwo) / 96  # 将像素转换为英寸
          height_in_inches <- as.numeric(input$textHeighttwo) / 96  # 将像素转换为英寸
          ggsave(file, plot = selected_data2$plotData, width = width_in_inches, height = height_in_inches, dpi = 300)
        })
      incProgress(1, detail = "Completed.")
    })
  })
  observeEvent(input$btnSubmit4, {
    withProgress(message = 'Processing...', value = 0, { 
      req(input$selectGenefour)  # 确保基因名已选择
      selected_data <- Phenotypic_correlation(input$selectGenefour)
      
      # 假设您有一个 plotOutput 和一个 dataTableOutput
      output$plot4 <- renderPlot({
        plot(selected_data$plotData)  # 生成并渲染图表
      },width = input$textWidthfour, height = input$textHeightfour)
      
      output$geneDataTable4 <- renderDataTable({
        datatable(selected_data$tableData)  # 生成并渲染数据表
      })
      output$downloadResultfour <- downloadHandler(
        filename = function() {
          paste("data-plot", Sys.Date(), ".png", sep="")
        },
        content = function(file) {
          # 假设图像已经存储在一个名为 plotVar 的 ggplot 对象中
          width_in_inches <- as.numeric(input$textWidthfour) / 96  # 将像素转换为英寸
          height_in_inches <- as.numeric(input$textHeightfour) / 96  # 将像素转换为英寸
          ggsave(file, plot = selected_data$plotData, width = width_in_inches, height = height_in_inches, dpi = 300)
        })
      output$downloadTablefour <- downloadHandler(
        filename = function() {
          paste("data-table", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          # 假设数据表已经存储在一个名为 dataTable 的变量中
          write.csv(selected_data$tableData, file)
        })
      incProgress(1, detail = "Completed.")
    })
  })
  observeEvent(input$btnSubmit6, {
    withProgress(message = 'Processing...', value = 0, { 
      req(input$selectGene6)  # 确保基因名已选择
      validate(
        need(length(input$selectGene6) > 1, "Please select more than one gene.")
      )
      immgeneresult <- immgene(input$selectGene6)  # 假设 immgene 函数返回含有 plotData 和 tableData 的列表
      # 分别渲染图表和数据表
      output$plot6 <- renderPlot({
        plot(immgeneresult$plotData)
      },width = input$textWidthsix, height = input$textHeightsix)
      
      output$geneDataTable6 <- renderDataTable({
        datatable(immgeneresult$tableData)  # 使用结果中的 tableData 生成数据表
      })
      output$downloadResultsix <- downloadHandler(
        filename = function() {
          paste("data-plot", Sys.Date(), ".png", sep="")
        },
        content = function(file) {
          width_in_inches <- as.numeric(input$textWidthsix) / 96  # 将像素转换为英寸
          height_in_inches <- as.numeric(input$textHeightsix) / 96  # 将像素转换为英寸
          # 假设图像已经存储在一个名为 plotVar 的 ggplot 对象中
          ggsave(file, plot = immgeneresult$plotData, width = width_in_inches, height = height_in_inches, dpi = 300)
        })
      output$downloadTablesix <- downloadHandler(
        filename = function() {
          paste("data-table", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          # 假设数据表已经存储在一个名为 dataTable 的变量中
          write.csv(immgeneresult$tableData, file)
        })
      incProgress(1, detail = "Completed.")
    })
  })
  values_c <- reactiveValues()
  # 4. 点击 Step 1 后筛选对应数据
  # 4. 点击 Step 1 后筛选对应数据
  observeEvent(input$btnSubmit_c, {
    req(input$selectDataset_c)
    values_c$selected_data_c <- group[grep(input$selectDataset_c, group$group), ]
  })
  
  # 5. 根据选择器更新 gene 选项
  observe({
    req(values_c$selected_data_c)
    gene_choices <- values_c$selected_data_c  # 假设这里有 gene 列
    updateSelectizeInput(session, "selectDataset_c.1", choices = gene_choices, server = TRUE)
  })
  
  observe({
    req(values_c$selected_data_c)
    selected_gene1 <- input$selectDataset_c.1
    gene_choices <- values_c$selected_data_c
    filtered_choices <- setdiff(gene_choices, selected_gene1)
    updateSelectizeInput(session, "selectDataset_c.2", choices = filtered_choices, server = TRUE)
  })
  observeEvent(input$btnSubmit_c.1, {
    withProgress(message = 'Processing...', value = 0, {
      req(input$selectDataset_c.1)
      req(input$selectDataset_c.2)
      req(input$selectDataset_c)
      req(input$Threshold_c)
      req(input$textHeight_c)
      req(input$textWidth_c)
      
      # 检查 selectDataset_c.1 和 selectDataset_c.2 中是否包含 input$selectDataset_c
      organ_keyword <- input$selectDataset_c
      if (!(grepl(organ_keyword, input$selectDataset_c.1) && grepl(organ_keyword, input$selectDataset_c.2))) {
        showNotification(paste0("❌ Error: One or both selected samples do not contain the organ keyword '", organ_keyword, "'."), type = "error", duration = NULL)
        return(NULL)
      }
      
      RNA_seq <- switch(organ_keyword,
                        "White.Adipose" = RNA_White,
                        "Lung" = RNA_Lung,
                        "Liver" = RNA_Liver,
                        "Kidney" = RNA_Kidney,
                        "Heart" = RNA_Heart,
                        "Gastrocnemius" = RNA_Gastrocnemius,
                        "cortex" = RNA_Cortex,
                        stop("❌ Unknown organ"))
      proteomics <- switch(organ_keyword,
                           "White.Adipose" = proteomics_white,
                           "Lung" = proteomics_lung,
                           "Liver" = proteomics_liver,
                           "Kidney" = proteomics_kidney,
                           "Heart" = proteomics_heart,
                           "Gastrocnemius" = proteomics_gastrocnemius,
                           "cortex" = proteomics_cortex,
                           stop("❌ Unknown organ"))
      
      result <- Multi_omics(
        Control = input$selectDataset_c.1,
        Experimental = input$selectDataset_c.2,
        Threshold = input$Threshold_c,
        RNA = RNA_seq,
        protein = proteomics
      )
      validate(
        need(!is.null(result), "❌ Analysis failed: Multi_omics returned NULL."),
        need(!is.null(result$plotData), "❌ Analysis failed: plotData is missing."),
        need(!is.null(result$table), "❌ Analysis failed: table is missing.")
      )
      
      output$plot_c <- renderPlot({ result$plotData },
                                  width = input$textWidth_c,
                                  height = input$textHeight_c)
      output$Table_c <- renderDataTable({
        DT::datatable(result$table)
      })
      output$downloadResult_c <- downloadHandler(
        filename = function() { paste0("multi_omics_", Sys.Date(), ".png") },
        content = function(file) {
          ggsave(file, plot = result$plotData,
                 width = input$textWidth_c / 96,
                 height = input$textHeight_c / 96,
                 dpi = 300)
        }
      )
      output$downloadTable_c <- downloadHandler(
        filename = function() { paste0("multi_omics_", Sys.Date(), ".csv") },
        content = function(file) {
          write.csv(result$table, file, row.names = TRUE)
        }
      )
      
      incProgress(1)
    })
  })
  observeEvent(input$btnSubmit8, {
    req(input$selectGeneeight)  # 确保已选择蛋白名称
    withProgress(message = 'Processing...', value = 0, {
      # 执行数据处理
      selected_data <- panprotein(input$selectGeneeight)
      
      # 更新进度条
      incProgress(0.5, detail = "Rendering plot and table...")
      
      # 渲染图表
      output$plot8 <- renderPlot({
        plot(selected_data$plotData)  # 假设 plotData 是一个图表对象
      },width = input$textWidtheight, height = input$textHeighteight)
      
      # 渲染数据表
      output$proteinDataTableeight <- renderDataTable({
        datatable(selected_data$tableData)  # 假设 tableData 是一个数据表对象
      })
      
      # 下载图像处理
      output$downloadResulteight <- downloadHandler(
        filename = function() {
          paste("data-plot", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
          width_in_inches <- as.numeric(input$textWidtheight) / 96  # 将像素转换为英寸
          height_in_inches <- as.numeric(input$textHeighteight) / 96  # 将像素转换为英寸
          ggsave(file, plot = selected_data$plotData, width = width_in_inches, height = height_in_inches, dpi = 300)
        }
      )
      
      # 下载表格处理
      output$downloadTableeight <- downloadHandler(
        filename = function() {
          paste("data-table", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(selected_data$tableData, file)
        }
      )
      
      # 更新进度条到完成
      incProgress(1, detail = "Completed.")
    })
  })
  observeEvent(input$btnSubmit9 , {
    # 根据选择的数据集名称加载数据
    selected_dataset <- switch(input$selectDatasetnine,
                               "cortex" = cortex,  # dataset1, dataset2等是预先加载的数据框
                               "gastrocnemius" = gastrocnemius,
                               "heart" = heart,
                               "kidney" = kidney,
                               "liver" = liver,
                               "lung" = lung,
                               "white-adipose" = white)
    withProgress(message = 'Processing...', value = 0, {
      # 提取并处理基因名
      name <- sub("_.*$", "", row.names(selected_dataset))
      name <- unique(name)
      req(input$selectGenenine)
      selected_gene <- input$selectGenenine
      
      # 检查选中的基因是否存在于数据集中
      if (!selected_gene %in% name) {
        # 如果基因不在数据集中，显示错误通知
        showNotification("Error: The protein is not found in the dataset.", type = "error")
        return()  # 终止后续操作
      }
      
      result9 <- panphosphoproteomic(selected_dataset,input$selectGenenine)
      output$plot9 <- renderPlot({
        plot(result9$plotData)  # 生成并渲染图表
      },width = input$textWidthnine, height = input$textHeightnine)
      output$geneDataTable9 <- renderDataTable({
        datatable(result9$tableData)  # 生成并渲染数据表
      })
      output$downloadResultnine <- downloadHandler(
        filename = function() {
          paste("data-plot", Sys.Date(), ".png", sep="")
        },
        content = function(file) {
          width_in_inches <- as.numeric(input$textWidthnine) / 96  # 将像素转换为英寸
          height_in_inches <- as.numeric(input$textHeightnine) / 96  # 将像素转换为英寸
          # 假设图像已经存储在一个名为 plotVar 的 ggplot 对象中
          ggsave(file, plot = result9$plotData, width = width_in_inches, height = height_in_inches, dpi = 300)
        })
      output$downloadTablenine <- downloadHandler(
        filename = function() {
          paste("data-table", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          # 假设数据表已经存储在一个名为 dataTable 的变量中
          write.csv(result9$tableData, file)
        })
      incProgress(1, detail = "Completed.") 
    })
  })
###########
  file_data_a <- reactive({
    req(input$file_upload_a)
    
    df <- tryCatch(
      read.csv(input$file_upload_a$datapath, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE),
      error = function(e) stop("❌ Error reading the file: ", e$message)
    )
    
    validate(
      need(ncol(df) == 2, "❌ Error: The file must contain exactly 2 columns."),
      need(identical(colnames(df), c("gene", "log2FC")),
           paste0("❌ Column names must be exactly: 'gene' and 'log2FC'. Found: ", paste(colnames(df), collapse = ", "))),
      need(is.numeric(df$log2FC), "❌ Error: 'log2FC' column must be numeric."),
      need(all(is.finite(df$log2FC)), "❌ Error: 'log2FC' must not contain NA, NaN, or Inf.")
    )
    
    return(df)
  })
  
  output$file_check_a <- renderPrint({
    tryCatch({
      df <- file_data_a()
      cat("✅ Success: File format is valid.\n\nPreview:\n")
      print(head(df, 5))
    }, error = function(e) {
      cat("❌", e$message)
    })
  })
  
  updateSelectizeInput(session, "selectterm_a", 
                       choices = terma$term, 
                       server = TRUE)
  
  observeEvent(input$btnSubmit_a, {
    withProgress(message = 'Processing...', value = 0, {
      df <- file_data_a()
      # 处理物种选择
      term2gene <- tryCatch({
        switch(input$selectDataset_a,
               "human" = human,
               "mouse" = mouse,
               stop("❌ Unknown species selection."))
      }, error = function(e) {
        showNotification(e$message, type = "error")
        return(NULL)
      })
      
      # 如果 term2gene 出错，则中止
      req(!is.null(term2gene))
      
      # 检查所选 term 是否在 term2gene
      if (!input$selectterm_a %in% term2gene$term) {
        showNotification("❌ Selected term is not found in the species data. Please check.", type = "error")
        return()
      }
      
      # 调用 sportgenerich，捕获可能出错
      result <- tryCatch({
        sportgenerich(df, term2gene, input$selectterm_a)
      }, error = function(e) {
        showNotification(paste("❌ GSEA calculation failed:", e$message), type = "error")
        return(NULL)
      })
      
      req(!is.null(result))
      
      # 正常绘图和表格输出
      output$plot_a <- renderPlot({ result$plotData },
                                  width = input$textWidth_a,
                                  height = input$textHeight_a)
      
      output$Table_a <- renderDataTable({
        DT::datatable(result$tableData)
      })
      
      output$downloadResult_a <- downloadHandler(
        filename = function() { paste0("GSEA_plot_", Sys.Date(), ".png") },
        content = function(file) {
          ggsave(file, plot = result$plotData,
                 width = input$textWidth_a / 96,
                 height = input$textHeight_a / 96,
                 dpi = 300)
        }
      )
      
      output$downloadTable_a <- downloadHandler(
        filename = function() { paste0("GSEA_table_", Sys.Date(), ".csv") },
        content = function(file) {
          write.csv(result$tableData, file, row.names = FALSE)
        }
      )
      
      incProgress(1)
    })
  })
########
  # 1. 上传并读取用户 CSV 文件
  file_data <- reactive({
    user_data <- tryCatch(
      read.csv(input$file_upload_b$datapath, 
               row.names = 1, 
               fileEncoding = "UTF-8-BOM", 
               stringsAsFactors = FALSE),
      error = function(e) stop("❌ Error reading the file: ", e$message)
    )
    
    # 校验格式
    required_cols <- c("log2FoldChange", "padj")
    validate(
      need(all(required_cols %in% colnames(user_data)),
           paste0("❌ The file must contain: ", paste(required_cols, collapse = ", "))),
      need(is.numeric(user_data$log2FoldChange), "❌ log2FoldChange must be numeric"),
      need(all(is.finite(user_data$log2FoldChange)), "❌ log2FoldChange has NA/NaN/Inf"),
      need(!any(duplicated(rownames(user_data))), "❌ gene (row names) must be unique")
    )
    
    return(user_data)
  })
  
  # 2. 显示文件上传检查
  output$file_check <- renderPrint({
    tryCatch({
      user_data <- file_data()
      cat("✅ Success: File format is valid.\n\nPreview:\n")
      print(head(user_data, 5))
    }, error = function(e) {
      cat("❌", e$message)
    })
  })
  
  # 3. 初始化器官选择列表
  values_b <- reactiveValues()
  observe({
    updateSelectizeInput(session, "selectDataset_b", 
                         choices = c("White.Adipose", "Lung", "Liver", "Kidney", 
                                     "Heart", "Gastrocnemius", "cortex"), 
                         server = TRUE)
  })
  
  # 4. 点击 Step 1 后筛选对应数据
  observeEvent(input$btnSubmit_b, {
    req(input$selectDataset_b)
    values_b$selected_data_b <- group[grep(input$selectDataset_b, group$group), ]
  })
  
  # 5. 根据选择器更新 gene 选项
  observe({
    req(values_b$selected_data_b)
    gene_choices <- values_b$selected_data_b  # 假设这里有 gene 列
    updateSelectizeInput(session, "selectDataset_b.1", choices = gene_choices, server = TRUE)
  })
  
  observe({
    req(values_b$selected_data_b)
    selected_gene1 <- input$selectDataset_b.1
    gene_choices <- values_b$selected_data_b
    filtered_choices <- setdiff(gene_choices, selected_gene1)
    updateSelectizeInput(session, "selectDataset_b.2", choices = filtered_choices, server = TRUE)
  })
  
  # 6. Step 2 分析逻辑
  observeEvent(input$btnSubmit_b.1, {
    withProgress(message = 'Processing...', value = 0, {
      req(input$selectDataset_b.1)
      req(input$selectDataset_b.2)
      req(input$selectDataset_b)
      
      RNA_seq <- switch(input$selectDataset_b,
                        "White.Adipose" = RNA_White,
                        "Lung" = RNA_Lung,
                        "Liver" = RNA_Liver,
                        "Kidney" = RNA_Kidney,
                        "Heart" = RNA_Heart,
                        "Gastrocnemius" = RNA_Gastrocnemius,
                        "cortex" = RNA_Cortex,
                        stop("❌ Unknown organ"))
      
      proteomics <- switch(input$selectDataset_b,
                           "White.Adipose" = proteomics_white,
                           "Lung" = proteomics_lung,
                           "Liver" = proteomics_liver,
                           "Kidney" = proteomics_kidney,
                           "Heart" = proteomics_heart,
                           "Gastrocnemius" = proteomics_gastrocnemius,
                           "cortex" = proteomics_cortex,
                           stop("❌ Unknown organ"))
      
      # 判断是否上传了用户数据
      has_upload <- !is.null(input$file_upload_b)
      
      result <- tryCatch({
        if (has_upload) {
          user_rich(
            Control = input$selectDataset_b.1,
            Experimental = input$selectDataset_b.2,
            RNA_upload = file_data(),
            RNA = RNA_seq,
            protein = proteomics
          )
        } else {
          user_rich(
            Control = input$selectDataset_b.1,
            Experimental = input$selectDataset_b.2,
            RNA = RNA_seq,
            protein = proteomics
          )
        }
      }, error = function(e) {
        showNotification(paste("❌ Error:", e$message), type = "error")
        return(NULL)
      })
      
      req(result)
      
      output$plot_b <- renderPlot({ result$plotData },
                                  width = input$textWidth_b,
                                  height = input$textHeight_b)
      
      output$Table_b <- renderDataTable({
        DT::datatable(result$table)
      })
      
      output$downloadResult_b <- downloadHandler(
        filename = function() { paste0("multi_omics_", Sys.Date(), ".png") },
        content = function(file) {
          ggsave(file, plot = result$plotData,
                 width = input$textWidth_b / 96,
                 height = input$textHeight_b / 96,
                 dpi = 300)
        }
      )
      
      output$downloadTable_b <- downloadHandler(
        filename = function() { paste0("multi_omics_", Sys.Date(), ".csv") },
        content = function(file) {
          write.csv(result$table, file, row.names = TRUE)
        }
      )
      
      incProgress(1)
    })
  })
  output$"eleven-01" <- downloadHandler(
    filename = function() {
      "Rapid search of human genes.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/Rapid search of human genes.pdf", file)
    }
  )
  output$"eleven-02" <- downloadHandler(
    filename = function() {
      "Rapid search of mouse genes.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/Rapid search of mouse genes.pdf", file)
    }
  )
  output$"eleven-03" <- downloadHandler(
    filename = function() {
      "Knowledge graph.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/Knowledge graph.pdf", file)
    }
  )
  output$"eleven-two" <- downloadHandler(
    filename = function() {
      "Exercise-Induced Gene Changes at Different Duration.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/Exercise-Induced Gene Changes at Different Duration.pdf", file)
    }
  )
  output$"eleven-four" <- downloadHandler(
    filename = function() {
      "Phenotype and Gene.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/Phenotype and Gene.pdf", file)
    }
  )
  output$"eleven-six" <- downloadHandler(
    filename = function() {
      "immune cell Gene Related.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/immune cell Gene Related.pdf", file)
    }
  )
  output$"eleven-eight" <- downloadHandler(
    filename = function() {
      "Proteomics.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 假设你的文件已经存在于某个路径下
      file.copy("/srv/shiny-server/myapp/www/Proteomics.pdf", file)
    }
  )
  output$"eleven-eight" <- downloadHandler(
    filename = function() {
      "Proteomics.pdf"  # 你想要下载的文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/Proteomics.pdf", file)
    }
  )
  output$"eleven-nine" <- downloadHandler(
    filename = function() {
      "Phosphoproteomic.pdf"  # 下载的文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/Phosphoproteomic.pdf", file)
    }
  )
  output$"eleven-c" <- downloadHandler(
    filename = function() {
      "Multi-omics.pdf"  # 下载的文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/Multi-omics.pdf", file)
    }
  )
  output$"eleven-a" <- downloadHandler(
    filename = function() {
      "SportEnrich.pdf"  # 文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/SportEnrich.pdf", file)
    }
  )
  output$"eleven-a_1" <- downloadHandler(
    filename = function() {
      "SportEnrich-GSEA Example File.csv"  # 文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/SportEnrich-GSEA test.csv", file)
    }
  )
  output$"eleven-b" <- downloadHandler(
    filename = function() {
      "Exercise-Drug Convergence.pdf"  # 文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/Exercise-Drug Convergence.pdf", file)
    }
  )
  output$"eleven-b_1" <- downloadHandler(
    filename = function() {
      "Exercise-Drug Convergence Example File.csv"  # 你想要下载的文件名
    },
    content = function(file) {
      # 路径下
      file.copy("/srv/shiny-server/myapp/www/Exercise-Drug Convergence test.csv", file)
    }
  )
}
