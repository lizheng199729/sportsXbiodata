ui <- dashboardPage(
  dashboardHeader(title = "sportXbiodata"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "zero", icon = icon("home")),
      menuItem("Rapid search", icon = icon("search"),
               # 子选项
               menuSubItem("Rapid search of human genes", tabName = "01"),
               menuSubItem("Rapid search of mouse genes", tabName = "02"),
               menuSubItem("Knowledge graph", tabName = "03"),
               menuSubItem("Paper-GENE", tabName = "04")
      ),
      menuItem("MoTrPAC Database", icon = icon("dna"),
               # 子选项
               menuSubItem("Gene Changes at Different Duration", tabName = "two"),
               menuSubItem("Phenotype and Gene", tabName = "four"),
               menuSubItem("immune cell Gene Related", tabName = "six"),
               menuSubItem("proteomics", tabName = "eight"),
               menuSubItem("phosphoproteomic", tabName = "nine"),
               menuSubItem("Multi-omics", tabName = "c")),
      menuItem("Data Upload", icon = icon("file-upload"),
               menuSubItem("SportEnrich-GSEA", tabName = "a"),
               menuSubItem("Exercise-Drug Convergence", tabName = "b"))
    )
  ),
  dashboardBody(
    tags$style(HTML("
      .custom-text {
        font-size: 30px;
        color: #3FA0C0;
        text-align: center;
      }
      .text2 {
        font-size: 20px;
        text-align: center;
      }
      .custom-box {
        background-color: #f5f5f5;
        border: none !important;
      }
      .text3 {
        font-size: 20px;
      }
    ")),
    tabItems(
      tabItem(tabName = "zero",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  class = "custom-box",
                  p(class = "custom-text", "sportXbiodata"),
                  p(class = "text2", "sportXbiodata is a comprehensive, manually curated multi-omics database dedicated 
                    to exercise biology. It contains 2,180 sequencing datasets and 3,600 publications covering exercise, metabolism, 
                    nutrition, and disease, with data from sources including NCBI, MoTrPAC, GEO, PubMed, 
                    and GeneCard. The database integrates multiple omics layers—RNA-seq, proteomics, 
                    phosphoproteomics—and offers information on gene expression, differential expression and annotation, 
                    integrated gene-protein profiles, phenotype-gene correlations, immune associations, and knowledge graphs. 
                    It features visualization options ready for publication, supports user-uploaded 
                    data for exercise-related gene enrichment analysis, and includes an exercise-drug homology module to 
                    identify key targets linking exercise and drug interventions. sportsXbiodata provides 
                    a powerful resource for advancing research on the molecular mechanisms of exercise and exercise-related diseases.")
                )
              ),
                fluidRow(
                box(
                  width = 12,
                  title = "Description of the sample",
                  status = "primary",
                  solidHeader = TRUE,
                  tags$img(src = "幻灯片1(1).png", 
                           alt = "Description of the sample1", 
                           width = "1200px")
                )
            ),
                fluidRow(
                box(
                  width = 12,
                  title = "Main functions",
                  status = "primary",
                  solidHeader = TRUE,
                  tags$img(src = "幻灯片2.PNG", 
                           alt = "Main functions1", 
                           width = "1200px")
                )
              ),
            fluidRow(
              box(
                width = 6,
                title = "Data processing",
                status = "primary",
                solidHeader = TRUE,
                tags$img(src = "幻灯片3.PNG", 
                         alt = "Data processing1", 
                         width = "600px")
              ),
              box(
                width = 6,
                title = "Website Structure",
                status = "primary",
                solidHeader = TRUE,
                tags$img(src = "sportXbiodata(1).png", 
                         alt = "Website Structure1", 
                         width = "600px")
              )
            )
          
      ),
      tabItem(tabName = "01",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "Human Gene Search"),
                    p(class = "text2", "Conveniently querying gene expression differences 
                      across various exercise-related tissues and experimental conditions.")
                )
              ),
              fluidRow(
                box(title = "Gene Input", status = "warning", solidHeader = TRUE, collapsible = TRUE,width = 12,
                    selectInput("selectGene01", "GENENAME", choices = NULL, multiple = T),
                    actionButton("btnSubmit01", "Submit", class = "btn-primary"),
                    downloadButton("eleven-01", "instruction manual")
                )
              ),
              fluidRow(
                box(title = "GEO Result", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    plotOutput("plot01_1"),
                    downloadButton("downloadResult01_1", "Download Plot")
                ),
                box(title = "GEO Data Table", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    dataTableOutput("geneDataTable01_1"),
                    downloadButton("downloadTable01_1", "Download Table")
                )
              ) 
      ),
      tabItem(tabName = "02",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "Mouse Gene Search"),
                    p(class = "text2", "Conveniently querying gene expression differences 
                      across various exercise-related tissues and experimental conditions.")
                )
              ),
              fluidRow(
                box(title = "Gene Input", status = "warning", solidHeader = TRUE, collapsible = TRUE,width = 12,
                    selectInput("selectGene02", "GENENAME:", choices = NULL, multiple = FALSE),
                    actionButton("btnSubmit02", "Submit", class = "btn-primary"),
                    downloadButton("eleven-02", "instruction manual")
                )
              ),
              fluidRow(
                box(title = "MoTrPAC Result", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    plotOutput("plot02_1"),
                    downloadButton("downloadResult02_1", "Download Plot")
                ),
                box(title = "MoTrPAC Data Table", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    dataTableOutput("geneDataTable02_1"),
                    downloadButton("downloadTable02_1", "Download Table")
                )
              ),
              fluidRow(
                box(title = "GEO Result", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    plotOutput("plot02_2", height = "500px", width = "100%"),
                    downloadButton("downloadResult02_2", "Download Plot")
                ),
                box(title = "GEO Data Table", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    dataTableOutput("geneDataTable02_2"),
                    downloadButton("downloadTable02_2", "Download Table")
                )
              )
      ),
      tabItem(tabName = "03",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "Knowledge graph"),
                    p(class = "text2", "By compiling a comprehensive collection of exercise-related literature and entering specific keywords, an exercise-gene 
                      knowledge graph can be automatically generated.")
                )
              ),
              fluidRow(
                box(title = "Enter keywords", status = "warning", solidHeader = TRUE, collapsible = TRUE,width = 12,
                    selectInput("selectkeywords03", "keywords:", choices = NULL, multiple = FALSE),
                    actionButton("btnSubmit03", "Submit", class = "btn-primary"),
                    downloadButton("eleven-03", "instruction manual")
                )
              ),
              fluidRow(
                box(title = "Knowledge graph Result", solidHeader = TRUE, collapsible = TRUE, width = 6,
                    plotOutput("plot03_1"),
                    downloadButton("downloadResult03_1", "Download Plot")
                ),
                box(title = "Knowledge graph Table", solidHeader = TRUE, collapsible = TRUE, width = 6,
                    dataTableOutput("geneDataTable03_1"),
                    downloadButton("downloadTable03_1", "Download Table")
                )
              ) 
      ),
      tabItem(tabName = "04",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "Paper-GENE"),
                    p(class = "text2", "Paper-GENE is a literature query module within the database, 
                      encompassing 3,600 publications and enabling rapid association of genes with research topics.")
                )
              ),
              fluidRow(
                box(width = 3, solidHeader = TRUE, status = "primary", title = "Filter Genes",
                    textInput("keyword", "Search gene name", "")
                ),
                box(width = 9, title = "Results Table", status = "success", solidHeader = TRUE,
                    DTOutput("gene_table"))
              )
      ),
      tabItem(tabName = "two",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  class = "custom-box",
                  p(class = "custom-text", "Gene Changes at Different Duration"),
                  p(class = "text2", "Investigate the dynamic changes of specific genes across various tissues at different exercise time points.e")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "GeneNAME input", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectGeneTwo", "GENENAME:",choices = NULL, multiple = FALSE),
                  selectInput("selectTissue", "Tissue:",
                              choices = c('blood', 'hippocampus', 'cortex', 'hypothalamus',
                                          'gastrocnemius', 'vastus', 'heart', 'kidney', 'adrenal', 'colon',
                                          'spleen', 'testes', 'ovaries', 'lung', 'liver', 'brown', 'white', 'vena')),
                  selectInput("Statistical2","Statistical method",
                              choices = c("wilcox.test","t.test")),
                  numericInput("textHeighttwo", "height:", 300),
                  numericInput("textWidthtwo", "width:", 300),
                  actionButton("btnSubmit2", "Submit", class = "btn-primary"),
                  downloadButton("eleven-two", "instruction manual")
                )),
              fluidRow(
                box(
                  title = "Result", solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  plotOutput("plot2", height = "100%", width = "100%"),
                  downloadButton("downloadResultTwo", "Download Results")
                )
          )
      ),
      tabItem(tabName = "four",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  class = "custom-box",
                  p(class = "custom-text", "Phenotype and Gene"),
                  p(class = "text2", "Evaluate the correlation between the expression of 
                    target genes across different organs and phenotypic data, which may indicate potential associations between gene expression and final 
                    phenotypic traits such as body weight.")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "GeneNAME input", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectGenefour", "GENENAME:",choices = NULL, multiple = FALSE),
                  numericInput("textHeightfour", "height:", 500),
                  numericInput("textWidthfour", "width:", 800),
                  actionButton("btnSubmit4", "Submit", class = "btn-primary"),
                  downloadButton("eleven-four", "instruction manual")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "Result", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot4", height = "100%", width = "100%"),
                  downloadButton("downloadResultfour", "Download Results")
                )
              ),
              fluidRow(
                box(title = "Data Table", 
                    status = "primary",
                    width = 12,
                    dataTableOutput("geneDataTable4"),
                    downloadButton("downloadTablefour", "Download Table")
                )
            )
      ),
      tabItem(tabName = "six",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  class = "custom-box",
                  p(class = "custom-text", "immune cell Gene Related"),
                  p(class = "text2", "Assess the correlation between immune cells 
                  and the input genes, with example data available for download.")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "Upload Data File GENE names", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectGene6", "GENENAME:", choices = NULL, multiple = TRUE),
                  numericInput("textHeightsix", "height:",500),
                  numericInput("textWidthsix", "width:",800),
                  actionButton("btnSubmit6", "Submit", class = "btn-primary"),
                  downloadButton("eleven-six", "instruction manual")
                )
              ),
              box(
                width = 12,
                title = "Result", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot6", height = "auto"),
                downloadButton("downloadResultsix", "Download Results")
              ),
              fluidRow(
                box(title = "Data Table", 
                    width = 12,
                    status = "primary", 
                    dataTableOutput("geneDataTable6"),
                    downloadButton("downloadTablesix", "Download Table")
                )
              )
      ),
      tabItem(tabName = "eight",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  class = "custom-box",
                  p(class = "custom-text", "proteomics"),
                  p(class = "text2", "To examine the expression of specific proteins 
                    in various organs pre- and post-exercise.")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "protein input", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectGeneeight", "proteinNAME:", choices = NULL, multiple = F),
                  numericInput("textHeighteight", "height:",500),
                  numericInput("textWidtheight", "width:",800),
                  actionButton("btnSubmit8", "Submit", class = "btn-primary"),
                  downloadButton("eleven-eight", "instruction manual")
                )
              ),
              fluidRow(
                box(
                  title = "Result", solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  plotOutput("plot8", height = "100%", width = "100%"),
                  downloadButton("downloadResulteight", "Download Results")
                )
              ),
              fluidRow(
                box(
                  title = "Table result",
                  width = 12,
                  dataTableOutput("proteinDataTableeight"),
                  downloadButton("downloadTableeight", "Download Table")
                )
              )
      ),
      tabItem(tabName = "nine",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  class = "custom-box",
                  p(class = "custom-text", "phosphoproteomic"),
                  p(class = "text2", "Analyze the changes in phosphorylation sites 
                    before and after exercise.")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "ProteinNAME input", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectGenenine", "ProteinNAME:", choices = NULL, multiple = F),
                  selectInput("selectDatasetnine", "Choose a dataset:", choices = c("cortex", "gastrocnemius", "heart", "kidney","liver","lung","white-adipose")),
                  numericInput("textHeightnine", "height:",500),
                  numericInput("textWidthnine", "width:",800),
                  actionButton("btnSubmit9", "Submit", class = "btn-primary"),
                  downloadButton("eleven-nine", "instruction manual")
                )
              ),
              # 新增的数据集选择框
              fluidRow(
                box(
                  title = "Result", solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  plotOutput("plot9", height = "100%", width = "100%"),
                  downloadButton("downloadResultnine", "Download Results")
                ),
              box(
                title = "Data Table", 
                status = "primary", 
                width = 12,
                dataTableOutput("geneDataTable9"),
                downloadButton("downloadTablenine", "Download Table" )
              )
            )
              
      ),
      tabItem(tabName = "c",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "Multi-omics"),
                    p(class = "text2", "This module visualizes multi-omics data from an eight-week exercise intervention in 
                      mice, comparing mRNA and protein expression changes across 12 organs with customizable experimental conditions.")
                )
              ),
              fluidRow(
                box(
                  title = "Difference analysis of datasets", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectDataset_c", "Choosing Tissues:",
                              choices = NULL,multiple = F),
                  actionButton("btnSubmit_c", "Step 1", class = "btn-primary"),
                  selectInput("selectDataset_c.1", "Control Group:",
                              choices = NULL,multiple = F),
                  selectInput("selectDataset_c.2", "Experimental Group:",
                              choices = NULL,multiple = F),
                  numericInput("Threshold_c", "Threshold logFC:", 0.5),
                  numericInput("textHeight_c", "height:", 500),
                  numericInput("textWidth_c", "width:", 500),
                  actionButton("btnSubmit_c.1", "Step 2", class = "btn-primary"),
                  downloadButton("eleven-c", "instruction manual")
                ),
                box(
                  title = "Result", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot_c", height = "auto"),
                  downloadButton("downloadResult_c", "Download Plot")
                ),
                box(
                  title = "Data Table",
                  status = "primary",
                  width = 12,
                  collapsible = TRUE,
                  dataTableOutput("Table_c"),
                  downloadButton("downloadTable_c", "Download Table")
                )
              )
      ),
      tabItem(tabName = "a",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "SportEnrich-GSEA"),
                    p(class = "text2", "This module performs GSEA using 120 exercise-related gene sets to 
                      reveal how different exercise types (e.g., aerobic or resistance training) influence molecular changes, automatically 
                      generating enrichment plots from uploaded sequencing data.")
                )
              ),
              fluidRow(
                box(title = "Upload the CSV file", status = "primary", solidHeader = TRUE,
                    width = 6,
                    fileInput("file_upload_a", "Select a file", accept = c(".csv")),
                    verbatimTextOutput("file_check_a"),
                    selectInput("selectDataset_a", "Species", choices = c("human", "mouse")),
                    selectInput("selectterm_a", "Enriched terms:", choices = NULL, multiple = F),
                    numericInput("textHeight_a", "Plot height (px):", 500),
                    numericInput("textWidth_a", "Plot width (px):", 500),
                    actionButton("btnSubmit_a", "Submit", class = "btn-primary"),
                    downloadButton("eleven-a", "instruction manual"),
                    downloadButton("eleven-a_1", "Example File")
                ),
                box(
                  title = "Result", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot_a", height = "auto"),
                  downloadButton("downloadResult_a", "Download Plot")
                ),
                box(
                  title = "Data Table",
                  status = "primary",
                  width = 12,
                  collapsible = TRUE,
                  dataTableOutput("Table_a"),
                  downloadButton("downloadTable_a", "Download Table")
                )
              )
      ),
      tabItem(tabName = "b",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    p(class = "custom-text", "Exercise-Drug Convergence"),
                    p(class = "text2", "This module integrates exercise and drug 
                      intervention multi-omics data to identify shared molecular 
                      pathways and potential exercise-drug interaction networks.")
                )
              ),
              fluidRow(
                box(
                  title = "Difference analysis of datasets", status = "warning",
                  solidHeader = TRUE, collapsible = TRUE,
                  selectInput("selectDataset_b", "Choosing Tissues:",
                              choices = NULL,multiple = F),
                  actionButton("btnSubmit_b", "Step 1", class = "btn-primary"),
                  selectInput("selectDataset_b.1", "Control Group:",
                              choices = NULL,multiple = F),
                  selectInput("selectDataset_b.2", "Experimental Group:",
                              choices = NULL,multiple = F),
                  fileInput("file_upload_b", "Select a file", accept = c(".csv")),
                  verbatimTextOutput("file_check"),
                  numericInput("textHeight_b", "height:", 250),
                  numericInput("textWidth_b", "width:", 250),
                  actionButton("btnSubmit_b.1", "Step 2", class = "btn-primary"),
                  downloadButton("eleven-b", "instruction manual"),
                  downloadButton("eleven-b_1", "Example File")
                ),
                box(
                  title = "Result", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot_b", height = "auto"),
                  downloadButton("downloadResult_b", "Download Plot")
                ),
                box(
                  title = "Data Table",
                  status = "primary",
                  width = 12,
                  collapsible = TRUE,
                  dataTableOutput("Table_b"),
                  downloadButton("downloadTable_b", "Download Table")
                )
              )
      )
      
    )
  )
)

