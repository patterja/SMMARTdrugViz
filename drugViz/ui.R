# Define the UI
ui <- fluidPage(
  titlePanel("Erlotinib response based on DUOX1 "),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("celllines", label = "Filter by cell line study", 
                         choices = list("CCLE" = "CCLE", "GDSC"="GDSC","GDSC1000" = "GDSC1000", "GCSI" = "GCSI", "FIMM"="FIMM", "CTRPV2"="CTRPV2"),
                         selected = c("CCLE","GDSC","GDSC1000","GCSI","FIMM","CTRPV2")),
      selectInput("SortByGeneExp", label = "Sort by increasing gene TPM or drug sensitivity",
                   choices = list("EGFR" = "egfr", "DUOX" = "duox", "median Sensitivity"="medianSensitivity"), 
                   selected = "duox"),
      hr(),
      downloadButton("report", "Generate report")),
    
      
    mainPanel(fluidPage(
      plotOutput("heatmap", height=310),
      plotOutput("scatter"),
      tableOutput("table")
    ))
    
  )
  
)