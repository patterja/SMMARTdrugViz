# Define the UI
ui <- fluidPage(
  titlePanel("Erlotinib response based on DUOX1 "),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("celllines", label = "Filter by cell line study", 
                         choices = list("CCLE" = "CCLE", "GDSC"="GDSC","GDSC1000" = "GDSC1000", "GCSI" = "GCSI", "FIMM"="FIMM", "CTRPV2"="CTRPV2"),
                         selected = c("CCLE","GDSC","GDSC1000","GCSI","FIMM","CTRPV2")),
      radioButtons("SortByGeneExp", label = "Sort by increasing TPM or Sensitivity",
                   choices = list("EGFR" = "egfr", "DUOX" = "duox", "median Sensitivity"="medianSensitivity"), 
                   selected = "duox"),
      hr()),
      
    mainPanel(fluidPage(
      plotOutput("heatmap", height=310),
      plotOutput("scatter"),
      tableOutput("table")
    ))
    
  )
  
)