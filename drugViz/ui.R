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
                   selected = "duox")
      ),
      mainPanel(
        plotOutput("heatmap"),
        plotOutput("scatter")
        )
    ),
    sidebarLayout(
      sidebarPanel(  
      selectInput("ttest", label = "Separate into two groups for T-TEST",
                  choices = list("DUOX expression" = "duox", "EGFR expression" = "egfr", "variant"="variant"), 
                  selected = "duox"),
      conditionalPanel(
        condition="input.ttest == 'duox'",
        sliderInput("tpm", "TPM:",
                    min = min(na.omit(exp89$CCLE)), max = max(na.omit(exp89$CCLE)),
                    value = mean(na.omit(exp89$CCLE))))
      ),
    
      mainPanel(
        plotOutput("bxplt"),
        verbatimTextOutput("stattable"),
        radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                     inline = TRUE),
        downloadButton("report", "Download")
    ))
  
)