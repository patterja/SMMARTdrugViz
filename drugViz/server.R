server <- function(input, output) {

  
duoxexp_filt = exp89[exp89$cellline %in% auc_erlo$cell_line,]
duox_order = duoxexp_filt$cellline[c(order(duoxexp_filt$median_duox))]

#filter egfr and sensitivities inside duox list
egfrexp_filt = expEgfr[expEgfr$cellline %in% duox_order,]
egfr_order = egfrexp_filt$cellline[c(order(egfrexp_filt$median_egfr))]
aucerlo_study = auc_erlo[auc_erlo$cell_line %in% duox_order,]
med=apply(aucerlo_study[,c(2:ncol(aucerlo_study))], 1, medianWithoutNA)
sens_order= aucerlo_study$cell_line[c(order(med))]



#Ordering
order <- reactive({
  switch(input$SortByGeneExp,
                 egfr = egfr_order,
                 duox = duox_order,
                 medianSensitivity= sens_order)
  })
#----------------------------------------

mel <- reactive({
  duoxexp_filt$cellline = factor(duoxexp_filt$cellline, levels=order())
  mel=melt(duoxexp_filt[,c(1:4)], id.vars="cellline", measure=c("GRAY","UHN","CCLE"))
  mel$value[which(mel$value<1e-10)]=1e-10
  mel
})

mee <- reactive({
  egfrexp_filt$cellline = factor(egfrexp_filt$cellline, levels=order())
  mee = melt(egfrexp_filt[,c(1:4)], id.vars="cellline", measure=c("GRAY","UHN","CCLE"))
  mee$value[which(mee$value<1e-10)]=1e-10
  mee
})

duox_erlo_heatmap <- reactive({
  ggplot(mel(),  aes(cellline, variable)) +
  geom_tile(aes(fill = value)) +
  labs(y = "DUOX1", x=element_blank()) +
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        plot.margin=unit(c(0,0,-1,1), "mm"), 
        legend.position = "top") +
  scale_fill_gradientn(name=expression(paste("log"[2],"TPM")), 
                       colours=c("navyblue","darkmagenta","darkorange1"), 
                       na.value="transparent", 
                       trans = 'log2',
                       breaks=c(9.313226e-10,9.536743e-07,0.0009765625,1),
                       labels=c("-30","-20","-10","1")) +
  coord_fixed(ratio = 1)
})

egfr_erlo_heatmap <- reactive({
  ggplot(mee(),  aes(cellline, variable)) +
  geom_tile(aes(fill = value)) +
  labs(y = "EGFR", x="cell line") +
  theme(panel.background = element_blank(),
        plot.margin=unit(c(-1,0,-1,1), "mm"), 
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  scale_fill_gradientn(name=expression(paste("log"[2],"TPM")), 
                       colours=c("navyblue","darkmagenta","darkorange1"), 
                       na.value="transparent", 
                       trans = 'log2',
                       breaks=c(9.313226e-10,9.536743e-07,0.0009765625,1),
                       labels=c("-30","-20","-10","1")) +
  coord_fixed(ratio = 1)
})


#-------------------------------------------------
#SCATTER
mec_auc = reactive({
  aucErlo_filt=aucerlo_study[,append(c("cell_line"), input$celllines)]
  aucErlo_filt$cell_line = factor(aucErlo_filt$cell_line, levels=order())
  mec_auc=melt(aucErlo_filt, id.vars="cell_line", measure=colnames(aucErlo_filt[,-1]))
  mec_auc
  })

erlo_scatter <- reactive({
  erlo_scatter <- ggplot(mec_auc(), aes(x=cell_line, y=value)) +
  geom_point(aes(colour=factor(variable))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill="white", colour = "black"),
        axis.text=element_text(colour="black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_line( size=.1, color="gray" ), 
        legend.position="bottom",
        plot.margin=unit(c(-1,0,0,1), "mm")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(y="AAC", x=paste0("cell lines\nordered by increasing ", toupper(input$SortByGeneExp))) +
  ylim(0, 0.3)
})
#-------------------------------------------------
#BOXPLOT



sepExpCCLE <- reactive({
  switch(input$ttest,
         duox = mel()[mel()$variable=="CCLE",],
         egfr = mee()[mee()$variable=="CCLE",],
         variant = c())
})


less <- reactive({
  l = sepExpCCLE()$cellline[sepExpCCLE()$value < input$tpm]
  na.omit(mec_auc()$value[mec_auc()$cell_line %in% l])
  
})

more <- reactive({
  m = sepExpCCLE()$cellline[sepExpCCLE()$value >= input$tpm]
  na.omit(mec_auc()$value[mec_auc()$cell_line %in% m])
  
})

#-------------------------------------------------
#ALL OUTPUT PLOTS
output$heatmap <- renderPlot({
  layout=matrix(c(1:3,3,3),5,1)
  pl=c(duox_erlo_heatmap(), egfr_erlo_heatmap(), erlo_scatter())
  grid.arrange(duox_erlo_heatmap(), egfr_erlo_heatmap(), erlo_scatter(), nrow=3, heights = c(0.4, 0.4,1))
}, height = 600)


output$bxplt <- renderPlot({
  boxplot(less(),more())
})

output$stattable <-renderPrint({
  t.test(less(), more())
})

#-------------------------------------------------
#DOESN'T WORK!!!!
output$report <- downloadHandler(
  filename = function() {
    paste('my-report', sep = '.', switch(
      input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
    ))
  },
  
  content = function(file) {
    src <- normalizePath('report.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    file.copy(src, 'report.Rmd', overwrite = TRUE)
    
    library(rmarkdown)
    out <- render('report.Rmd', switch(
      input$format,
      PDF = pdf_document(), HTML = html_document(), Word = word_document()
    ))
    file.rename(out, file)
  }
)


  }
