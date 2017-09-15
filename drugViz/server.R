server <- function(input, output) {

  
duoxexp_filt = exp89[exp89$cellline %in% auc_erlo$cell_line,]
duox_order = duoxexp_filt$cellline[c(order(duoxexp_filt$median_duox))]

#filter egfr and sensitivities inside duox list
egfrexp_filt = expEgfr[expEgfr$cellline %in% duox_order,]
egfr_order = egfrexp_filt$cellline[c(order(egfrexp_filt$median_egfr))]
aucerlo_study = auc_erlo[auc_erlo$cell_line %in% duox_order,]
med=apply(aucerlo_study[,c(2:ncol(aucerlo_study))], 1, medianWithoutNA)
sens_order= aucerlo_study$cell_line[c(order(med))]

# Trying to input reactivity basedon observation. Getting an error with this step `aucErlo_filt$cell_line[order(medianSens)])`. 
# Error is something like Error in order: unused argument (medianSens)
# Fixed above by just sorting by median of all sensitivities combined and not the study filtered median
# sens_order <- reactive({
#   aucErlo_filt = aucerlo_study[,append(c("cell_line"), input$celllines)]
#   medianSens = apply(aucErlo_filt[,c(2:ncol(aucErlo_filt))], 1, medianWithoutNA)
#   print(c(aucErlo_filt$cell_line[order(medianSens)]))
# })
# sens_order <- reactiveValues(data = NULL)
# 
# 
# observeEvent(input$celllines, {
#              aucErlo_filt = aucerlo_study[,append(c("cell_line"), input$celllines)]
#              medianSens = apply(aucErlo_filt[,c(2:ncol(aucErlo_filt))], 1, medianWithoutNA)
#              sens_order$data <-aucErlo_filt$cell_line[c(order(medianSens))]
# })
# sens_order <- eventReactive(input$celllines, {
#   aucErlo_filt = aucerlo_study[,append(c("cell_line"), input$celllines)]
#   medianSens= apply(aucErlo_filt[,c(2:ncol(aucErlo_filt))], 1, medianWithoutNA)
#   sord=aucErlo_filt$cell_line[c(order(medianSens))]
#   print(sord)
#   })

#Ordering
order <- reactive({
  switch(input$SortByGeneExp,
                 egfr = egfr_order,
                 duox = duox_order,
                 medianSensitivity= sens_order)
  })
  


output$heatmap <- renderPlot({
  duoxexp_filt$cellline = factor(duoxexp_filt$cellline, levels=order())
  egfrexp_filt$cellline = factor(egfrexp_filt$cellline, levels=order())
  mel=melt(duoxexp_filt[,c(1:4)], id.vars="cellline", measure=c("GRAY","UHN","CCLE"))
  mel$value[which(mel$value<1e-10)]=1e-10
  mee <- melt(egfrexp_filt[,c(1:4)], id.vars="cellline", measure=c("GRAY","UHN","CCLE"))
  mee$value[which(mee$value<1e-10)]=1e-10
  
  duox_erlo_heatmap <- ggplot(mel,  aes(cellline, variable)) +
    geom_tile(aes(fill = value)) +
    labs(y = "DUOX1", x=element_blank()) +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color="black"),
          plot.margin=unit(c(0,0,-20,1), "mm"), 
          legend.position = "top") +
    scale_fill_gradientn(name=expression(paste("log"[2],"TPM")), 
                         colours=c("navyblue","darkmagenta","darkorange1"), 
                         na.value="transparent", 
                         trans = 'log2',
                         breaks=c(9.313226e-10,9.536743e-07,0.0009765625,1),
                         labels=c("-30","-20","-10","1")) +
    coord_fixed(ratio = 1)
  
  egfr_erlo_heatmap <- ggplot(mee,  aes(cellline, variable)) +
    geom_tile(aes(fill = value)) +
    labs(y = "EGFR", x="cell line") +
    theme(panel.background = element_blank(),
          plot.margin=unit(c(-20,0,0,1), "mm"), 
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
  print(grid.arrange(duox_erlo_heatmap, egfr_erlo_heatmap, ncol=1, nrow=2))
  })

#-------------------------------------------------
#Filter for duox order


output$scatter <- renderPlot({
  aucErlo_filt = aucerlo_study[,append(c("cell_line"), input$celllines)]
  aucErlo_filt$cell_line = factor(aucErlo_filt$cell_line, levels=order())
  
  
  mec_auc=melt(aucErlo_filt, id.vars="cell_line", measure=colnames(aucErlo_filt[,-1]))
  erlo_scatter <- ggplot(mec_auc, aes(x=cell_line, y=value)) +
    geom_point(aes(colour=factor(variable))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.background = element_rect(fill="white", colour = "black"),
          axis.text=element_text(colour="black"),
          legend.title = element_blank(),
          panel.grid.major.x = element_line( size=.1, color="gray" ), 
          legend.position="bottom") +
    scale_colour_brewer(palette = "Dark2") +
    labs(y="AAC", x="cell lines\nordered by increasing DUOX1 expression") +
    ylim(0, 0.3)
  print(erlo_scatter)
})

}