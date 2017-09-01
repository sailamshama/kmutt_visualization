curvy <- list(7, 74, 115, 228, 254)
normal <- list(358, 361, 367, 372, 392)
plist =list()

normal_names = list("Regulation of Acin Cytoskeleton", "G13 signaling athway","Glycolysis and Gluconeogenesis", "Complement and Coagulation Cascades KEGG", "Toll-like receptor signaling pathway (for GENMAPP)")
curvy_names = list("Galactose metabolism", "Pyruvate metabolism", "ABC transporters-general", "IL-2 signaling", "Translation Factors")
names <- read.csv("E:/Visualization_material/NewPathwayAPI_withName.txt", sep="\t")
#names <- read.csv("F:/Visualization_material/NewPathwayAPI_withName.txt", sep="\t")
names = names[, 1:2]


for (i in curvy){
  
  directory <- paste("E:/Visualization_material/lean_data/",i ,"/", sep="")
  #directory <- paste("F:/Visualization_material/lean_data/",i ,"/", sep="")
  filename <- paste(directory,i,"_sig-non.txt", sep="")
  myGraph <- read.table(filename, header = T)	
  #suppressWarnings(myGraph$Weight <- as.numeric(as.character(myGraph$Weight)))
  myGraph$Weight <- as.numeric(levels(myGraph$Weight))[myGraph$Weight]
  #geneset_name = as.character(names[which(names[,1] == i), 2])
  geneset_name <- curvy_names[which(curvy == i)][[1]]
  
  plot <- ggplot(myGraph, aes(x=No, y=Weight, color = Significant))+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_point()+
    # xlab("Gene-pair rank by weight")+
    #  ylab("Weight")+
    #        ggtitle(paste(geneset_name, ": Weights",sep=""))+
    #ggtitle(curvy_names[which(curvy == i)])+
    ggtitle(geneset_name)+
    theme(plot.title = element_text(face = "bsold", hjust = 0.5), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          legend.position="none", 
          #legend.direction = "horizontal",
          legend.title = element_blank()
    )
  plist[[i]] <- plot
  
}

#pdf(paste("E:/Visualization_material/","poster_plots_one_legend",i,".pdf", sep=""))
#grid.arrange(plot1, plot5)
#dev.off()
 
plist <- list(p1, p2, p3, p4, p5)

grid.arrange(grobs = plist, ncol = 1) ## display plot
ggsave(file = "E:/Visualization_material/poster_plots_loop_test.pdf", arrangeGrob(grobs = plist, ncol = 1))  ## save plot

