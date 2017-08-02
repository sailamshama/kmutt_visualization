library(ggplot2)

names <- read.csv("/Users/saila/Desktop/Visualization_material/NewPathwayAPI_withName.txt", sep="\t")
names = names[, 1:2]
  
for (i in 3:10){
  directory <- paste("/Users/saila/Desktop/Visualization_material/lean_data/",i , sep="")
  filename <- paste("/Users/saila/Desktop/Visualization_material/lean_data/",i,"/",i,"_sig-non.txt", sep="")
  if (file.exists(filename)){
    geneset_name = as.character(names[which(names[,1] == i), 2])
    myGraph <- read.table(filename, header = T)	##Set path
    suppressWarnings(myGraph$Weight <- as.numeric(as.character(myGraph$Weight)))
    
    png(paste(directory, "/profile",i,".png", sep=""))
    ggplot(myGraph, aes(x=No, y=Weight))+
      #scale_x_continuous(limits=c(0,600))
      geom_point(aes(colour = Significant), show.legend=TRUE)+
      xlab("gene-pair rank by weight")+
      ylab("weight")+
      ggtitle(paste(geneset_name, ": Weights",sep=""))+
      theme(plot.title = element_text(hjust = 0.5))
    dev.off()
                
    png(paste(directory, "/sig_pairs",i,".png", sep="")) 
    ggplot(myGraph, aes(x=Significant, y=Weight))+
      geom_point(aes(colour=Significant), show.legend = TRUE)+
      xlab("")+
      ylab("weight")+
      ggtitle(paste(geneset_name, ": Dissimilarity between gene-pair",sep=""))+
      theme(plot.title = element_text(hjust = 0.5))
    dev.off()
  }
}



#========================================================================
#Noon's graphing


##To count the number of each group (frequency table)
#table(myGraph$Significant)

# weights <- myGraph$Weight
# pair_rank <- myGraph$No
# sig_status <- myGraph$Significant


# p <- qplot(Significant, Weight, colour = Significant, data = head(myGraph), 
#            xlab="Gene-pair", ylab="Weight")  
# p+theme(axis.title=element_text(face="bold", size="10", color="black"), legend.position="top")
# 
# 
# qplot(No, Weight, colour = Significant, data = myGraph[1:500,], 
#       main="Gene-pair profile", xlab="Gene-pair", ylab="Weight") 	