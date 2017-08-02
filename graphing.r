library(ggplot2)

i  = 0
directory <- paste("/Users/saila/Desktop/Visualization_material/lean_data/",i , sep="")
setwd(directory)

filename <- paste("/Users/saila/Desktop/Visualization_material/lean_data/",i,"/",i,"_sig-non.txt", sep="")
myGraph <- read.table(filename, header = T)	##Set path
myGraph$Weight <- as.numeric(as.character(myGraph$Weight))

png(paste(directory, "/profile",i,".png", sep=""))
ggplot(myGraph, aes(x=No, y=Weight))+
  geom_point(aes(colour = Significant), show.legend=TRUE)
dev.off()
            
png(paste(directory, "/sig_pairs",i,".png", sep="")) 
ggplot(myGraph, aes(x=Significant, y=Weight))+
  trio_graph+geom_point(aes(colour=Significant), show.legend = TRUE)
dev.off()

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