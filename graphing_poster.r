library(ggplot2)
library(gridExtra)
library(extrafont)


names <- read.csv("/media/philipqlu/NO NAME/Visualization_material/NewPathwayAPI_withName.txt", sep="\t")
#names <- read.csv("F:/Visualization_material/NewPathwayAPI_withName.txt", sep="\t")
names = names[, 1:2]
  
for (i in list(64, 174)){
  if(i==192){}
  else{
    directory <- paste("/media/philipqlu/NO NAME/Visualization_material/lean_data/",i ,"/", sep="")
    #directory <- paste("F:/Visualization_material/lean_data/",i ,"/", sep="")
    filename <- paste(directory,i,"_sig-non.txt", sep="")
    if (file.exists(filename)){
      #geneset_name = as.character(names[which(names[,1] == i), 2])
      if (i == 64){geneset_name <- "Glycerophospholipid Metabolism"}
      else{geneset_name <- "Type II Diabetes Mellitus"}
      myGraph <- read.table(filename, header = T)	
      #suppressWarnings(myGraph$Weight <- as.numeric(as.character(myGraph$Weight)))
      myGraph$Weight <- as.numeric(levels(myGraph$Weight))[myGraph$Weight]
      
      ############## GENE PROFILE - to plot all points ###############################
#      png(paste(directory, "profile",i,".png", sep=""))
      plot1 <- ggplot(myGraph, aes(x=No, y=Weight, color = Significant))+
        scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        geom_point()+
        xlab("Gene-pair rank by weight")+
        ylab("Weight")+
#        ggtitle(paste(geneset_name, ": Weights",sep=""))+
        ggtitle(geneset_name)+
        theme(plot.title = element_text(face = "bold", hjust = 0.5), 
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"), 
              legend.position="none", 
#              legend.direction = "horizontal",
              legend.title = element_blank()
              )
      
      # tiff(paste("/Users/saila/Desktop/", geneset_name, "plot1.tiff", sep=""), width = 600)
      # print(plot1)
      # dev.off()
        
#      print(plot1) #important for saving while in loop
#      dev.off()
      
      ############ GENE PROFILE -  to plot gene-pairs ranked below 600####################
      # partGraph <- myGraph[1:600,]
      # partGraph$Weight <- as.numeric(levels(partGraph$Weight))[partGraph$Weight]
      # 
      # png(paste(directory, "profile",i,"_part.png", sep=""))
      # plot2 <- ggplot(partGraph, aes(x=No, y=Weight))+
      #   #scale_x_continuous(limits=c(0,600))
      #   geom_point(aes(colour = Significant), show.legend=TRUE)+
      #   xlab("gene-pair rank by weight")+
      #   ylab("weight")+
      #   ggtitle(paste(geneset_name, ": Weights",sep=""))+
      #   theme(plot.title = element_text(hjust = 0.5))
      # print(plot2) #important for saving while in loop
      # dev.off()
      
      ############ Box plot for each categories w/ MEDIAN  ####################
 #     myGraph$No <- as.factor(myGraph$No)
      medians <- aggregate(Weight ~ Significant, myGraph, median)
      medians$Weight <- round(medians$Weight, 2)
      
#       #plot median with boxplot
#       plot3 <- ggplot(myGraph, aes(x=Significant, y=Weight))+
#         geom_jitter(aes(colour=Significant), show.legend = FALSE, position=position_jitter(width=.2, height=0))+
#         geom_boxplot(outlier.shape = NA, width=0.2, show.legend = FALSE)+
# #        stat_summary(fun.y=median, geom="point", shape=8, size=4)+
#         geom_text(data = medians, aes(label = Weight, y = Weight - 0.02), show.legend = FALSE)+
#         xlab("")+
#         ylab("weight")+
#         ggtitle(paste(geneset_name, ": Dissimilarity in gene-pair",sep=""))+
#         theme(plot.title = element_text(hjust = 0.5))
#       
    
      ############ Strip plot for three categories w/ MEAN and STD. DEV.  ####################  
      #plot mean with std. dev. 
      means <- aggregate(Weight ~ Significant, myGraph, mean)
      means$Weight <- round(means$Weight, 2)
      
      stddevs <- aggregate(Weight ~ Significant, myGraph, sd)
      stddevs$Weight <- round(stddevs$Weight,2)
      
      # plot4 <- ggplot(myGraph, aes(x=Significant, y=Weight))+
      #   geom_jitter(aes(colour=Significant), show.legend = TRUE, position=position_jitter(width=.2, height=0))+
      #   stat_summary(fun.data="mean_sdl", mult=1, geom="pointrange", color="black")+
      #   geom_text(data = means, aes(label = Weight, y = Weight + 0.03))+
      #   xlab("")+
      #   ylab("weight")+
      #   ggtitle(paste(geneset_name, ": Dissimilarity in gene-pair",sep=""))+
      #   theme(plot.title = element_text(hjust = 0.5))
      # 
      
      #################### BASIC JITTER STRIP PLOT #################################
 #     png(paste(directory, "/sig_pairs",i,"_stat_summary.png", sep=""))
      plot5 <- ggplot(myGraph, aes(x=Significant, y=Weight, color = Significant))+
        #geom_boxplot(outlier.shape = NA, position="dodge")+
        geom_jitter(position=position_jitter(width=.2, height=0))+
        scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        stat_summary(fun.data="mean_sdl", mult=1, 
                     geom="pointrange", color="red")+
        geom_text(data = means, aes(label = Weight, y = Weight, fontface="bold"), nudge_x = 0.3, show.legend = FALSE, color="red")+
      
        stat_summary(fun.y=median, geom="point", shape=18,
                     size=3, color="black")+
        geom_text(data = medians, aes(label = Weight, y = Weight, fontface="bold"), nudge_x = 0.3, show.legend = FALSE, color="black")+
      
        geom_text(data = stddevs, aes(label = Weight, y = means$Weight - Weight), nudge_x = 0.3, show.legend = FALSE, color="red")+
        xlab("Gene pair type by significance")+
        ylab("Weight")+
#        ggtitle(paste(geneset_name, ": Dissimilarity in gene-pair",sep=""))+
        theme(#plot.title = element_text(face = "bold", hjust = 0.5), 
              axis.title.x = element_text(),
              axis.title.y = element_text(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"), 
              legend.position = "bottom",
              #legend.position=c(0.8, 0.1), 
              #legend.direction = "horizontal",
              legend.title = element_blank()
        )
      
      # tiff(paste("/Users/saila/Desktop/", geneset_name, "plot2.tiff", sep=""), width = 600)
      # print(plot5)
      # dev.off()
      
      # png("/Users/saila/Desktop/plot_draft.png", res =300)
      # grid.arrange(plot5)
      # dev.off()
      
  #    print(plot2)
  #    dev.off()
      
      pdf(paste("/media/philipqlu/NO NAME/Visualization_material/","poster_plots_one_legend",i,".pdf", sep=""))
      grid.arrange(plot1, plot5)
      dev.off()
    }
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