## -------- START --------
#######################################################
#inputFile = GCN network file having one and two column is gene pair and third column is weight, without header
#specify GCN tsv file
#create loop here for all tsv files for (i in 0:397) {}
i = 0
networkInput <- as.matrix(read.csv(paste("/Users/saila/Desktop/Visualization_material/GCN_353/", i,".tsv", sep=""), header= F, sep = "\t"))     ##Set path
#REMOVE HEADER
netInCol1  <- as.matrix(paste(":", networkInput[,1],":", sep=""))
netInCol2  <- as.matrix(paste(":", networkInput[,2],":", sep=""))
netWeight <- as.data.frame(cbind(netInCol1, netInCol2, networkInput[,3]))

#######################################################
#Significant genes of each gene-set record of each dataset = all genes member of each gene-set record are significant gene
#no need to update path
inputSigList <- scan("/Users/saila/Desktop/Visualization_material/All_sig_gene_10072.txt", what="", sep="\n")	##Set path
#Separate elements by one or more whitepace
myList <- strsplit(inputSigList, "[[:space:]]+")
#Extract the first vector element and set it as the list element name
names(myList) <- sapply(myList, `[[`, 1)
#Remove the first vector element from each list element
myList <- lapply(myList, `[`, -1)


###TO SPECIFY which gene-set will be used to draw a correlation graph

#######################################################
#must be same as file chosen in first path
drawList <- i
myL <- myList[]
nameGS<-names(myL)
matName <- charmatch(drawList, nameGS)
if(!is.na(matName)){
	ml<- myList[[matName]]
	mlm<-as.matrix(ml)
}          
sigG    <- as.matrix(paste(":", mlm[,1],":", sep=""))
#######################################################

##To specify the type of relationship of each gene-pair
FindExp <- function(sigG, netWeight,...){
geneMatch <- NULL
geneMatchSub <- NULL	
	for (j in 1:nrow(netWeight)){	
    match1 <- charmatch(netWeight[j,1], sigG[,1])
		if(!is.na(match1)){				
			match2 <- charmatch(netWeight[j,2], sigG[,1])	
			if(!is.na(match2)){
				a<-cbind(netWeight[j,],significant ="Sig-sig")
				geneMatchSub <- rbind(geneMatchSub, a)				
			}else{
				a<-cbind(netWeight[j,],significant ="Sig-non")
				geneMatchSub <- rbind(geneMatchSub, a)	
			}			
		}else{
			match2 <- charmatch(netWeight[j,2], sigG[,1])
			if(!is.na(match2)){
				a<-cbind(netWeight[j,],significant ="Sig-non")  #Should label "Non-sig" but it is categorized into same group as "Sig-non" 
				geneMatchSub <- rbind(geneMatchSub, a)				
			}else{
				a<-cbind(netWeight[j,],significant ="Non-non")
				geneMatchSub <- rbind(geneMatchSub, a)	
			}
		}				
	}
	return(geneMatchSub)
}
fg<-FindExp(sigG, netWeight)
fg <- as.matrix(fg)
fgc <-gsub(":", "", fg)
#create new file containing lean_data
## maybe create new directory
#mainDir <- "/Users/saila/Desktop/Visualization_material/lean_data/"
#subDir <- "64"
#dir.create(file.path(mainDir, subDir))
# https://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
mainDir <- "/Users/saila/Desktop/Visualization_material/lean_data/"
subDir <- drawList
dir.create(file.path(mainDir, subDir))
write.table(fgc, file = paste("/Users/saila/Desktop/Visualization_material/lean_data/", drawList, "/", drawList, "_sig-non.txt", sep =""),
	row.names = F, col.names = c("Gene1", "Gene2", "Weight", "Significant"), sep ="\t")

#######################################################

##### MUST DO WORK HERE #####
filename <- paste("/Users/saila/Desktop/Visualization_material/lean_data/",drawList,"/",drawList,"_sig-non.txt", sep="")
gene_pair_file <- as.matrix(read.csv(filename, sep="\t"))

## SORT WEIGHT COLUMN ##
weight_column <- gene_pair_file[, "Weight"]
sorted_weight_column <- sort(weight_column)
gene_pair_file[, "Weight"] <- sorted_weight_column

## ADD "RANK" COLUMN ##
new_column <- seq(1, nrow(gene_pair_file))

gene_pair_file_chunk1 <- gene_pair_file[, 1:2]
gene_pair_file_chunk2 <- gene_pair_file[, 3:4]

new_gene_pair_file <- cbind(gene_pair_file_chunk1, new_column, gene_pair_file_chunk2)
##name rank column as "No"
colnames(new_gene_pair_file) <- c("Gene1",	"Gene2","No",	"Weight",	"Significant")

write.table(new_gene_pair_file, filename, sep="\t")

#######################################################

###Manually add "No" column to specify the order of gene-pair based on weight that ascending sort


## -------- END preparation file for graph plotting --------



##Draw graph using "ggplot2" package in R
library(ggplot2)

##Read file to draw network
myGraph <- read.table(filename, header = T)	##Set path
##To count the number of each group (frequency table)
table(myGraph$Significant)


############ GET gene set name corresponding to gene set number ###############
##Graph with title which name depend on the input
#you may add pathname from NewPathwayAPI_withName.txt
#######################################################
#=========================================================================================

qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of MAPK signaling pathway", xlab="Gene-pair", ylab="Weight")
	
qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of ErbBsignaling pathway", xlab="Gene-pair", ylab="Weight")  
	
	
qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of MAPKsignaling pathway", xlab="Gene-pair", ylab="Weight")  
	
	
qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Non-small cell lungcancer pathway", xlab="Gene-pair", ylab="Weight")  
	

qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Endometrial cancer pathway", xlab="Gene-pair", ylab="Weight")  

	
qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Prostate cancer pathway", xlab="Gene-pair", ylab="Weight")  
	

qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Regulation of Actin Cytoskeleton pathway", xlab="Gene-pair", ylab="Weight")  

## Graph with title =========================================================================================	




## Plot the graph based on category according to Sig-sig, Sig-non, non-non group(128_sig-non.txt should be ordered by Significant column)
#Pattern I
myGraph <- read.table("/Users/saila/Desktop/Visualization_material/lean_data/192/192_sig-non.txt", header = T)
qplot(Significant, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Prostate cancer pathway", xlab="Gene-pair", ylab="Weight")  

#Pattern II
myGraph <- read.table("/Users/saila/Desktop/Visualization_material/lean_data/192/192_sig-non.txt", header = T)
p <- qplot(Significant, Weight, colour = Significant, data = myGraph, 
	 xlab="Gene-pair", ylab="Weight")  
p+theme(axis.title=element_text(face="bold", size="10", color="black"), legend.position="top")



### formatting graph ###

## Plot the graph based on order of the gene-pair by weights (128_sig-non.txt should be ordered by weight (No column))
#Pattern I
myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\128\\128_sig-non.txt", header = T)
qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Prostate cancer pathway", xlab="Gene-pair", ylab="Weight")  

#Pattern II
myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\128\\128_sig-non.txt", header = T)
p <- qplot(No, Weight, colour = Significant, data = myGraph, 
	 xlab="Gene-pair", ylab="Weight")  
p+theme(axis.title=element_text(face="bold", size="10", color="black"), legend.position="top", 
	legend.text=element_text(size=12))



# With graph title	
myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\128\\128_sig-non.txt", header = T)
p <- qplot(No, Weight, colour = Significant, data = myGraph, 
	main="MAPKsignaling pathway", xlab="Gene-pair", ylab="Weight")  
p+theme(axis.title=element_text(face="bold", size="10", color="black"), legend.position="bottom")
	

################################################################
	
myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\210\\210_sig-non.txt", header = T)
p <- qplot(Significant, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Apoptosis signaling", xlab="Gene-pair", ylab="Weight") 
p+theme(legend.text=element_text(size=12))

################################################################
	
	
	

## Read file to draw network 
myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\192\\192_sig-non.txt", header = T)
table(myGraph$Significant)


myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\210\\210_sig-non.txt", header = T)
qplot(No, Weight, colour = Significant, data = myGraph,  
	main=" Dissimilarity between gene pair of GCN of Apoptosis signaling", xlab="Gene-pair", ylab="Weight")  
	

	
myGraph <- read.table("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\74\\74_sig-non.txt", header = T)
qplot(Significant, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Pyruvate metabolism pathway", xlab="Gene-pair", ylab="Weight")  

qplot(No, Weight, colour = Significant, data = myGraph, 
	main=" Dissimilarity between gene pair of GCN of Pyruvate metabolism pathway", xlab="Gene-pair", ylab="Weight") 	

	

#When specify which gene-pair presented as the biomarker from GSNFS method	
myGraph[,'Biomarker'] <- as.factor(myGraph[,'Biomarker'])	
qplot(No, Weight, colour = Significant, data = myGraph, shape = Biomarker, 
	main=" Weight of GCN of Pyruvate metabolism", xlab="Gene-pair", ylab="Weight")  
	
	


# (NOt use) For adjust graph 
p + theme_bw()
p + theme(axis.title=element_text(face="bold.italic", 
   size="12", color="brown"), legend.position="top")
