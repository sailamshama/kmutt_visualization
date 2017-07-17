## To construct a gene-set-based gene co-expression network
## Each gene-set will be used to construct a fully-connected network
## "NewPathwayAPI_forGCN.txt" is gene-set data
## "GSE10072_prepro_nonredun_var_noProbe_z.txt" is preprocessed expression data

## Gene Co-expression Network for GSNFS method
## Necessary package for GCN construction
install.packages("WGCNA")
install.packages("knitr", dependencies = TRUE)
library(WGCNA)
install.packages("igraph")
library(igraph)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()




# Read each gene-set from PathwayAPI file
# Read in the data
inputList <- scan("NewPathwayAPI_forGCN.txt", what="", sep="\n")
# Separate elements by one or more whitepace
myList <- strsplit(inputList, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(myList) <- sapply(myList, `[[`, 1)
#names(myList) <- sapply(myList, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
myList <- lapply(myList, `[`, -1)
#myList <- lapply(myList, function(x) x[-1]) # same as above
## use lapply() to return list object, put ; symbol to gene name to let the search term specific 100% not just subset
geneSet <- lapply(myList,function(x) paste(";",x, ";", sep = ""))            


# From .txt file store information of which ID connect to which ID gene
expFile  <- as.matrix(read.table("GSE10072_prepro_nonredun_var_noProbe_z.txt", header =T, sep= "\t"))     
expFile1 <- as.matrix(paste(";", expFile[,1],";", sep=""))
expFiles <- as.matrix(cbind(expFile1,expFile[,2:ncol(expFile)]))


SeparateExp <-function(geneSet, expFiles){
	for (i in 1:length(geneSet)){	
		lgn <- geneSet[i]
		nameGS<-names(lgn)
		lg <- geneSet[[i]]
		sigMatch <- NULL
		sigMatchNet  <- NULL		
		
		for (j in 1:length(lg)){	
			sigMatch <- expFile[grep(lg[j], expFiles[,1]),]
			sigMatchNet <- rbind(sigMatchNet, sigMatch)
			}		
		#write.table(sigMatchNet, file = paste("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\GeneCoNet\\GeneSetExp\\" ,i, "_GSallExp_10072.txt", sep=""), row.names =F, quote=F, sep = "\t")
		write.table(sigMatchNet, file = "GSallExp_10072.txt", row.names =F, quote=F, sep = "\t")
		
		### Construct GCN of each gene-set
		inputFile = as.data.frame(read.table("GSallExp_10072.txt", header = T, sep = "\t"));  
		datExpr = as.matrix(t(inputFile[,2:ncol(inputFile)]));

		names(datExpr) = inputFile$GeneName;  # To set column name as GeneName
		nGenes = ncol(datExpr)
		nSamples = nrow(datExpr)
		powers = c(c(1:10), seq(from = 12, to=30, by=2))
		sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

		
		# To draw and save graph 
		
		sizeGrWindow(9, 5)
		par(mfrow = c(1,2));
		cex1 = 0.9;
		
		png(file = paste("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\GeneCoNet\\GeneCoNet_eachGeneSet\\withWeight\\graph\\",nameGS, "_R^2.png", sep=""))		
		plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
		     main = paste("Scale independence"));
		text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		     labels=powers,cex=cex1,col="red");
		abline(h=0.90,col="red")
		dev.off()
		
		png(file = paste("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\GeneCoNet\\GeneCoNet_eachGeneSet\\withWeight\\graph\\",nameGS, "_MeanCon.png", sep=""))
		plot(sft$fitIndices[,1], sft$fitIndices[,5],
		     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
		     main = paste("Mean connectivity"))
		text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
		dev.off()
		
		
		# Store output from pickSoftThreshold into outsft
		outsft<- as.matrix(sft$fitIndices)
		# Find Power (for softPower) given the maximum value of SFT.R.sq from outsft
		sp <- outsft[match(max(outsft[,2]),outsft[,2]),1]
		
		print ("sp = ")
		print (sp)

		####Change number of softpower
		softPower = sp;   # change the number according to the appropriate power
		adjacency = adjacency(datExpr, power = softPower);
		TOM = TOMsimilarity(adjacency);
		dissTOM = 1-TOM
		colnames(dissTOM) = inputFile$GeneName
		rownames(dissTOM) = inputFile$GeneName

		#Construct network
		geneCoNet = graph.adjacency(dissTOM, weighted=TRUE, mode="upper") 
		
		#write graph
		# write file that store pair of interaction (node id start from zero not 1)
		write.graph(geneCoNet, file="gCoNet_10072.txt", format="edgeList")   
		geneCoNetWeight <- E(geneCoNet)$weight    
		geneCoNetNodeName <- V(geneCoNet)$name
		# weight of the pair of interaction in network file (geneCoNet) start from the first pair
		write.csv(geneCoNetWeight, file = "gCoNetWeight_10072.txt")        
		# list of node (according to node ID)
		write.csv(geneCoNetNodeName, file = "gCoNetNodeName_10072.csv")	


		# To construct pairwise gene-gene interaction from gcNet.txt file and gcNodeName.csv 
		# .csv file contains ID and gene name
		nodeNameIn  <- as.matrix(read.table("gCoNetNodeName_10072.csv", header = T, sep= ","))       	
		nodeNameIn  <- gsub(" ", "", nodeNameIn)
		## put ; symbol to gene name to let the search term specific 100% not just subset
		nodeNameAdd <- as.matrix(paste(";", nodeNameIn[,1],";", sep=""))     
		nodeName    <- as.matrix(cbind(nodeNameAdd,nodeNameIn[,2]))
		# From .txt file -- store information of which ID connect to which ID gene
		gcNet  <- as.matrix(read.table("gCoNet_10072.txt", header =F, sep= " "))      
		gcNet <- apply(gcNet, 1:2, function(x) x+1)
		gcNet1 <- as.matrix(paste(";", gcNet[,1],";", sep=""))
		gcNet2 <- as.matrix(paste(";", gcNet[,2],";", sep=""))
		gcNets <- as.matrix(cbind(gcNet1,gcNet2))

		FindGene  <- function(nodeName, gcNets,...){
		nodeMatch <- NULL
		nodeMatchNet  <- NULL	
			for (i in 1:nrow(gcNets)){	
			nodeMatch1 <- nodeName[grep(gcNets[i,1], nodeName[,1]),2]
			nodeMatch2 <- nodeName[grep(gcNets[i,2], nodeName[,1]),2]
			nodeMatch  <- cbind(nodeMatch1,nodeMatch2)
			nodeMatchNet <- rbind(nodeMatchNet, nodeMatch)
			}
			return(nodeMatchNet)
		}
		findGC <- FindGene(nodeName,gcNets)
		# Add weight column to the gcNet
		gcWeight <- as.matrix(read.table("gCoNetWeight_10072.txt", header = T, sep= ",")) 
		findG <- cbind(findGC, gcWeight[,2]) 
		
		write.table(findG, file = paste("D:\\GNFS_modi\\CoExpress\\WGCNA_file\\VIS\\NewGeneSet\\" ,nameGS, ".txt",, sep=""), 
			row.names =F, col.names = F, quote=F, sep = "\t")
				
			
	## End of GCN construction	
		print ("---------------------------------------------")
		print (nameGS)
		}
		
}
fgs<-SeparateExp(geneSet, expFiles)
	


	