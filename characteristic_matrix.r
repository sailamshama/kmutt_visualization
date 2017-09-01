device_directory <- "F:/Visualization_material/"


names <- read.csv(paste(device_directory,"NewPathwayAPI_withName.txt", sep=""), sep="\t")
names = names[, 1:2]


characteristic_file <- read.csv(paste(device_directory, "clusters/curvy/complete_plots/curvy_genesets_characteristics.csv", sep=""))
curvy_geneset_nos <- characteristic_file[,"geneSetNo"]
  
for (i in curvy_geneset_nos){
  
  gene_set_row_num <- which(characteristic_file[,"geneSetNo"] == i)
  
  #read data for gene-set
  gene_set_data <- read.csv(paste(device_directory, "lean_data/", i, "/", i,"_sig-non.txt", sep=""), sep="\t")
  gene_set_data$Weight <- as.numeric(levels(gene_set_data$Weight))[gene_set_data$Weight]
  
  #get characteristics
  gene_pair_count <- nrow(gene_set_data)
  
  counts <- table(gene_set_data[,"Significant"])
  
  sigsig_count <- counts["Sig-sig"][1]
  signon_count <- counts["Sig-non"][1]
  nonnon_count <- counts["Non-non"][1]
  
  medians <- aggregate(Weight ~ Significant, gene_set_data, median)
  
  sigsig_median = as.numeric(medians[which(medians[,"Significant"] == "Sig-sig"), ][2])
  signon_median = as.numeric(medians[which(medians[,"Significant"] == "Sig-non"), ][2])
  nonnon_median = as.numeric(medians[which(medians[,"Significant"] == "Non-non"), ][2])
  
  means <- aggregate(Weight ~ Significant, gene_set_data, mean)
  
  sigsig_mean = as.numeric(means[which(means[,"Significant"] == "Sig-sig"), ][2])
  signon_mean = as.numeric(means[which(means[,"Significant"] == "Sig-non"), ][2])
  nonnon_mean = as.numeric(means[which(means[,"Significant"] == "Non-non"), ][2])
  
    
  sds <- aggregate(Weight ~ Significant, gene_set_data, sd)

  sigsig_sd = as.numeric(sds[which(sds[,"Significant"] == "Sig-sig"), ][2])
  signon_sd = as.numeric(sds[which(sds[,"Significant"] == "Sig-non"), ][2])
  nonnon_sd = as.numeric(sds[which(sds[,"Significant"] == "Non-non"), ][2])
  
  characteristics <- list(gene_pair_count, sigsig_count, sigsig_mean, sigsig_median, sigsig_sd, signon_count, signon_mean, signon_median, signon_sd, nonnon_count, nonnon_mean, nonnon_median, nonnon_sd)
  characteristic_file[gene_set_row_num, 3:15] = characteristics
  
}


write.csv(characteristic_file, paste(device_directory,"/clusters/curvy/complete_plots/characteristics.csv", sep=""))

####### RECYCLE BIN ########
# sigsig_data <- gene_set_data[which(gene_set_data[,"Significant"] == "Sig-sig"), ]
# 
# #for Sig-Sig  
# medians <- aggregate(Weight ~ Significant, sigsig_data, median)
# 
# means <- aggregate(Weight ~ Significant, sigsig_data, mean)
# 
# stddevs <- aggregate(Weight ~ Significant, sigsig_data, sd)
# 
# 
# #for Non-Sig
# signon_data <- gene_set_data[which(gene_set_data[,"Significant"] == "Sig-non"), ]
# 
# #for Non-Non
# nonnon_data <- gene_set_data[which(gene_set_data[,"Significant"] == "Non-non"), ]
