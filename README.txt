GSE10072_prepro_nonredun_var_noProbe_z.txt 	= Preprocessed_expression_input
NewPathwayAPI_forGCN.txt 			= Gene-set file for constructing GCN network of each gene-set (first column is gene-set ID)

NewPathwayAPI_withName.txt 			= Gene-set file with biological process name of each gene-set ID
All_sig_gene_10072.txt 				= Gene-set file with all the gene member is significant genes based-on expression value from GSE10072 dataset

genesetCoNet_vis.r 	= R code for constructing GCN of each gene-set 
buildGraphGCNweight.r 	= R code for preparation file by specify gene-pair (sig-sig, sig-non, non-non) and draw graph