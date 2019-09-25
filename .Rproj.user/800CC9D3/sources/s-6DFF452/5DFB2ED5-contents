library(VennDiagram)
library(stringr)
library(ggplot2)

## Load the data into R
dataSet <- read.table("./data.txt", sep="\t", stringsAsFactors = FALSE)
#View(dataSet)
head(dataSet)

# replace thymys to thymus
names(dataSet) <- str_replace(names(dataSet), pattern = 'Thymys',replacement = 'Thymus')

#find all the organs  
organs <- unique((sapply(strsplit(names(dataSet), "_"), function(x)x[[1]])))

sex_specific_organ <- str_starts(names(dataSet[-1]), pattern = "Vescicular|Uterus|Ovary|Testis")
sex_specific_data <- dataSet[,c(TRUE, sex_specific_organ)]
non_sex_specific_data <- dataSet[,c(TRUE, !sex_specific_organ)]

get_gene_count <- function(organ,data){
  all_samples <- data[,str_starts(names(data), pattern = paste("gene", organ, sep = "|"))] # select sample of an organ and the gene_type column
  gene_expressed_row_sum <- rowSums(all_samples[,-1]>0.1) # for each row we take every sample check if the value is greate than 0.1 and make the sum of boolean values for each row
  gene_expressed_check <- gene_expressed_row_sum == dim(all_samples[,-1])[2] # check if the gene is expressed by taking the row sum and see if it is equal to the number of samples for the organ
  gene_expressed <- all_samples[gene_expressed_check,] # use the codition to filter all samples and get only expressed genes
  gene_expressed <- table(gene_expressed$gene_type) # we use the table function to count the number of lincRNA and coding_protein in expressed genes
  gene_count <- as.data.frame(gene_expressed) # we get the gene type count as a data frame
  names(gene_count) <- c("gene_type",organ) # we the data frame the appropriate names 
  gene_count 
}

# Get the count of expressed genes for each organs

Adrenal <- get_gene_count("Adrenal",non_sex_specific_data)
Brain <- get_gene_count("Brain", non_sex_specific_data)
Forestomach <- get_gene_count("Forestomach",non_sex_specific_data)
Heart <- get_gene_count("Heart",non_sex_specific_data)
Kidney <- get_gene_count("Kidney",non_sex_specific_data)
Liver <- get_gene_count("Liver",non_sex_specific_data)
Large <- get_gene_count("Large",non_sex_specific_data)
Lung <- get_gene_count("Lung",non_sex_specific_data)
Muscle <- get_gene_count("Muscle",non_sex_specific_data)
Small <- get_gene_count("Small",non_sex_specific_data)
Spleen <- get_gene_count("Spleen",non_sex_specific_data)
Stomach <- get_gene_count("Stomach",non_sex_specific_data)
Thymus <- get_gene_count("Thymus",non_sex_specific_data)

# bind all the column of non sex specific organs but select gene_type column only for the first organ to avoid duplication of the gene_type column
sex_exp_gene_count <- cbind(Adrenal,Brain[2],Forestomach[2],Heart[2],Kidney[2],Liver[2],Large[2],Lung[2],
                        Muscle[2],Small[2],Spleen[2],Stomach[2])
sex_exp_gene_count


Vescicular <- get_gene_count("Vescicular",sex_specific_data)
Uterus <- get_gene_count("Uterus",sex_specific_data)
Ovary <- get_gene_count("Ovary",sex_specific_data)
Testis <- get_gene_count("Testis",sex_specific_data)

# bind all the column sex specific organs but select gene_type column only for the first organ to avoid duplication of the gene_type column
non_sex_exp_gene_count <- cbind(Ovary,Thymus[2],Uterus[2],Vescicular[2])
non_sex_exp_gene_count



#  Venn Diagramm for the three organs : kidney, heart and  Adrenal_gland

get_expressed_genes <- function(organ){
  all_samples <- dataSet[,str_starts(names(dataSet), pattern = paste("gene", organ, sep = "|"))] # select sample of an organ and the gene_type column
  gene_expressed_row_sum <- rowSums(all_samples[,-1]>0.1) # for each row we take every sample check if the value is greate than 0.1 and make the sum of boolean values for each row
  gene_expressed_check <- gene_expressed_row_sum == dim(all_samples[,-1])[2] # check if the gene is expressed by taking the row sum and see if it is equal to the number of samples for the organ
  gene_expressed <- all_samples[gene_expressed_check,] # use the codition to filter all samples and get only expressed genes
  gene_expressed$gene_name <- rownames(gene_expressed) # we store the names of the genes in a new column named gene_names
  gene_expressed <- gene_expressed[c("gene_type", "gene_name")]  # We select the gene_type and the gene_name column
  gene_expressed 
}



kidney_genes <- get_expressed_genes("Kidney") # get genes expressed in kidney
heart_genes <- get_expressed_genes("Heart") # get genes expressed in Heart
adrenal_genes <- get_expressed_genes("Adrenal") # get genes expressed in Adrenal

# ploting the venn diagrams

venn_titles <- c("Heart", "kidney", "Adrenal")

venn.diagram(x = list(heart_genes$gene_name, kidney_genes$gene_name, adrenal_genes$gene_name), 
             category.names = venn_titles,
             imagetype = "png",
             filename = 'Venn_diagram_1.png', 
             col = c("red", "blue", "green"))

venn.diagram(x = list(heart_genes[heart_genes$gene_type=="protein_coding","gene_name"], 
                      kidney_genes[kidney_genes$gene_type=="protein_coding","gene_name"], 
                      adrenal_genes[adrenal_genes$gene_type=="protein_coding","gene_name"]), 
             category.names = venn_titles,
             imagetype = "png",
             filename = 'Venn_diagram_2.png', 
             col = c("red", "blue", "green"))

venn.diagram(x = list(heart_genes[heart_genes$gene_type=="lincRNA","gene_name"], 
                      kidney_genes[kidney_genes$gene_type=="lincRNA","gene_name"], 
                      adrenal_genes[adrenal_genes$gene_type=="lincRNA","gene_name"]), 
             category.names = venn_titles,
             imagetype = "png",
             filename = 'Venn_diagram_3.png', 
             col = c("red", "blue", "green"))

# Make the barplot


  liver_samples <- dataSet[,str_starts(names(dataSet), pattern ="Liver")] # select sample of the liver organ and the gene_type column
  gene_expressed_row_sum <- rowSums(liver_samples>0.1) # for each row we take every sample check if the value is greate than 0.1 and make the sum of boolean values for each row
  gene_expressed_check <- gene_expressed_row_sum == dim(liver_samples)[2] # check if the gene is expressed by taking the row sum and see if it is equal to the number of samples for the organ
  liver_genes <- liver_samples[gene_expressed_check,] # we filter the data and select expressed genes 
  liver_genes$row_mean <- rowMeans(liver_genes) # we conpute the row mean using the for samples 
  liver_genes_ordered <- liver_genes[order(liver_genes$row_mean, decreasing = TRUE),] # we order the expressed data according to the row_mean
  liver_20_genes <- head(liver_genes_ordered, n=20) # select the 20 most expressed genes 
  
  
# Drawing the barplot for each sample of data for the liver organ and also using the rownmae as the label for the x axis
  
ggplot(liver_20_genes, aes(rownames(liver_20_genes), Liver_Female_1)) +
    labs(title = "Liver female 1 20 most expressed genes", x="genes names", y="gene expression level") +
    geom_bar(stat = "identity",fill="pink") + theme(axis.text.x = element_text(angle = 90))
  
ggplot(liver_20_genes, aes(rownames(liver_20_genes), Liver_Female_2)) + 
    labs(title = "Liver female 2 20 most expressed genes", x="genes names", y="gene expression level") +
    geom_bar(stat = "identity", fill="red") + theme(axis.text.x = element_text(angle = 90))
  
ggplot(liver_20_genes, aes(rownames(liver_20_genes), Liver_Male_1)) + 
    labs(title = "Liver male 1 20 most expressed genes", x="genes names", y="gene expression level") +
    geom_bar(stat = "identity", fill="lightgreen") + theme(axis.text.x = element_text(angle = 90))
  
ggplot(liver_20_genes, aes(rownames(liver_20_genes), Liver_Male_2)) + 
    labs(title = "Liver male 1 20 most expressed genes", x="genes names", y="gene expression level") +
    geom_bar(stat = "identity", fill="darkorange") + theme(axis.text.x = element_text(angle = 90))
