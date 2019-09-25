library(dplyr)
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
organs <- lapply(strsplit(names(dataSet), "_"), "[[", 1) %>% unlist() %>% unique()
organs <- organs[which(organs != "gene")]
# or use this code 
# sapply(strsplit(names(dataSet[-1]), "_"), get_name) %>% unique()



# This function get gene_count compute the number of expressed genes for a given organ
# It thene return the a dataframe containing the number of lincRNA, the number of coding_protein
# and the total number of expressed genes

get_gene_count <- function(organ){
  dataSet %>% select("gene_type",starts_with(organ)) %>% 
    mutate(check=rowSums(.[,-1]>0.1)==dim(.[,-1])[2])  %>%
   filter(check == TRUE) %>% count(gene_type) %>% 
    as.data.frame() %>% select(gene_type,n) %>% add_row(gene_type="total", n=sum(.[,2]))
}

# We then apply the get_gene_count function to all the organs and then bind the columns of the data 
# frame returned by the get_gene_cout function for each organ. We then remove duplicated gene_type
# columns to finally get the reslut
  
geneCount <- lapply(organs, get_gene_count) %>% bind_cols() %>% select(-matches("gene_type[0-9]"))
names(geneCount) <- c("gene_type",organs)
view(geneCount)

# Venn diagram of expressed_gene (coding_protein_lincrna) for organs kidney, heart and  Adrenal_gland

# The get_expressed_genes function take an organ and returns the a dataframe containing the gene_type
# and the names of the genes that are expressed in th the organ

get_expressed_genes <- function(organ){
  dataSet %>% select("gene_type",starts_with(organ)) %>% 
    mutate(check=rowSums(.[,-1]>0.1)==dim(.[,-1])[2], gene_names=rownames(.))  %>%
    filter(check == TRUE) %>% select(gene_type, gene_names)
}

# Using the function above we can get the expressed genes for the three choosen organs 
heart_exp_genes <- get_expressed_genes("Heart")
kidney_exp_genes <- get_expressed_genes("Kidney")
adrenal_exp_genes <- get_expressed_genes("Adrenal")

venn_titles <- c("Heart", "kidney", "Adrenal")

# To draw the venn diagram we passe to the function venn.diagram a list containing expressed genes for 
# ech organ
venn.diagram(x = list(heart_exp_genes$gene_names, kidney_exp_genes$gene_names, adrenal_exp_genes$gene_names), 
             category.names = venn_titles,
             imagetype = "png",
             filename = 'Venn_diagram_1.png', 
             col = c("red", "blue", "green"))

venn.diagram(x = list(heart_exp_genes %>% filter(gene_type=='protein_coding') %>% select(gene_names) %>% unlist(),
                      kidney_exp_genes %>% filter(gene_type=='protein_coding') %>% select(gene_names) %>% unlist(), 
                      adrenal_exp_genes %>% filter(gene_type=='protein_coding') %>% select(gene_names) %>% unlist()), 
             category.names = venn_titles,
             imagetype = "png",
             filename = 'Venn_diagram_2.png', 
             col = c("red", "blue", "green"))

venn.diagram(x = list(heart_exp_genes %>% filter(gene_type=='lincRNA') %>% select(gene_names) %>% unlist(),
                      kidney_exp_genes %>% filter(gene_type=='lincRNA') %>% select(gene_names) %>% unlist(), 
                      adrenal_exp_genes %>% filter(gene_type=='lincRNA') %>% select(gene_names) %>% unlist()), 
             category.names = venn_titles,
             imagetype = "png",
             filename = 'Venn_diagram_3.png', 
             col = c("red", "blue", "green"))

# Barplot for liver 20 most expressed genes

# This code below takes the liver organ and compute 20 most expressed genes in all 4 samples

liver_20_most_exp_genes <- dataSet %>% select("gene_type",starts_with("Liver")) %>% 
      mutate(check=rowSums(.[,-1]>0.1)==dim(.[,-1])[2], gene_names=rownames(.))  %>%
      filter(check == TRUE) %>% mutate(row_means= rowMeans(.[,-c(1, dim(.)[2])])) %>%
      arrange(desc(row_means)) %>% head(n=20) %>%
      select(starts_with("Liver"), gene_names)

liver_20_most_exp_genes

# Drawing the barplot for each sample of data for the liver organ

ggplot(liver_20_most_exp_genes, aes(gene_names, Liver_Female_1)) +
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

ggplot(liver_20_most_exp_genes, aes(gene_names, Liver_Female_2)) + 
  geom_bar(stat = "identity",  fill="red") + theme(axis.text.x = element_text(angle = 90))

ggplot(liver_20_most_exp_genes, aes(gene_names, Liver_Male_1)) + 
  geom_bar(stat = "identity", fill="darkblue") + theme(axis.text.x = element_text(angle = 90))

ggplot(liver_20_most_exp_genes, aes(gene_names, Liver_Male_2)) + 
  geom_bar(stat = "identity", fill="darkorange") + theme(axis.text.x = element_text(angle = 90))
