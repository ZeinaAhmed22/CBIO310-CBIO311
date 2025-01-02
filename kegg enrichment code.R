
#4.2. Gene Set Enrichment usig KEGG ----

#Prepare Input
#no need to repeat the following steps if done before (from reading: name the vetor)
# reading in data from deseq2
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
df_up <- read_table("up_res_df.txt")
View(df_up)
names(df_up)[names(df_up) == "X"] <- "GeneSymbol"

df_down <- read_table("down_res_df.txt")
View(df_down)
names(df_down)[names(df_down) == "X"] <- "GeneSymbol"

# we want the log2 fold change 
original_gene_list_up <- df_up$log2FC_1
original_gene_list_down <- df_down$log2FC_2
# name the vector
names(original_gene_list_up) <- df_up$Gene
names(original_gene_list_down) <- df_down$Gene

kegg_gene_list_up = df_up$log2FC_1
kegg_gene_list_down = df_down$log2FC_2

# Name vector with ENTREZ ids
names(kegg_gene_list_up) = df_up$Gene
names(kegg_gene_list_down) = df_down$Gene

# omit any NA values 
kegg_gene_list_up = na.omit(kegg_gene_list_up)
kegg_gene_list_down = na.omit(kegg_gene_list_down)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list_up = sort(kegg_gene_list_up, decreasing = TRUE)
kegg_gene_list_down = sort(kegg_gene_list_down, decreasing = TRUE)

# Set the desired KEGG organism
kegg_organism = "hsa"
# for reproducibility
set.seed(50) 
# Perform Gene Set Enrichment Analysis using KEGG
library(clusterProfiler)
kegg1 = gseKEGG(geneList     = kegg_gene_list_up,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               by = "fgsea",
               eps = 1e-10,
               keyType = "ncbi-geneid")

# Display the first few rows of the KEGG Gene Set Enrichment Analysis results
head(kegg1)

kegg2 = gseKEGG(geneList     = kegg_gene_list_down,
                organism     = kegg_organism,
                minGSSize    = 3,
                maxGSSize    = 800,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                by = "fgsea",
                eps = 1e-10,
                keyType = "ncbi-geneid")

# Display the first few rows of the KEGG Gene Set Enrichment Analysis results
head(kegg2)

dotplot(kegg1, showCategory = 10, title = "Enriched Pathways", split = ".sign") + facet_grid(. ~ .sign)
dotplot(kegg2, showCategory = 10, title = "Enriched Pathways", split = ".sign") + facet_grid(. ~ .sign)

library(enrichplot)

# Generate term similarity matrix
kegg_sim_1 <- pairwise_termsim(kegg1)
kegg_sim_2 <- pairwise_termsim(kegg2)

# Create the enrichment map plot
emapplot(kegg_sim_1)
emapplot(kegg_sim_2)

cnetplot(kegg1, categorySize = "pvalue", foldChange = kegg_gene_list_up)
cnetplot(kegg2, categorySize = "pvalue", foldChange = kegg_gene_list_down)

library(ggridges)
library(ggplot2)
ridgeplot(kegg1) + labs(x = "enrichment distribution")
gseaplot(kegg1, by = "all", title = kegg1$Description[1], geneSetID = 1)
ridgeplot(kegg2) + labs(x = "enrichment distribution")
gseaplot(kegg2, by = "all", title = kegg2$Description[1], geneSetID = 1)


# Step 4: Visualize a specific KEGG pathway
pathway_id <- "hsa04130"  # Example pathway ID for Homo sapiens
BiocManager::install("pathview")
library(pathview)
pathview(gene.data = kegg_gene_list_up, pathway.id = pathway_id, species = kegg_organism)
pathview(gene.data = kegg_gene_list_down, pathway.id = pathway_id, species = kegg_organism)

# Save enriched pathways
write.csv(as.data.frame(kegg1), file = "KEGG_GSEA_up_Results.csv")
write.csv(as.data.frame(kegg2), file = "KEGG_GSEA_down_Results.csv")


# Extract the pathway descriptions and save them to a CSV file
pathways_kegg_1 <- as.data.frame(kegg1@result[["Description"]])
write.csv(as.matrix(pathways_kegg_1),file="pathways_kegg_1.csv",quote = F,row.names = T)
pathways_kegg_2 <- as.data.frame(kegg2@result[["Description"]])
write.csv(as.matrix(pathways_kegg_2),file="pathways_kegg_2.csv",quote = F,row.names = T)

