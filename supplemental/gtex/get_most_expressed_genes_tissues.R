expression_data = read.csv("supplemental/gtex/AvgExpr.csv", stringsAsFactors = FALSE)

expression_data_subset = expression_data[, names(expression_data) %in% c("X", "Adipose...Subcutaneous",
                                                                    "Artery...Tibial",
                                                                    "Heart...Left.Ventricle",
                                                                    "Lung",
                                                                    "Muscle...Skeletal",
                                                                    "Nerve...Tibial",
                                                                    "Skin...Sun.Exposed..Lower.leg.",
                                                                    "Thyroid")]


# for(i in 2:dim(expression_data_subset)[2]){
#   temp = data.frame(gene = expression_data_subset[, 1], tissue = expression_data_subset[, i], stringsAsFactors = FALSE)
#   genes_temp = temp[order(temp[, 2], decreasing = TRUE), 1]
#   genes_temp = genes_temp[1:10]
#   genes = c(genes, genes_temp)
# }

temp = apply(as.matrix(expression_data_subset[, -1]), 1, min)
expression_data_subset = expression_data_subset[order(temp, decreasing = TRUE)[1:50], ]

runinfo = read.csv("supplemental/gtex/rna_seq_runinfo.csv", stringsAsFactors = FALSE)

runinfo_subset = runinfo[runinfo$body_site_s %in% c("Adipose - Subcutaneous",
                                                    "Artery - Tibial",
                                                    "Heart - Left Ventricle",
                                                    "Lung",
                                                    "Muscle - Skeletal",
                                                    "Nerve - Tibial",
                                                    "Skin - Sun Exposed (Lower leg)",
                                                    "Thyroid"), ]

known_genes = read.delim("supplemental/gtex/knownGene.txt.gz", header = FALSE, stringsAsFactors = FALSE)

gene_names = read.delim("supplemental/gtex/protein-coding_gene.txt", header = TRUE, stringsAsFactors = FALSE)

gene_names_ensembl = gsub("[.].*$", "", expression_data_subset[, 1])

gene_names_ucsc = gene_names$ucsc_id[gene_names$ensembl_gene_id %in% gene_names_ensembl]

gene_info = known_genes[known_genes$V1 %in% gene_names_ucsc, 1:5]

sample_list = runinfo_subset$Run_s
gene_region = paste0(gsub("^chr", "", gene_info$V2), ":", gene_info$V4, "-", gene_info$V5)

write(sample_list, "supplemental/gtex/sample_list.txt")
write(gene_region, "supplemental/gtex/gene_region.txt")
