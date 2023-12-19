bulk <- read.table("data_RNA_Seq_v2_expression_median.txt", header = T, sep = '\t')
bulk$sum <- rowSums(bulk[,3:ncol(bulk)])
bulk <- bulk %>% arrange(desc(sum))
bulk <- bulk[!duplicated(bulk$Hugo_Symbol),]
bulk <- bulk[!is.na(bulk$Hugo_Symbol),]
rownames(bulk) <- bulk$Hugo_Symbol
bulk <- bulk[,c(-1,-2,-ncol(bulk))]
bulk <- t(bulk)
bulk <- as.data.frame(bulk)

clinical <- read.csv("tcga_chol_io_v22.8.2_u160F.csv")
clinical$Sample_ID <- str_replace_all(clinical$Sample_ID, "-",".")
clinical <- clinical %>% dplyr::select(Sample_ID, immune_phenotype)
clinical$immune_phenotype <- ifelse(clinical$immune_phenotype=="Inflamed","Inflamed", "Others")

bulk <- merge( clinical, bulk, by.y="row.names", by.x="Sample_ID")
rownames(bulk) <- bulk$Sample_ID
bulk <- bulk[,-1]
res <- data.frame(matrix(ncol = 3, nrow = 1, NA))

inflamed <- apply(bulk[bulk$immune_phenotype=="Inflamed", -1],2, mean)
Others <- apply(bulk[bulk$immune_phenotype=="Others", -1], 2,mean)
inflamed <- data.frame(gene=names(inflamed), inflamed=inflamed)
Others <- data.frame(gene=names(Others), Others=Others)

res <- merge(inflamed, Others, by="gene")
res$log2FC <- log2(res$inflamed/res$Others)
res <- res %>% dplyr::arrange(desc(log2FC))
res <- res[!res$log2FC %in% c(0, Inf, -Inf),]

GSEA_symbol_list <- res$log2FC
names(GSEA_symbol_list) <- res$gene
GSEA_symbol_list <- sort(GSEA_symbol_list, decreasing = T)

gse_GO <-gseGO(geneList = GSEA_symbol_list,# ont = "BP", maxGSSize = 10000,
                            keyType = "SYMBOL", #nPerm=10000,
                            pvalueCutoff = 0.05, verbose = T,
                            OrgDb = get("org.Hs.eg.db"), pAdjustMethod = "fdr")
