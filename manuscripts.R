###### differential gene analysis was performed on each dataset separately
library(limma)
fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
library(dplyr)
diffgene <- allDiff %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)

######## integrated differential gene
########
library(RobustRankAggreg)
ups=aggregateRanks(glist = glist, N = length(unique(unlist(glist))))
tmp=as.data.frame(table(unlist(glist)))
ups$Freq=tmp[match(ups$Name,tmp[,1]),2]
gs=ups[ups$Score < 0.05,1]

##############
## STRING web tool and Cytoscape software
##
library(UpSetR) 

####### logistic regression analysis
library(rms) 
fit <- glm(status~FN1+CD44+C1QB+C1QA, data=bab3, family = binomial(link="logit"), x=T)

################# Correlation analysis
dd = cor.test(genedata,as.numeric(data[,i]),method="pearson")

########################## GO enrichment analysis
library(clusterProfiler)
go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")

######### KEGG pathway function enrichment analyses
KEGG <- enrichKEGG(gene = gene$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   use_internal_data =T)

###GSEA
kegg <- read.gmt("data/c2.cp.kegg.v7.4.symbols.gmt")
y <- GSEA(geneList,TERM2GENE =kegg)
##########
