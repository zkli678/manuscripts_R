###download the data from GEO database and raw data processing for each dataset
###https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE30528&format=file
###https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE47183&format=file
###https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE104948&format=file
###### differential gene analysis was performed on each dataset separately
library(limma)
fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
library(dplyr)
diffgene <- allDiff %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)
library(pheatmap)
library(ggplot2)
pheatmap(heatdata,
         cluster_rows = T,
         cluster_cols = F,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = F,
         scale = "row",
         color =colorRampPalette(c("green", "black","red"))(100),
         cellwidth = 15, 
         cellheight = 0.8,
         fontsize = 10
)

ggplot(data=data, aes(x=logFC, y =-log10(P.Value))) +
  geom_point(data=subset(data,abs(data$logFC) <= 1),aes(size=abs(logFC)),color="grey",alpha=0.1) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 1),aes(size=abs(logFC)),color="red",alpha=0.2) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -1),aes(size=abs(logFC)),color="green",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')
######## integrated differential gene
########
library(RobustRankAggreg)
ups=aggregateRanks(glist = glist, N = length(unique(unlist(glist))))
tmp=as.data.frame(table(unlist(glist)))
ups$Freq=tmp[match(ups$Name,tmp[,1]),2]
gs=ups[ups$Score < 0.05,1]

##############
## STRING web tool and Cytoscape software
####### logistic regression analysis
library(rms) 
fit <- glm(status~FN1+CD44+C1QB+C1QA, data=data, family = binomial(link="logit"), x=T)
################# Correlation analysis
library(dplyr)
gene <- "FN1"
genedata <- as.numeric(data[,gene])

for(i in 1:length(genelist)){
  print(i)
  dd = cor.test(genedata,as.numeric(data[,i]),method="pearson")
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("FN1","symbol","cor","p.value")


cor_data <- correlation %>% 
  
  filter(p.value < 0.05) %>% 
  
  arrange(desc(abs(cor)))%>% 
  
  dplyr::slice(1:50)
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
