###download the data from GEO database and raw data processing for each dataset
###https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE30528&format=file
###https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE47183&format=file
###https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE104948&format=file
library(oligo) 
dir='/GSE_RAW/'
setwd(dir)
celFiles <- list.celfiles(listGzipped = T)
celFiles
affyRaw <- read.celfiles( celFiles )
eset <- rma(affyRaw)
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
##
library(UpSetR) 
upset(data, sets =c("MCC","MNC","Degree","EPC",         
                    "Bottleneck","EcCentricity","Closeness" ,   "Radiality",   
                    "Betweenness","Stress" ), nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), 
      queries =list(list(query=intersects, params=list("MCC","MNC","Degree","EPC",         
                                                       "Bottleneck","EcCentricity","Closeness" ,   "Radiality",   
                                                       "Betweenness","Stress" ), color="red", active=TRUE)),
      decreasing = c(TRUE,FALSE))
####### logistic regression analysis
library(rms) 
fit <- glm(status~FN1+CD44+C1QB+C1QA, data=data, family = binomial(link="logit"), x=T)

################# Correlation analysis
library(dplyr)
dd = cor.test(genedata,as.numeric(data[,i]),method="pearson")
##
ggplot(bb,aes(COL6A3,FN1))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  stat_cor(method = "pearson",digits = 3,size=5)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))

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
