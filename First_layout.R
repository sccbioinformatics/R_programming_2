library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridarrange)

pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k")
# write a dense matrix
write.table(pbmc.data,"PBMC_full_matrix.txt")

fd <- as.matrix(read.delim("../PBMC_full_matrix.txt",header=T,row.names=1,sep="\t"))

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

features <- apply(fd >0, 2, sum)
nCount <- apply(fd,2,sum)
perc.mt <- 100*apply(fd[grep("^MT-",rownames(fd)),],2,sum)/apply(fd,2,sum)


min.cls <- apply(fd >0, 1, sum)

qc.df <- data.frame(features=features,nCount=nCount,perc.mt=perc.mt)

qc.df.mlt <- melt(qc.df)

qc.feat <- data.frame(features=features)
qc.nCount <- data.frame(nCount=nCount)
qc.perc.mt <- data.frame(perc.mt=perc.mt)

p <- ggplot(melt(qc.feat), aes(variable,y=value))
p <- p + geom_violin() +  geom_jitter(shape=".", position=position_jitter(0.2)) + facet_grid(rows = vars(variable),scales="free")

q <- ggplot(melt(qc.nCount), aes(variable,y=value))
q <- q + geom_violin() +  geom_jitter(shape=".", position=position_jitter(0.2)) + facet_grid(rows = vars(variable),scales="free")

r <- ggplot(melt(qc.perc.mt), aes(variable,y=value))
r <- r + geom_violin() +  geom_jitter(shape=".", position=position_jitter(0.2)) + facet_grid(rows = vars(variable),scales="free")

####


fill_color <- "lightcoral"

# Assuming your data frame is called qc.df
# Create separate violin plots for each variable
plots <- lapply(names(qc.df), function(variable) {
  ggplot(qc.df, aes(x = 1, y = qc.df[[variable]], fill = fill_color)) +
    geom_violin() +
    labs(x = NULL, y = variable) +
    scale_fill_manual(values = fill_color, name = variable) +  # Set the fill color
    theme_minimal() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    geom_jitter(shape=".", position=position_jitter(0.2))
})

# Arrange the plots in a single row
grid.arrange(grobs = plots, ncol = length(plots))

p1 <- ggplot(qc.df, aes(x=nCount, y=perc.mt)) + geom_point(color = fill_color)
p2 <- ggplot(qc.df, aes(x=nCount, y=features)) + geom_point(color = fill_color)
p1|p2
## make the same plots as seurat scatters

#subset the data.frame:
d.filt <- qc.df %>% filter(features  <= 2500) %>% filter(features  >= 200) %>% filter(perc.mt <= 5)

count.filt <- fd[names(which(min.cls>2)),rownames(d.filt)]

#make a list


#######
# Use the package


library(SCAP)

pbmc.data <- readRDS("/home/shamit/Git/R_programming_2/PBMC_data.rds")
pbmc <- CreateMySCO(as.matrix(pbmc.data))
pbmc <- CalcMitoPct(pbmc,"^MT-")
#MakeQCPlots(pbmc)
pbmc <- FilterData(pbmc,sub="features > 200 & features < 2500 & perc.mt < 5",min.cells=3)
pbmc <- NormaliseData(pbmc,10000)
pbmc <- FindHVGs(pbmc,2000)
PlotHVGs(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- CalcPCs(pbmc)
pbmc <- ClusterCells(pbmc,nPC=10,nK=30)
pbmc <- MakeUMAP(pbmc,nPC=10)


dd <- pbmc@data.scale[pbmc@hvgs,]
dd.pr <- prcomp(t(dd),center=FALSE)$x


snn <- RANN::nn2(dd.pr[,1:10], k=30)$nn.idx
adjacency_matrix <- matrix(0L, nrow(dd.pr), nrow(dd.pr))

rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- colnames(dd)

for(ii in 1:nrow(dd.pr)) {
    adjacency_matrix[ii,rownames(dd.pr)[snn[ii,]]] <- 1L
}
#check that rows add to k
sum(adjacency_matrix[1,]) == 30
table(apply(adjacency_matrix, 1, sum))

clus <- leiden(adjacency_matrix)

plot(dd.pr[,1:2])


points(dd.pr[which(clus==1),1:2],col="red")


ump <- umap(dd.pr[,1:10])
#ump.gex <- umap(t(pbmc@data.scale[pbmc@hvgs,]))


plot(ump$layout)

tcols <- terrain.colors(max(clus))

points(ump$layout[which(clus==1),1:2],col="red")

for(i in 1:max(clus)){

  points(ump$layout[which(clus==i),1:2],col=tcols[i])

}
