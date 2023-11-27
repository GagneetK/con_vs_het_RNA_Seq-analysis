fly = read.csv("/Users/gagneetkaur/Desktop/final_project/modified.txt",header =T, sep = "\t", row.names = 1)
head(fly)

orthologs=read.csv("/Users/gagneetkaur/Downloads/dmel-orthologs-in-drosophila-species-fb-2021-02 2.tsv", header=T, sep="\t")
orthosearch=function(x){orthologs[orthologs$Ortholog_FBgn_ID==x,1]}
orthosearchsingle=function(x){orthosearch(x)[1]}
head(orthologs)

orthosearchinverse=function(x){orthologs[orthologs$FBgn_ID ==x,6]}

## Load databases
library(edgeR) 
library(tidyverse) 
library(Glimma)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ggrepel)



colnames(fly) %>% str_subset(pattern = "C|H") -> newCOls
fly_mating <- dplyr::select(fly, all_of(newCOls))

colnames(fly_mating) %>% str_replace("[:digit:]*_[:digit:]", "")-> groups
groups

#colnames(fly_mating) %>% str_replace("_[:digit:]", "")-> groups2
#groups2

#?substr

y <- DGEList(counts = fly_mating, group = groups)



keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]


y <- calcNormFactors(y)

design <- model.matrix(~ 0 + groups)
colnames(design) %>%  str_replace("groups", "") -> colnames(design)
y <- estimateDisp(y, design)


data = log2(y$counts+1)
data = as.matrix(data) # required for base R functions
sample_cor = cor(data)
pheatmap(sample_cor)



glMDSPlot(y, groups = y$samples$group)


plotBCV(y)


fit <- glmQLFit(y, design)


C.v.H <- makeContrasts(C-H, levels=design)
#fH.v.tm <- makeContrasts(H6-H12, levels=design)



qlf.C.v.H  <- glmQLFTest(fit, contrast = C.v.H)
qlf.C.v.H.tTags <- topTags(qlf.C.v.H, n = NULL) # Need to speicy n as NULL tp get all genes
qlf.C.v.H.tTags.table <- qlf.C.v.H.tTags$table



### Check the number of DE genes:
#summary(decideTests(qlf.fH.v.cm))
summary(decideTests(qlf.C.v.H))



### get rid of row names:
#qlf.fH.v.cm.tTags.table <- qlf.fH.v.cm.tTags.table %>% rownames_to_column(var = "gene") %>% as_tibble()

qlf.C.v.H.tTags.table <- qlf.C.v.H.tTags.table %>% rownames_to_column(var = "gene") %>% as_tibble()



#qlf.fH.v.cm.tTags.table %>% 
#mutate(significant = if_else((logFC > 1 & FDR < 0.05) | (logFC < -1 & FDR < 0.05), #"yes", "no")) -> qlf.fH.v.cm.tTags.table

qlf.C.v.H.tTags.table %>% 
  mutate(significant = if_else((logFC > 1 & FDR < 0.05) | (logFC < -1 & FDR < 0.05), "yes", "no")) -> qlf.C.v.H.tTags.table



qlf.C.v.H.tTags.table %>% 
  ggplot(aes(x = logFC, y = -log10(PValue), colour = significant)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(data = filter(qlf.C.v.H.tTags.table, (-log10(PValue) > 5) & (significant=="yes")), aes(label = gene), size = 2) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "red"))



cpm = cpm(fly)
cpm=data.frame(cpm)
cpm_mating <- dplyr::select(cpm, all_of(newCOls))

cpm_mating2 <- cpm_mating %>% rownames_to_column(var = "gene") %>% as_tibble()

cpm_mating2 %>% gather(key = sample, value = cpm, -gene) %>% separate(sample, c("sample","replicate"))%>% separate(sample, c("condition", "time"), sep=1) -> cpm_mating_gathered




### For heatmaps, we often have to perform some transformations of the CPM values

# make a vector of significant genes
qlf.C.v.H.tTags.table %>%  filter(significant == "yes") %>% pull(gene) -> sigGenes_C.v.H
# select significant genes only
data = subset(cpm_mating, rownames(cpm_mating) %in% sigGenes_C.v.H)
# log transform the TPM values
data = log2(data + 1)
# Calculate the median centered TPM value
data = as.data.frame(t(scale(t(data), scale = F)))
# limit the maximum/minimum log2 value displayed to 4
data[data < -2] = -2
data[data > 2] = 2


pheatmap(mat = data, 
         border_color = NA, 
         show_colnames = TRUE, 
         show_rownames = F, 
         drop_levels = TRUE, 
         annotation_names_row = F, 
         fontsize = 8)


plotGene <- function(id) {
  cpm_mating_gathered %>% 
    filter(gene == id) %>% 
    group_by(gene, condition,time) %>% 
    summarize(mean = mean(cpm), n = n(), se = sd(cpm)/sqrt(n)) %>% 
    ggplot(aes(fill = time)) +
    geom_col(aes(condition, mean), position = "dodge") +
    geom_errorbar(aes(x = condition, ymin = mean - se, ymax = mean + se), position = position_dodge(0.9), width = 0.3) +
    geom_point(data = filter(cpm_mating_gathered, gene == id), aes(condition, cpm), position = position_dodge2(width = 0.5)) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> myPlot ### must make a plot object for a function to work, then use the "return command"
  return(myPlot)
}



orthosearch("FBgn0209590")##defensin
plotGene("FBgn0209590")

orthosearch("FBgn0201305")## noannotation

orthosearch("FBgn0197572")#socs36e - involved in JAK-STAT Pathway
plotGene("FBgn0197572")

orthosearch("FBgn0197687")

plotGene("FBgn0210579")##Relish




## Create a vector of gene names for enrichment analysis:
as_tibble(qlf.C.v.H.tTags.table) %>% filter(significant == "yes") %>% pull(gene) -> as.data.frame(sigGenes_C.v.H)
as.data.frame(sigGenes_C.v.H) ->sigGenes_C.v.H

colnames(sigGenes_C.v.H) <- "significant genes"
print(sigGenes_C.v.H)

### For GO enrichment analysis,
jointdataset <-merge (sigGenes_C.v.H, orthologs, by.x = 'significant genes', by.y ='Ortholog_FBgn_ID')
joint <- (jointdataset[2])
print(joint)


orthosearch=function(x){orthologs[orthologs$Ortholog_FBgn_ID==x,1]}
orthosearchsingle=function(x){orthosearch(x)[1]}
orthosearchinverse=function(x){orthologs[orthologs$FBgn_ID ==x,6]}




sigGenes_C.v.H_mel=lapply(sigGenes_C.v.H, orthosearchsingle)
sigGenes_C.v.H_mel= sigGenes_C.v.H_mel[!is.na(sigGenes_C.v.H_mel)]
sigGenes_C.v.H_mel= unlist(sigGenes_C.v.H_mel)
print(sigGenes_C.v.H_mel)

allGenes_mel=lapply(rownames(fly_mating), orthosearchsingle)
allGenes_mel=allGenes_mel[!is.na(allGenes_mel)]
allGenes_mel= unlist(allGenes_mel)
print(allGenes_mel)

orthologs[1]


### For GO enrichment analysis, we need to extract the test set and the background set
mapIds(org.Dm.eg.db, keys = sigGenes_C.v.H_mel, column = "ENTREZID", keytype = "FLYBASE") -> sigGenes_entrez
mapIds(org.Dm.eg.db, keys = allGenes_mel, column = "ENTREZID", keytype = "FLYBASE") -> allGenes_entrez


ego <- enrichGO(gene          = sigGenes_entrez,
                universe      = allGenes_entrez,
                keyType       = 'ENTREZID',
                OrgDb         = org.Dm.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, 
                readable      = TRUE)

dotplot(ego)

ego %>% as_tibble()

ego %>% as_tibble() %>%
  dplyr::filter(p.adjust<3, Count>2) %>%
  separate(GeneRatio, c("geneDE", "background")) %>%
  mutate(GeneRatio= as.numeric(geneDE)/as.numeric(background))%>%
  ggplot( aes(x=GeneRatio, y=fct_reorder(Description, Count)))+
  geom_point(aes(colour=-log(p.adjust), size=Count))+
  facet_grid(rows= vars(ONTOLOGY), scale="free", space ="free")+
  ylab(NULL)


