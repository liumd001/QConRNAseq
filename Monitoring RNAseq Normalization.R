
#####################################################
##  Data Description  GSE60424  
#####################################################
require(tidyverse)
require(openxlsx)
require(DESeq2)
require(org.Hs.eg.db)
require(pheatmap)
require(limma)
require(flextable)
library(org.Hs.eg.db)
require(pheatmap)
require(biomaRt)

recode = dplyr::recode
filter = dplyr::filter
rename = dplyr::rename
select = dplyr::select

dtafolder = "/Volumes/home/Drive/learning_lmd/R/deconvolution"
dtafile = file.path(dtafolder,"dataset","GSE60424_GEOSubmit_FC1to11_normalized_counts.txt.gz")
sampleAnnfile = file.path(dtafolder,"dataset","sampleAnnotation.xlsx")

dta = read.table(gzfile(dtafile), header = TRUE)
ann = read.xlsx(sampleAnnfile) %>%
  filter(!is.na(Donorid)) %>%
  mutate(Cellcount = as.numeric(Cellcount),
         Donorid = factor(Donorid))

#characteristic table 
ann %>% 
  select(Donorid, Gender, Diseasestatus, Age, Collectiondate, Race, Smoker) %>% 
  distinct() %>% 
  arrange(Diseasestatus) %>% 
  flextable() %>% 
  autofit()

#################################################
##    Housekeeping genes  
#################################################
hk = as.data.frame(rbind(c("ACTB","Cellular cytoskeleton"),
                         c("TBP","Transcription/Translation"),
                         c("RPL13A","Transcription/Translation"),
                         c("RPL27","Transcription/Translation"),
                         c("GAPDH","Transcription/Translation"),
                         c("SDHA","Transcription/Translation"),
                         c("HPRT1","Transcription/Translation"),
                         c("PGK1","Transcription/Translation"),
                         c("ATP5PB","Transcription/Translation"),
                         c("TBP","Transcription/Translation"),
                         c("UBC","Protein Degradation"),
                         c("B2M","MHC"),
                         c("RER1","Structural protein")))
names(hk) = c("hgnc_symbol","description")

mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mart = useDataset("hsapiens_gene_ensembl", mart)
annotLookup = getBM(
  mart = mart,
  attributes = c("hgnc_symbol","ensembl_gene_id"),
  uniqueRows = TRUE
)

hk = annotLookup[annotLookup$hgnc_symbol %in% hk$hgnc_symbol,]
write.csv(hk, file.path(dtafolder,"HKgenes.csv"), row.names = FALSE)

######################################################
##    Normal ranges of blood cell types  
######################################################
require(rvest)
require(DESeq2)
require(tidyverse)
require(limma)
require(ggplot2)
require(biomaRt)

url = "https://www.bio-rad-antibodies.com/flow-cytometry-cell-frequency.html"
page = read_html(url)
table = html_table(page, fill = TRUE)[[3]]
normal_range = data.frame(celltype = table[[5]][-1],
                          range = table[[6]][-1])  %>%
  dplyr::filter(celltype != "")
normal_range$low = as.numeric(do.call(rbind,strsplit(normal_range$range, split = "-"))[,1])
normal_range$high = as.numeric(do.call(rbind,strsplit(normal_range$range, split = "-"))[,2])
normal_range
write.csv(normal_range, file.path(dtafolder,"results","NormalRange_BloodComponents.csv"), row.names = FALSE)

####################################################
##    Samples  
####################################################
#total cell numbers
ann = ann %>%
  left_join(ann %>%
              group_by(Donorid) %>%
              summarise(totcell = sum(Cellcount, na.rm = TRUE)), by = "Donorid") %>%
  mutate(Cellfreq = Cellcount/totcell)


dta1 = ann %>%
  mutate(Cellcount = as.numeric(Cellcount),
         Donorid = factor(Donorid)) %>%
  filter(Celltype != "Whole Blood") 

g1 = dta1 %>%
  ggplot(aes(x = Celltype, y = Cellcount,color = Donorid, group = Donorid)) +
  geom_line() +
  theme_bw() +
  scale_y_log10() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3)) +
  labs(x  = "Cell Types", y = "Cell Counts", color = "Donor")

temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g1, width = 10, height = 8, units = "in")
file.copy(temp_file, file.path(dtafolder, "results","TotalCells_acrossSamples.png"))

####################################################
##    Sequencing Depth  
####################################################
#total reads
tot = colSums(dta[,-1])

dta2 = data.frame(Title = names(tot),
                  totseq = tot) %>%
  left_join(ann %>%
              select(Title, Donorid, Celltype), by = "Title")

# dta2 %>%
#   ggplot(aes(x = Celltype, y = tot,color = Donorid, group = Donorid)) +
#   geom_line() +
#   theme_bw() +
#   scale_y_log10() +
#   theme(legend.position = "bottom") +
#   guides(color = guide_legend(nrow = 3)) +
#   labs(x  = "Cell Types", y = "Cell Counts", color = "Donor")

dta3 = ann %>%
  dplyr::select(Title, Donorid, Celltype, totcell) %>%
  left_join(dta2 %>%
              dplyr::select(Title, Donorid, Celltype, totseq))

g2 = dta3 %>%
  ggplot(aes(x = totcell, y = totseq, color = Celltype)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Total Cells", y = "Total Read Counts", color = "Cell Types")

#reads by gene
totgene = data.frame(counts = rowSums(log2(dta[,-1] + 0.5)))
g3 = totgene %>%
  ggplot(aes(x = counts)) +
  geom_histogram() +
  theme_bw() +
  ylab("Total Counts By Gene") +
  xlab("Genes")


ft = ann %>%
  select(Donorid, Celltype, Age, Gender, Diseasestatus) %>%
  group_by(Donorid, Celltype, Age, Gender, Diseasestatus) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "Celltype", values_from = "n") %>%
  arrange(Diseasestatus) %>%
  flextable() %>%
  autofit() %>%
  save_as_image(path = file.path(dtafolder,"results","SampleCharacteristics.png"))

##############################################
##   Normalization  
##############################################
#log transform
#dta_raw = dta
#dta[,-1] = apply(dta[,-1],2, function(x) log2(x + 0.5))

#deviation
sum = rowSums(log2(dta[,-1] + 0.5), na.rm = TRUE)
med = apply(dta[,-1], 1, function(x) median(log2(x + 0.5), na.rm = TRUE))
sd = apply(dta[,-1],1,function(x) sd(log2(x + 0.5), na.rm = TRUE))

require(limma)
lowessfit = limma::weightedLowess(x = med, y = sd)

dd1 = data.frame(geneid = dta[,1],
                 sum = sum,
                 med = med,
                 sd = sd,
                 fittedsd = lowessfit$fitted,
                 geneid = dta[,1])

#med vs sd
g4 = dd1 %>%
  ggplot(aes(x = med, y = sd)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(x = med, y = fittedsd), color = "red") +
  geom_smooth(method = "loess", formula = y ~ x, na.rm = TRUE, color = "red")

g4

#why there are so many -1s in med
filteredgenes = dd1 %>%
  filter(med == -1) %>%
  pull(geneid)

t1 = dta %>%
  filter(genenames %in% filteredgenes)

#median vs sum
quantiles_med = quantile(dd1$med, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE)
quantiles_tot = quantile(dd1$sum, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE)

#filter
filter = dd1$med >= quantiles_med[8] & dd1$sum >= quantiles_tot[8]


#set filter at 70% percentile
thelabel = data.frame(x = c(1,10,2.5,10),
                      y = c(-100,1650,1500,350),
                      label = c(paste0("n = ",sum(!filter)),
                                paste0("n = ",sum(filter)),
                                "Quantile 80",
                                "Quantile 80"))
g5 = dd1 %>%
  ggplot(aes(x = med, y = sum)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "red") +
  geom_vline(xintercept = quantiles_med[8], color = "green") +
  geom_hline(yintercept = quantiles_tot[8], color = "purple") +
  theme_bw() +
  geom_label(data = thelabel, aes(x = x, y = y, label = label), col = "blue") +
  labs(x = "Median Counts", y = "Total Counts") 
g5

require(cowplot)
g6 = plot_grid(g1, g2, g3, g5, ncol = 2)

temp_file <- tempfile(fileext = ".png")
save_plot(temp_file, g6, base_height = 6, base_width = 8)
file.copy(temp_file, file.path(dtafolder, "results","CompositePlot_DataPreP.png"), overwrite = TRUE)

######################################################
##    Normalization using DESeq2  
######################################################
require(DESeq2)
dta5 = apply(as.data.frame(dta[filter,-1]),2,as.numeric)
rownames(dta5) = dta[filter,1]

#dta5 = apply(dta5, 2, as.numeric)
se =  SummarizedExperiment(dta5)
rownames(ann) = ann$Title
colData(se) <- DataFrame(ann)

dds = DESeqDataSet(se, design = ~ 1)

#size factors
dds = DESeq(dds)

data.frame(Title = names(sizeFactors(dds)),
           sizefactor = sizeFactors(dds),
           colsum = colSums(counts(dds))) %>%
  left_join(ann %>%
              select(Title, Celltype)) %>%
  ggplot(aes(x = sizefactor, y = colsum, color= Celltype)) +
  geom_point() +
  scale_y_log10() +
  theme_bw() +
  xlab("Size Factors") +
  ylab("Library Size")


#normalized = TRUE let normlizing on size factors
logcounts = log2(counts(dds, normalized = TRUE) + 1)

plot(logcounts[,1], logcounts[,2], cex = .1)

#stablizing variannce
rld = rlog(dds)

#pca
pc = prcomp(t(assay(rld)))
pc$x %>%
  as.data.frame() %>%
  mutate(Celltype = as.data.frame(ann)$Celltype,
         id = as.data.frame(ann)$Donorid) %>%
  ggplot(aes(x = PC1, y = PC2, color = Celltype)) +
  geom_point()

####################################################
##    Differential expression  using DEseq2
####################################################

#remove whole blood
wb = ann$Title[ann$Celltype == "Whole Blood"]

dta6 = dta5[,which(!colnames(dta5) %in% wb)]
dta6 = apply(dta6, 2, as.numeric)
rownames(dta6) = rownames(dta5)

ann1 = ann[!ann$Title %in% wb,]
ann1 = ann1 %>%
  mutate(neut = ifelse(Celltype == "Neutrophils",1,0),
         mono = ifelse(Celltype == "Monocytes",1,0),
         bcells = ifelse(Celltype == "B-cells",1,0),
         cd4 = ifelse(Celltype == "CD4",1,0),
         cd8 = ifelse(Celltype == "CD8",1,0),
         nk = ifelse(Celltype == "NK",1,0))

se1 = SummarizedExperiment(dta6)
colData(se1) = DataFrame(ann1)

#neut
dds1 = DESeqDataSet(se1, design = ~ neut)
dds1 = DESeq(dds1)
res = results(dds1)
resBigFC <- results(dds1, lfcThreshold=1, altHypothesis="greaterAbs")
#plotMA(resBigFC, ylim = c(-5,5))
resSort <- res[order(res$pvalue),]
#head(resSort)
#Count for the first gene: the unnormalized count
k <- counts(dds1)[rownames(resSort)[1],]
#Make a stripchart by combining time and protocol
cond <- with(colData(se1), factor(neut))
#stripchart(log2(k + 1) ~ cond, method="jitter", vertical=TRUE, las=2)

#The first 20 genes according to the lowest p-value
library(org.Hs.eg.db)
require(pheatmap)
require(BiomaRt)
#mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

geneinfo <- biomaRt::select(org.Hs.eg.db, 
                            keys=rownames(resSort)[1:20],
                            columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                            keytype="ENSEMBL")

geneinfo  = geneinfo %>%
  mutate(Celltype = "neut")

#mono
dds1 = DESeqDataSet(se1, design = ~ mono)
dds1 = DESeq(dds1)
res = results(dds1)
resBigFC <- results(dds1, lfcThreshold=1, altHypothesis="greaterAbs")
resSort <- res[order(res$pvalue),]
k <- counts(dds1)[rownames(resSort)[1],]
cond <- with(colData(se1), factor(neut))
geneinfo <- geneinfo %>%
  bind_rows(biomaRt::select(org.Hs.eg.db, keys=rownames(resSort)[1:20],
                            columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                            keytype="ENSEMBL") %>%
              mutate(Celltype = "mono"))

#bcells
dds1 = DESeqDataSet(se1, design = ~ bcells)
dds1 = DESeq(dds1)
res = results(dds1)
resBigFC <- results(dds1, lfcThreshold=1, altHypothesis="greaterAbs")
resSort <- res[order(res$pvalue),]
k <- counts(dds1)[rownames(resSort)[1],]
cond <- with(colData(se1), factor(neut))
geneinfo <- geneinfo %>%
  bind_rows(biomaRt::select(org.Hs.eg.db, keys=rownames(resSort)[1:20],
                            columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                            keytype="ENSEMBL") %>%
              mutate(Celltype = "bcells"))
#cd4
dds1 = DESeqDataSet(se1, design = ~ cd4)
dds1 = DESeq(dds1)
res = results(dds1)
resBigFC <- results(dds1, lfcThreshold=1, altHypothesis="greaterAbs")
resSort <- res[order(res$pvalue),]
k <- counts(dds1)[rownames(resSort)[1],]
cond <- with(colData(se1), factor(neut))
geneinfo <- geneinfo %>%
  bind_rows(biomaRt::select(org.Hs.eg.db, keys=rownames(resSort)[1:20],
                            columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                            keytype="ENSEMBL") %>%
              mutate(Celltype = "cd4"))
#cd8
dds1 = DESeqDataSet(se1, design = ~ bcells)
dds1 = DESeq(dds1)
res = results(dds1)
resBigFC <- results(dds1, lfcThreshold=1, altHypothesis="greaterAbs")
resSort <- res[order(res$pvalue),]
k <- counts(dds1)[rownames(resSort)[1],]
cond <- with(colData(se1), factor(neut))
geneinfo <- geneinfo %>%
  bind_rows(biomaRt::select(org.Hs.eg.db, keys=rownames(resSort)[1:20],
                            columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                            keytype="ENSEMBL") %>%
              mutate(Celltype = "cd8"))

#nk
dds1 = DESeqDataSet(se1, design = ~ bcells)
dds1 = DESeq(dds1)
res = results(dds1)
resBigFC <- results(dds1, lfcThreshold=1, altHypothesis="greaterAbs")
resSort <- res[order(res$pvalue),]
k <- counts(dds1)[rownames(resSort)[1],]
cond <- with(colData(se1), factor(neut))
geneinfo <- geneinfo %>%
  bind_rows(biomaRt::select(org.Hs.eg.db, keys=rownames(resSort)[1:20],
                            columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                            keytype="ENSEMBL") %>%
              mutate(Celltype = "nk"))

write.csv(geneinfo, file.path(dtafolder,"results","DEGeneList.csv"), row.names = FALSE)
geneinfo = read.csv(file.path(dtafolder,"results","DEGeneList.csv"))


#check the genelist
tmp = as.data.frame(apply(dta6, 2, function(x) log2(x + 0.5)))
rownames(tmp) = rownames(dta6)
tmp = tmp %>%
  filter(rownames(.) %in% geneinfo$ENSEMBL) %>%
  as.data.frame()

rownames(ann1) = ann1$Title
pheatmap(tmp, 
         annotation_col = as.data.frame(ann1[,c("Donorid","Celltype")]),
         show_colnames = FALSE,
         show_rownames = FALSE)

############################################################
##    Normalization using TMM
############################################################
#dta5 is prefiltered count data
require(edgeR)
dge = DGEList(counts = dta5)
dge = calcNormFactors(dge, method = "TMM")
dta_tmm = cpm(dge, log = TRUE)

dge$samples = dge$samples %>%
  rownames_to_column(var = "Title") %>%
  left_join(ann %>%
              select(Title, Celltype)) 

#limma-voom
Celltype = dge$samples$Celltype
mm <- model.matrix(~ 0 + Celltype)
y = voom(dge, mm, plot = TRUE)

fit = lmFit(y, mm)
head(coef(fit))
contr = makeContrasts("CelltypeB-cells" - "CelltypeCD4", levels = colnames(coef(fit)))
     
##########################################
##    TPM  
##########################################

dta4 = dta[filter,] 
rownames(dta4) = dta[filter,1]
dta4 = dta4[,-1]
ensembl_list = rownames(dta4)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_coords = getBM(attributes = c("hgnc_symbol",
                                   "ensembl_gene_id",
                                   "transcript_length",
                                   "start_position",
                                   "end_position"),
                    filter = "ensembl_gene_id",
                    mart = human,
                    values = ensembl_list)

gene_coords = gene_coords %>%
  mutate(size = end_position - start_position)

dta5 = dta4 %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  left_join(gene_coords %>%
              select(ensembl_gene_id, size) %>%
              distinct()) %>%
  distinct() 

#calculate TPM
dta5[,2:135] = apply(dta5[,2:135], 2, function(x) 1e6*x/sum(x))
rownames(dta5) = dta5$ensembl_gene_id
dta5 = dta5[,-1]

###########################################################################
## restore sample correlation - how is preprocessing and normalization?  
###########################################################################
require(pheatmap)
set.seed(1234)

hkgenes = read.csv(file.path(dtafolder,"HKgenes.csv"))
geneinfo = read.csv(file.path(dtafolder,"results","DEGeneList.csv"))

dta_ori = log2(dta[,-1] + 0.5)
rownames(dta_ori) = dta$genenames
dta_prefilter = log2(dta4 + 0.5)
dta_vst = log2(assay(rld))
dta_tpm = log2(dta5 + 0.5)
dta_tpm = dta_tpm[,1:134]

#function to get correlation by sampling
getCorbySampling = function(data, nsampling, nsize, method = "pearson"){
  n = 1
  dta = data[sample(1:nrow(data), nsize, replace = FALSE),]
  t1a = as.matrix(cor(dta, method = method, use = "complete.obs")) 
  t1a = reshape2::melt(t1a)  
  colnames(t1a)[3] = paste0("sampling_",n)
  n = 2
  while (n < nsampling){
    dta = data[sample(1:nrow(data), nsize, replace = FALSE),]
    t1 = as.matrix(cor(dta, method = method, use = "complete.obs")) 
    t1 = reshape2::melt(t1) 
    colnames(t1)[3] = paste0("sampling_",n)
    t1a  = t1a %>%
      left_join(t1)
    n = n + 1
  }
  
  t1a = t1a %>%
    filter(Var1 != Var2)  
  t1a = t1a %>%
    mutate(dist = rowMeans(t1a[,3:ncol(t1a)], na.rm = TRUE) ) %>%
    select(Var1, Var2, dist)
  
  #remove duplicates
  t3 = data.frame(Var1 = t(combn(unique(t1a$Var1), 2))[,1],
                  Var2 = t(combn(unique(t1a$Var1), 2))[,2]) %>%
    left_join(t1a)
  
  return(t3)
}
getCorForGenePanel = function(data, method = "pearson", log = FALSE,genepanel){
  data = as.data.frame(data[rownames(data) %in% genepanel,])
  cormat = as.data.frame(t(combn(colnames(data),2)))
  
  # if (log){
  #   corm = reshape2::melt(cor(log2(data + 0.5)))
  # }else{
  #   corm = reshape2::melt(cor(data))
  # }
  # cormat = cormat %>%
  #   left_join(as.data.frame(corm), by = c("V1" = "Var1", "V2" ="Var2"))
  
  cormat$dist = NA
  for (i in 1:nrow(cormat)){
    if (log){
      cormat$dist[i] = cor(log2(data[,cormat$V1[i]] + 1), log2(data[,cormat$V2[i]] + 1),use = "complete.obs", method = method)
    }else{
      #print(data[,cormat$V1[i]])
      cormat$dist[i] = cor(x = data[,cormat$V1[i]], y = data[,cormat$V2[i]], method = "pearson", use = "complete.obs")
    }
  }
  cormat = cormat %>%
    rename(Var1 = V1, 
           Var2 = V2)
  return(cormat)
}

# add back celltype from ann
addAnnotation = function(d){
  d = d %>%
    left_join(ann %>%
                select(Title, Celltype, Donorid), by = c("Var1" = "Title")) %>%
    rename(Celltype1 = Celltype,
           Donorid1  = Donorid) %>%
    left_join(ann %>%
                select(Title, Celltype, Donorid), by = c("Var2" = "Title")) %>%
    rename(Celltype2 = Celltype,
           Donorid2 = Donorid) %>%
    mutate(sametype = ifelse(Celltype1 == Celltype2, "YES","NO"),
           sameid = ifelse(Donorid1 == Donorid2, "YES","NO"))
  
  d = d %>%
    filter(Celltype2 != "Whole Blood") %>%
    group_by(Celltype1, sametype, sameid) %>%
    summarise(par = mean(dist, na.rm = TRUE),
              sd =  sd(dist, na.rm = TRUE))
  return(d)
}

random = lapply(list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm), getCorbySampling, nsampling = 1000, nsize = 100, method = "pearson")
names(random) = c("ori","prefilter","VST","TMM","TPM")

hk = vector(mode = "list", length = 5)
names(hk) = c("ori","prefilter","VST","TMM","TPM")
genepanel = hkgenes$ensembl_gene_id
for (i in 1:5){
  dd = list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm)[[i]]
  dd = as.data.frame(dd[rownames(dd) %in% genepanel,])
  cormat = as.data.frame(t(combn(colnames(dd),2)))
  corm = reshape2::melt(cor(dd))
  cormat = cormat %>%
    left_join(as.data.frame(corm), by = c("V1" = "Var1", "V2" ="Var2")) %>%
    rename(dist = value,
           Var1 = V1,
           Var2 = V2)
  hk[[i]] <- cormat
}

mono = lapply(list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm), getCorForGenePanel, method = "pearson", log = FALSE, genepanel = geneinfo$ENSEMBL[geneinfo$Celltype == "mono"])
names(mono) = c("ori","prefilter","VST","TMM","TPM")

d2_random = do.call(rbind, lapply(lapply(random, addAnnotation), function(x){
  x = x %>%
    filter(Celltype1 == "Monocytes")
  return(x)
}))
d2_random$ds = rep(c("ori","prefilter","vst","tpm","tmm"), each = 3)

d2_hk = do.call(rbind, lapply(lapply(hk, addAnnotation), function(x){
  x = x %>%
    filter(Celltype1 == "Monocytes")
  return(x)
}))
d2_hk$ds = rep(c("ori","prefilter","vst","tpm","tmm"), each = 3)

d2_mono = do.call(rbind, lapply(lapply(mono, addAnnotation), function(x){
  x = x %>%
    filter(Celltype1 == "Monocytes")
  return(x)
}))
d2_mono$ds = rep(c("ori","prefilter","vst","tpm","tmm"), each = 3)

d2 = d2_random %>%
  mutate(genelist = "random") %>%
  bind_rows(d2_hk %>%
              mutate(genelist = "hk")) %>%
  bind_rows(d2_mono %>%
              mutate(genelist = "mono")) %>%
  mutate(Celltype = case_when(sametype == "NO" & sameid == "NO" ~ "Diff ID, Diff Celltype",
                              sametype == "YES" & sameid == "NO" ~ "Diff ID, Same Celltype",
                              sametype == "NO" & sameid == "YES" ~ "Same ID, Diff Celltype"))


cor_plot = d2 %>%
  mutate(ds = factor(ds, levels = c("ori","prefilter","vst","tpm","tmm"))) %>%
  ggplot(aes(x = ds, y = par, group = genelist, color = genelist)) +
  geom_line() +
  geom_errorbar(aes(ymin = par - sd, ymax = par + sd), 
                width = 0.1) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Data Sources", y = "Correlation Coefficients") +
  facet_grid(~ Celltype)

temp_file = tempfile(fileext = ".png")
ggsave(temp_file, cor_plot, width = 8, height = 6, units = "in")
file.copy(temp_file, file.path(dtafolder,"results","Lineplot_correationBySteps.png"), overwrite = TRUE)

#########################################################################
###   Application 1:  Assessing DE genelists 
#########################################################################

#function to extract top 20 DE genes with different dataset
ann2 = ann %>%
  filter(Celltype != "Whole Blood") %>%
  mutate(mono = ifelse(Celltype == "Monocytes", 1, 0)) %>%
  select(Title, Donorid, mono)

getDEGenes = function(data, log = FALSE, annotation){
  require(nlme)
  require(lmerTest)
  data = data[,annotation$Title]
  pvalues = apply(data,1, function(x) {
    d = data.frame(gx = as.numeric(x),
                   group = annotation$mono,
                   id = annotation$Donorid) %>%
      mutate(id = factor(id))
    
    #mixed effect with nlme
    myfit = try(lme(gx ~ group, random  = ~1|id,data = d))
    pv = ifelse(class(myfit) != "try-error", anova(myfit)$'p-value'[2], NA)
    direction = ifelse(class(myfit) != "try-error", ifelse(fixef(myfit)[2] < 0, "Negative","Positive"), NA)
    
    return(list(pvalues = pv, direction = direction))
  })
  
  dd = as.data.frame(do.call(rbind, pvalues)) %>%
    rownames_to_column(var = "gene") %>%
    arrange(pvalues) 
  return(dd)
}

de_prefilter = getDEGenes(data = dta_prefilter,ann = ann2)
de_vst = getDEGenes(data = dta_vst, annotation = ann2)
de_tpm = getDEGenes(data = dta_tpm, annotation = ann2)
de_tmm = getDEGenes(data = dta_tmm, annotation = ann2)

DEList = lapply(list(de_prefilter, de_vst, de_tpm, de_tmm), function(x) {
  x$gene[1:20]
})
names(DEList) = c("prefilter","vst","tpm", "tmm")

cor_de_prefilter = lapply(list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm), getCorForGenePanel, method = "pearson", log = FALSE, genepanel = DEList[[1]])
names(cor_de_prefilter) = c("ori","prefilter","VST","TMM","TPM")

cor_de_vst = lapply(list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm), getCorForGenePanel, method = "pearson", log = FALSE, genepanel = DEList[[2]])
names(cor_de_vst) = c("ori","prefilter","VST","TMM","TPM")

cor_de_tpm = lapply(list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm), getCorForGenePanel, method = "pearson", log = FALSE, genepanel = DEList[[3]])
names(cor_de_tpm) = c("ori","prefilter","VST","TMM","TPM")

cor_de_tmm = lapply(list(dta_ori, dta_prefilter, dta_vst,dta_tpm,dta_tmm), getCorForGenePanel, method = "pearson", log = FALSE, genepanel = DEList[[4]])
names(cor_de_tmm) = c("ori","prefilter","VST","TMM","TPM")

summarizeCorrelation = function(cordata){
  cordata = do.call(rbind, lapply(lapply(cordata, addAnnotation), function(x){
    x = x %>%
      filter(Celltype1 == "Monocytes")
    return(x)
  }))
  cordata$ds = rep(c("ori","prefilter","vst","tpm","tmm"), each = 3)
  return(cordata)
  
}

cor_de = do.call(rbind,lapply(list(cor_de_prefilter, cor_de_vst, cor_de_tpm, cor_de_tmm), summarizeCorrelation))
cor_de$genelist = rep(c("prefilter","vst","tpm","tmm"), each = 15)
cor_de = cor_de %>%
  mutate(Celltype = case_when(sametype == "NO" & sameid == "NO" ~ "Diff ID, Diff Celltype",
                              sametype == "YES" & sameid == "NO" ~ "Diff ID, Same Celltype",
                              sametype == "NO" & sameid == "YES" ~ "Same ID, Diff Celltype")) %>%
  mutate(ds = factor(ds, levels = c("ori","prefilter","vst","tpm", "tmm")))


pd = position_dodge(0.2)
cor_plot_eval = cor_de %>%
  ggplot(aes(x = ds, y = par, group = genelist, color = genelist)) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = par - sd, ymax = par + sd), width = 0.2, position = pd) +
  theme_bw() +
  facet_grid(  ~ Celltype)

temp_file = tempfile(fileext = ".png")
ggsave(temp_file, cor_plot_eval, width = 8, height = 6, units = "in")
file.copy(temp_file, file.path(dtafolder,"results","Lineplot_correation_evaluation.png"), overwrite = TRUE)

#venn diagram
#devtools::install_github("yanlinlin82/ggvenn")
require(ggvenn)
venn_degenes = ggvenn(DEList)
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, venn_degenes, width = 8, height = 6, units = "in")
file.copy(temp_file, file.path(dtafolder,"results","VennDiagram_correation_evaluation.png"), overwrite = TRUE)


#######################################################################
##    Application 2: optimizing the DE gene panel
#######################################################################
degenes_tmm = de_tmm$gene[1:400]

#getCorForGenePanel = function(data, method = "pearson", log = FALSE,genepanel)

set.seed(1234)
sampledgenes = cor_sampledgenes = vector(mode = "list", length = 1000)
for (i in 1:1000){
  sampledgenes[[i]] = sample(degenes_tmm, 20, replace = FALSE)
  d = getCorForGenePanel(data = dta_tmm,genepanel = sampledgenes[[i]])
  d1 = addAnnotation(d) %>%
    filter(Celltype1 == "Monocytes") %>%
    mutate(sn = i)
  cor_sampledgenes[[i]] = d1
}
  
cor_selected = do.call(rbind, cor_sampledgenes)

g_cor_selected = cor_selected %>%
  mutate(Celltype = case_when(sametype == "NO" & sameid == "NO" ~ "Diff ID, Diff Celltype",
                              sametype == "YES" & sameid == "NO" ~ "Diff ID, Same Celltype",
                              sametype == "NO" & sameid == "YES" ~ "Same ID, Diff Celltype")) %>%
  ggplot(aes(x = par)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Correlation Coefficients") +
  ggh4x::facet_grid2(~ Celltype, scales = "free_y", independent = "y")

temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g_cor_selected, width = 8, height = 6, units = "in")
file.copy(temp_file, file.path(dtafolder,"results","Histogram_correation_optimizingGeneList.png"), overwrite = TRUE)

#get the best panel
m = 2
selected_sn1 = cor_selected %>%
  filter(sametype == "NO", sameid == "NO") %>%
  filter(par > 0, par < 0.1) 
gp1 = sampledgenes[[selected_sn1$sn[m]]]

selected_sn2 = cor_selected %>%
  filter(sametype == "NO", sameid == "NO") %>%
  filter(par > 0.64, par < 0.65) 
gp2 = sampledgenes[[selected_sn2$sn[m]]]

venn_selectedgenes = ggvenn(list(better = gp1, worse = gp2, top20 = de_tmm$gene[1:20]))
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, venn_selectedgenes, width = 8, height = 6, units = "in")
file.copy(temp_file, file.path(dtafolder,"results","Venndiagram_correation_selectedPanel.png"), overwrite = TRUE)


ddd1 = dta_tmm[rownames(dta_tmm) %in% gp1,colnames(dta_tmm) %in% ann2$Title]
ddd2 = dta_tmm[rownames(dta_tmm) %in% gp2,colnames(dta_tmm) %in% ann2$Title]
ddd3 = dta_tmm[rownames(dta_tmm) %in% de_tmm$gene[1:20],colnames(dta_tmm) %in% ann2$Title]
coef_ddd1 = c()
for (i in 1:nrow(ddd1)){
  coef_ddd1 = c(coef_ddd1,lm(ddd1[i,] ~ ann2$mono)$coefficients[2])
}

coef_ddd2 = c()
for (i in 1:nrow(ddd2)){
  coef_ddd2 = c(coef_ddd2,lm(ddd2[i,] ~ ann2$mono)$coefficients[2])
}

coef_ddd3 = c()
for (i in 1:nrow(ddd3)){
  coef_ddd3 = c(coef_ddd3,lm(ddd3[i,] ~ ann2$mono)$coefficients[2])
}


ddd4 = data.frame(cor = c(coef_ddd1, coef_ddd2, coef_ddd3), ann = rep(c("better", "worse","top20"), each = 20)) %>%
  mutate(direction = ifelse(cor >=0,"Positive","Negative"))

boxplot_LFC_selectedPanels = ddd4 %>%
  ggplot(aes(x = ann, y = cor, color = direction)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(0.2)) +
  theme_bw() +
  labs(x = "Gene Panels", y = "Log-fold Change")
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, boxplot_LFC_selectedPanels, width = 8, height = 6, units = "in")
file.copy(temp_file, file.path(dtafolder,"results","Boxplot_LFC_selectedPanel.png"), overwrite = TRUE)



require(cowplot)
gg1 = plot_grid(venn_selectedgenes, boxplot_LFC_selectedPanels)
temp_file = tempfile(fileext = ".png")
save_plot(temp_file, gg1, base_width = 8, base_height = 6)
file.copy(temp_file, file.path(dtafolder,"results","Composite_selectedPanels_optimizingGeneList.png"), overwrite = TRUE)



###################################################################
##  Considering positive and negative correlation
###################################################################
combineMeanSD = function(m1, m2, n1, n2, sd1, sd2){
  mean = (n1*m1 + n2*m2)/(n1 + n2)
  sd = sqrt(((n1 - 1)*sd1*sd1 + (n2 - 1)*sd2*sd2 + n1*n2/(n1 + n2)*(m1^2 + m2^2 - 2*m1*m2))/(n1 + n2 - 1))
  return(data.frame(mean = mean, sd = sd))
}

getCorForGenePanelByDirection = function(data, method = "pearson", log = FALSE,genepanel){
  genepanel = as.data.frame(genepanel)
  positive_panel = genepanel$gene[genepanel$direction == "Positive"]
  negative_panel = genepanel$gene[genepanel$direction == "Negative"]
  
  
  cormat = lapply(list(positive_panel, negative_panel),function(x) getCorForGenePanel(data = data, log = log, genepanel = x))
  
  cormat = addAnnotation(cormat[[1]]) %>%
    rename(par_positive = par,
           sd_positive = sd) %>%
    left_join(addAnnotation(cormat[[2]]) %>%
                rename(par_negative = par,
                       sd_negative = sd) )  
  cormat = cormat %>%
    mutate(n_positive = length(positive_panel),
           n_negative = length(negative_panel))
  
  cormat = cbind(cormat, combineMeanSD(m1 = cormat$par_positive,
                              n1 = cormat$n_positive,
                              m2 = cormat$par_negative,
                              n2 = cormat$n_negative,
                              sd1 = cormat$sd_positive,
                              sd2 = cormat$sd_negative))
  return(cormat)
}

cor_prefilter2 = getCorForGenePanelByDirection(data = dta_prefilter,
                                         genepanel = de_tmm[1:20,]) %>%
  filter(Celltype1 == "Monocytes") 

cor_vst2 = getCorForGenePanelByDirection(data = dta_vst,
                                         genepanel = de_tmm[1:20,]) %>%
  filter(Celltype1 == "Monocytes") 


cor_tmm2 = getCorForGenePanelByDirection(data = dta_tmm,
                                        genepanel = de_tmm[1:20,]) %>%
  filter(Celltype1 == "Monocytes") 

cor_tpm2 = getCorForGenePanelByDirection(data = dta_tpm,
                                         genepanel = de_tmm[1:20,]) %>%
  filter(Celltype1 == "Monocytes") 

cor2 = cor_prefilter2 %>%
  mutate(ds = "prefilter") %>%
  bind_rows(cor_vst2 %>%
              mutate(ds = "vst")) %>%
  bind_rows(cor_tpm2 %>%
              mutate(ds = "tpm")) %>%
  bind_rows(cor_tmm2 %>%
              mutate(ds = "tmm"))

cor2 %>%
  pivot_longer(cols = c("par_positive","par_negative","mean"), names_to = "measurement", values_to = "pvalue") %>%
  filter(sametype == "NO", sameid == "NO") %>%
  mutate(ds = factor(ds, levels = c("prefilter","vst","tpm","tmm"))) %>%
  ggplot(aes(x = ds, y =pvalue, group = measurement, color = measurement)) +
  geom_line()
  

######################################################
##    save project
######################################################
save.image(file = file.path(dtafolder,"Monitoring RNAseq Normalization","Monitoring RNAseq Normalization.RData"))
load(file = file.path(dtafolder,"Monitoring RNAseq Normalization","Monitoring RNAseq Normalization.RData"), verbose = TRUE)






t3 = cor_hk_ori %>%
  dplyr::rename(hk_ori = corr) %>%  
  full_join(cor_hk_prefilter %>%
              dplyr::rename(hk_prefilter = corr)) %>%
  full_join(cor_hk_normalized %>%
              dplyr::rename(hk_normalized = corr)) %>%
  full_join(cor_hk_tpm %>%
              dplyr::rename(hk_tpm = corr)) %>%
  full_join(cor_mono_ori %>%
              dplyr::rename(mono_ori = corr)) %>%
  full_join(cor_mono_prefilter %>%
              dplyr::rename(mono_prefilter = corr)) %>%
  full_join(cor_mono_normalized %>%
              dplyr::rename(mono_normalized = corr)) %>%
  full_join(cor_mono_tpm %>%
              dplyr::rename(mono_tpm = corr)) %>%
  filter(is.finite(mono_ori)) %>%
  full_join(cor_random_ori %>%
              dplyr::rename(V1 = Var1,
                            V2 = Var2,
                            random_ori = dist)) %>%
  full_join(cor_random_prefilter %>%
              dplyr::rename(V1 = Var1,
                            V2 = Var2,
                            random_prefilter = dist)) %>%
  full_join(cor_random_normalization %>%
              dplyr::rename(V1 = Var1,
                            V2 = Var2,
                            random_normalized = dist)) %>%
  full_join(cor_random_tpm %>%
              dplyr::rename(V1 = Var1,
                            V2 = Var2,
                            random_tpm = dist)) 


t3 = t3 %>%
  left_join(alldta$ann %>%
              select(Title, Donorid, Celltype), by = c("V1" = "Title")) %>%
  rename(Donorid1 = Donorid,
         Celltype1 = Celltype) %>%
  left_join(alldta$ann %>%
              select(Title, Donorid, Celltype), by = c("V2" = "Title")) %>%
  rename(Donorid2 = Donorid,
         Celltype2 = Celltype)  %>%
  filter(Celltype1 != "Whole Blood") %>%
  filter(Celltype2 != "Whole Blood") %>%
  mutate(sameid = Donorid1 == Donorid2,
         samecelltype = Celltype1 == Celltype2) %>%
  mutate(comp = case_when(
    sameid & !samecelltype ~ "Same ID, Diff Celltype",
    !sameid & samecelltype ~ "Diff ID, Same Celltype",
    !sameid & !samecelltype ~ "Diff ID, Diff Celltype"
  ))

write.csv(t3, file.path(dtafolder,"results","CorForQC.csv"), row.names = FALSE)
t3 = read.csv(file.path(dtafolder,"results","CorForQC.csv"))

t4 = t3 %>%
  filter(!(Celltype1 == "Monocytes" & Celltype2 == "Monocytes")) %>%
  group_by(comp) %>%
  summarise(hk_ori = mean(hk_ori,na.rm = TRUE),
            hk_prefilter = mean(hk_prefilter, na.rm = TRUE),
            hk_normalized = mean(hk_normalized, na.rm = TRUE),
            hk_tpm = mean(hk_tpm, na.rm = TRUE),
            mono_ori = mean(mono_ori,na.rm = TRUE),
            mono_prefilter = mean(mono_prefilter, na.rm = TRUE),
            mono_normalized = mean(mono_normalized, na.rm = TRUE),
            mono_tpm = mean(mono_tpm, na.rm = TRUE),
            random_ori = mean(random_ori, na.rm = TRUE),
            random_prefilter = mean(random_prefilter,na.rm = TRUE),
            random_normalized = mean(random_normalized,na.rm = TRUE),
            random_tpm = mean(random_tpm,na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(comp))

t5 = t4 %>%
  pivot_longer(cols = 2:13,names_to = "name") %>%
  mutate(step = ifelse(grepl("ori", name), "ori", ifelse(grepl("normalized", name), "normalized", ifelse(grepl("tpm", name),"tpm", "prefilter"))),
         panel = ifelse(grepl("hk", name), "hk", ifelse(grepl("mono", name), "mono", ifelse(grepl("random", name),"random", NA_character_)))) %>%
  mutate(step = factor(step, levels = c("ori","prefilter","normalized","tpm")))

p1 = t5 %>%
  ggplot(aes(x = step, y = value, group = panel, color = panel)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggh4x::facet_grid2(~ comp)




t4 = t3 %>%
  left_join(ann %>%
              select(Title, Donorid, Celltype), by = c("V1" = "Title")) %>%
  rename(Donorid1 = Donorid,
         Celltype1 = Celltype) %>%
  left_join(ann %>%
              select(Title, Donorid, Celltype), by = c("V2" = "Title")) %>%
  rename(Donorid2 = Donorid,
         Celltype2 = Celltype)  %>%
  mutate(sameid = Donorid1 == Donorid2,
         samecelltype = Celltype1 == Celltype2) %>%
  mutate(comp = case_when(
    sameid & !samecelltype ~ "Same ID, Diff Celltype",
    !sameid & samecelltype ~ "Diff ID, Same Celltype",
    !sameid & !samecelltype ~ "Diff ID, Diff Celltype"
  ))

t5 = t4 %>%
  #filter(ds != "genelist") %>%
  group_by(comp) %>%
  summarise(ori = mean(ori, na.rm = TRUE),
            prefilter = mean(prefilter, na.rm = TRUE),
            normalized = mean(normalized,na.rm = TRUE),
            normalized_hk = mean(normalized_hk,na.rm = TRUE),
            ori_hk = mean(ori_hk,na.rm = TRUE),
            normalized_genelist = mean(normalized_genelist,na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(comp))

t5 = t5 %>%
  pivot_longer(cols = 2:7, names_to = "steps", values_to = "correlation")


g1 = t5 %>%
  ggplot(aes(x = steps, y = correlation, group = comp, color = comp)) +
  geom_point() +
  geom_line() +
  theme_bw()



mutate(ds = factor(ds, levels = c("ori","prefilter","normalization","genelist"))) %>%
  ggplot(aes(x = ds, y = cor, group = comp, color = comp)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") %>%
  labs(x = "Data Processing", y = "Pearson Correlation", color = "Comparisons") +
  ylim(c(0,1))

# t4_full = t4 %>%
#   group_by(comp) %>%
#   summarise(meancor = mean(dist),
#             sdcor = sd(dist))

t6 = t4 %>%
  filter(!sameid) %>%
  #mutate(ct = paste0(Celltype1, "-",Celltype2)) %>%
  group_by(ds, Celltype1, Celltype2) %>%
  summarise(cor = mean(dist, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ds = factor(ds, levels = c("ori","prefilter","normalization","genelist")))

t6 %>%
  #pivot_wider(names_from = Celltype2, values_from = cor) %>%
  ggplot(aes(x = Celltype1, y = Celltype2, fill = cor)) +
  geom_tile() +
  facet_wrap(~ ds)


g2 = t6 %>%
  ggplot(aes(x = ct, y = dist)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2), alpha = 0.2) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)) 
g2

combn(unique(t6$Celltype1))














```

##    Constrained Least square regression

```{r echo = FALSE, message = FALSE}
#install.packages("quadprog")
library("quadprog")

d = lapply(unique(alldta$ann$Celltype), function(m){
  d = dta5[,colnames(dta5) %in% alldta$ann$Title[alldta$ann$Celltype == m]]
  d = log2(d + 1)  
  d = apply(d,1, median)
  d1 = data.frame(feature = dta5$ensembl_gene_id,
                  value = d) 
  colnames(d1)[2] = m
  return(d1)
})
names(d) = unique(alldta$ann$Celltype)
X = d[[2]] %>%
  left_join(d[[3]]) %>%
  left_join(d[[4]]) %>%
  left_join(d[[5]]) %>%
  left_join(d[[6]]) %>%
  left_join(d[[7]])  %>%
  select(-feature) %>%
  as.matrix()
rownames(X) = rownames(dta5)

Y = dta5[,colnames(dta5) %in% alldta$ann$Title[alldta$ann$Celltype == "Whole Blood"]] 
Y = log2(Y + 1)
rownames(Y) = rownames(dta5)

Rinv <- solve(chol(t(X) %*% X));

C <- cbind(rep(1,ncol(X)), diag(ncol(X)))
b <- c(1,rep(0,ncol(X)))

i = 1
d <- t(Y[,i]) %*% X  
fit = solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
frac = data.frame(Celltype = colnames(d), fraction = round(fit$solution,2))
colnames(frac)[2] = colnames(Y)[i]

for (i in 2:ncol(Y)){
  d <- t(Y[,i]) %*% X  
  fit = solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  d1 = data.frame(Celltype = colnames(d), fraction = round(fit$solution,2))
  colnames(d1)[2] = colnames(Y)[i]
  frac = frac %>%
    left_join(d1)
}

frac = as.data.frame(t(frac))
colnames(frac) = frac[1,]
frac = frac[-1,]  




```
##    Deconvolution with granulator

```{r echo = FALSE,message = FALSE}
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("granulator")
library(granulator)
# load_ABIS()
# bulkRNAseq_ABIS[1:5, 1:5]
# sigMatrix_ABIS_S0[1:5, 1:5]
# groundTruth_ABIS[1:5, 1:5]
# 
# 
# sigList = list(
#   ABIS_S0 = sigMatrix_ABIS_S0,
#   ABIS_S1 = sigMatrix_ABIS_S1, 
#   ABIS_S2 = sigMatrix_ABIS_S2, 
#   ABIS_S3 = sigMatrix_ABIS_S3)
# 
# plot_similarity(sigMatrix=sigList)
# 
# decon <- deconvolute(m = bulkRNAseq_ABIS, sigMatrix = sigList)
# decon$proportions$svr_ABIS_S0[1:5, 1:5]
# plot_proportions(deconvoluted = decon, method = 'svr', signature = 'ABIS_S0')
# plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)

methods <- c('ols','nnls',"qprogwc",'qprog','rls','svr')
#decon <- deconvolute(bulkRNAseq_ABIS, list(ABIS_S2 = sigMatrix_ABIS_S2), methods)
#plot_deconvolute(decon, scale = TRUE, labels = FALSE)
decon = deconvolute(as.matrix(Y), X, methods)
plot_deconvolute(decon, scale = TRUE, labels = FALSE)

transformDta = function(data,method){
  data = data %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = colnames(data), names_to = "Celltype", values_to = "Fraction") %>%
    mutate(methods = method)
  return(data)
}

names(decon$proportions) = methods
frac2 = lapply(1:length(decon$proportions), function(x) transformDta(data = decon$proportions[[x]], method = names(decon$proportions)[x]))


frac2  = do.call(rbind, frac2)


```
##    Deconvolution with immunedeconv

```{r echo = FALSE, message = FALSE}

```



alldta = list(rld = rld,
              dta = dta,
              filter = filter,
              ann = ann)
save(alldta, file = file.path(dtafolder,"dataset","CompiledData.RData"))
load(file.path(dtafolder,"dataset","CompiledData.RData"), verbose = TRUE)

#random genes that are not in HK or DE
random = sample(rownames(alldta$dta_tpm), 90, replace = FALSE)
random =random[!random %in% c(hk$ensembl_gene_id, geneinfo$ENSEMBL)]
random = random[1:80]

de = unique(geneinfo$ENSEMBL)

genes = c(hk$ensembl_gene_id, de, random)
rowann = data.frame(ann = c(rep("HK", nrow(hk)),
                            rep("DE",length(de)),
                            rep("Random",length(random))))
rownames(rowann) = genes

dta5a = alldta$dta_tpm[rownames(alldta$dta_tpm) %in% genes,-135]

colann = data.frame(cn = colnames(dta5a)) %>%
  left_join(alldta$ann %>%
              select(Title, Celltype), by = c("cn" = "Title"))
rownames(colann) = colann$cn
colann = colann %>% select(Celltype)

pheatmap(log2(dta5a + 0.5),
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         annotation_row = rowann,
         annotation_col = colann,
         legend = FALSE)

```



#########################################################
##    Sample correlation in TCGA and GTEX  
#########################################################
require(TCGAbiolinks)

projects = getGDCprojects()
projects = projects %>%
  filter(grepl("TCGA",id))
tcgafolder = "/Volumes/home/Drive/learning_lmd/dataset/TCGA/RNAseq"

for (i in 1:nrow(projects)){
  print(projects$id[i])
  destfile = file.path(tcgafolder,paste0(projects$id[i],".RData"))
  if(file.exists(destfile)) next
  
  query <- GDCquery(project = projects$id[i], 
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "STAR - Counts")
  GDCdownload(query)
  data <- GDCprepare(query)
  save(data, file = destfile)
}

##GTEX

require(rvest)
require(openxlsx)
gtexfolder = "/Volumes/home/Drive/learning_lmd/dataset/GTEX"
ann = read.xlsx(file.path(gtexfolder,"tissuetypes.xlsx"))
tissues = gsub(" - ","_",ann$Tissue)
tissues = gsub("\\(|\\)","",tissues)
tissues = tolower(gsub(" ","_",tissues))

#counts
urls = paste0("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_",
              tissues,
              ".gct.gz")

lapply(1:length(urls), function(x) {
  print(tissues[x])
  destfile = paste0("gene_reads_",tissues[x],".gct.gz")
  destfile = file.path(gtexfolder,"counts",destfile)
  if (!file.exists(destfile))  download.file(urls[x], destfile = destfile)
})

#tpm
urls = paste0("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_",
              tissues,
              ".gct.gz")

lapply(1:length(urls), function(x) {
  print(tissues[x])
  destfile = paste0("gene_reads_",tissues[x],".gct.gz")
  destfile = file.path(gtexfolder,"tpm",destfile)
  if(!file.exists(destfile) ) download.file(urls[x], destfile = destfile)
})
```

##    Processing GTEX   

```{r echo = FALSE, message = FALSE}
gtexfolder = "/Volumes/home/Drive/learning_lmd/dataset/GTEX"
ctfiles = list.files(file.path(gtexfolder,"counts"))

#prefiltering function
prefilterData = function(dtafile){
  dtafile1 = file.path(gtexfolder,"counts",dtafile)
  dta = read.delim(gzfile(dtafile1), skip = 2)
  rownames(dta) = dta$Name
  dta = dta[,-c(1,2,3)]
  totbygene = rowSums(dta)
  quan30 = quantile(totbygene, probs = 0.3)
  
  filter = totbygene >= quan30
  dta = dta[filter,]
  
  filename = gsub("gene_reads","prefiltered", dtafile)
  filename = gsub("gct.gz","RData", filename)
  
  save(dta, file = file.path(gtexfolder,"prefilter",filename))
}

lapply(ctfiles, prefilterData)

#normalize using rld in DEseq2
pfiles = list.files(file.path(gtexfolder,"prefilter"))

for (i in 1:length(pfiles)){
  load(file.path(gtexfolder,"prefilter",pfiles[i]), verbose = TRUE)
  se =  SummarizedExperiment(as.matrix(dta))
  ann = data.frame(sample = colnames(dta),
                   tissue = gsub("prefiltered_|\\.RData","",pfiles[i]))
  rownames(ann) = ann$sampple
  colData(se) <- DataFrame(ann)
  dds = DESeqDataSet(se, design = ~ 1)
  rm(dta)
  rm(se)
  #dds = DESeq(dds)
  #stablizing variannce
  rld = vst(dds)
  fn = gsub("prefiltered","normalized",pfiles[i])
  save(rld, file = file.path("/Users/mingdongliu/Documents/normalized",fn))
}

```

##    DE with sampling  

```{r echo = FALSE, message = FALSE}


```

