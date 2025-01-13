########## DADA2 object building ############
library(dada2)
path <- "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/Cassiopea_Microbial_seqs_KMM_FLK_Aug2021" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
#Forward primer 17, reverse 21
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(285,230),
                     maxN=0, maxEE=c(4,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE,trimLeft = c(17,21)) # On Windows set multithread=FALSE
head(out)
saveRDS(out, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/out.rds")
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/errF.rds")
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errR, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/errR.rds")

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/dadaFs.rds")
saveRDS(dadaRs, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/dadaRs.rds")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
saveRDS(mergers, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/mergers.rds")
saveRDS(seqtab, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/seqtab.rds")
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 380:458]
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/seqtabnochim.rds")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#0.982436
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/kmmuf/Downloads/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
######## Tracking
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input","filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
saveRDS(track, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/track.rds")

########## Library load in
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(decontam)
library(dplyr)
library(cowplot)

### Build object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
sample_names(ps)
sampledata<-as.data.frame(read.csv("C:/Users/kmmuf/Downloads/CASS_METADATA.csv"))
rownames(sampledata)<-(sampledata$Sample_NAME)
samdat=sample_data(sampledata)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna, samdat)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/ps.rds")
psnomito<-subset_taxa(ps,Family!="Mitochondria"|is.na(Family))
psnochlor<-subset_taxa(psnomito,Class!="Chloroplast"|is.na(Class))
psnochlor2<-subset_taxa(psnochlor,Order!="Chloroplast"|is.na(Order))
psnochlor2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 34853 taxa and 96 samples ]:
#   sample_data() Sample Data:        [ 96 samples by 21 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 34853 taxa by 7 taxonomic ranks ]:
#   refseq()      DNAStringSet:       [ 34853 reference sequences ]
# taxa are columns
########## Contaminant removal #############
contamdf.prev <- isContaminant(psnochlor2, method="prevalence", neg="Neg")
table(contamdf.prev$contaminant)
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(psnochlor2, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Neg == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Neg == "FALSE", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
w<-ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, psnochlor2)
#47 ASVs trimmed out of 35443 ASVs in decontam
ps.noblanks<-subset_samples(ps.noncontam, Location!="NA")
ps.noblanks
psnomito.noblnk.ETOH<-subset_samples(ps.noblanks,Buffer=="ETOH")
psnomito.noblnk.DESS<-subset_samples(ps.noblanks,Buffer=="DESS")
psnomito.noblnk.DESS.Cx<-subset_samples(psnomito.noblnk.DESS,Species_HOST=="CX"|is.na(Species_HOST))
psnomito.noblnk.DESS.Cx
sum(psnomito.noblnk.DESS.Cx@otu_table)
# 6470034
saveRDS(psnomito.noblnk.DESS.Cx, "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/psnomito.noblnk.DESS.Cx.rds")
#34806 taxa
########### Rarefaction #################
library(vegan); packageVersion("vegan")
library("DESeq2")
veg_nb<-psotu2veg(psnomito.noblnk.DESS.Cx)
rarplot<-rarecurve(veg_nb, step=50, cex=0.5)
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/rarplot.png", plot = rarplot, width = 5, height = 5, units = "in", dpi = 300)

######## Alpha diversity #############
alphaPS.Cx<-psotu2veg(psnomito.noblnk.DESS.Cx)
alphaPS.Cx.s<-pssd2veg(psnomito.noblnk.DESS.Cx)
alphaPS.Cx.s$shannon<-vegan::diversity(x= alphaPS.Cx,index = "shannon")
kruskal.test(shannon~Body_cav,data = alphaPS.Cx.s)
dunn.test::dunn.test(alphaPS.Cx.s$shannon, alphaPS.Cx.s$Body_cav, method="holm")

# data: x and group
# Kruskal-Wallis chi-squared = 48.8388, df = 3, p-value = 0
# 
# 
# Comparison of x by group                            
# (Holm)                                     
# Col Mean-|
#   Row Mean |       Bell        GVC        Sub
# ---------+---------------------------------
#   GVC |   5.974667
# |    0.0000*
#   |
#   Sub |  -1.447190  -5.324066
# |     0.0739    0.0000*
#   |
#   WTR |   2.080691  -1.796184   2.807567
# |     0.0562     0.0725    0.0100*
#   
#   alpha = 0.05
# Reject Ho if p <= alpha/2

kruskal.test(shannon~Location,data = alphaPS.Cx.s)

xla<-plot_richness(psnomito.noblnk.DESS.Cx, x= "Location", measures=c("Observed","Shannon"), color="Body_cav")
yla<-xla + geom_boxplot()+ scale_color_manual(values=c("darkgreen","coral",'brown4','lightblue'))+theme_linedraw()+
  labs(color="Sample Type")+labs(title = "B.",y = " ")
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/Cass_shan_inv_sim_plot.png", plot = yla, width = 8, height = 5, units = "in", dpi = 300)
yla
sta<-plot_richness(psnomito.noblnk.DESS.Cx, x="Body_cav", measures=c("Observed","Shannon"), color="Body_cav")+
  geom_boxplot()+ scale_color_manual(values=c("darkgreen","coral",'brown4','lightblue'))+theme_linedraw()+
  labs(x="Sample Type",title = "A.")+guides(colour = "none")
sta

richal<-plot_grid(sta,yla,rel_widths = c(1,2))
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/Cass_shan_inv_sim_plot_F1.png", plot = richal, width = 9, height = 4, units = "in", dpi = 300)

rich = estimate_richness(psnomito.noblnk.DESS.Cx)
dunn.test::dunn.test(rich$Observed, sample_data(psnomito.noblnk.DESS.Cx)$Body_cav2, method = "BH")
# (Benjamini-Hochberg)                              
# Col Mean-|
#   Row Mean |       Bell        GVC        Sub
# ---------+---------------------------------
#   GVC |   4.913675
# |    0.0000*
#   |
#   Sub |  -1.461909  -4.650323
# |     0.0863    0.0000*
#   |
#   WTR |   2.713281  -0.475132   3.322710
# |    0.0050*     0.3173    0.0009*
#   
#   alpha = 0.05
# Reject Ho if p <= alpha/2

psnomito.noblnk.DESS.P = transform_sample_counts(psnomito.noblnk.DESS.Cx, function(x) 1 * x/sum(x))

set.seed(123)
BiocManager::install("DECIPHER")
library(DECIPHER)
#not using tree based approach


###Faiths D ####
library(btools)
library(MiscMetabar)
library(phangorn)
psnomito.noblnk.DESS.Cx.most<-prune_taxa(taxa_sums(psnomito.noblnk.DESS.Cx) > 25, psnomito.noblnk.DESS.Cx)
df_tree<-build_phytree_pq(psnomito.noblnk.DESS.Cx.most)
psnomito.noblnk.DESS.Cx.most2<-merge_phyloseq(psnomito.noblnk.DESS.Cx.most,df_tree$ML$tree)
saveRDS(psnomito.noblnk.DESS.Cx.most2,"C:/Users/kmmuf/Downloads/Cassiopea_tree.RDS")
library("ape")
library(metagMisc)
library(PhyloMeasures)
fd_cas<-phyloseq_phylo_div(
  psnomito.noblnk.DESS.Cx.most2,
  measures = c("PD"))

sample_data(psnomito.noblnk.DESS.Cx)$Faith<-fd_cas$PD
alphaPS.Cx.s$Faith<-fd_cas$PD
boxplot(alphaPS.Cx.s$Faith~alphaPS.Cx.s$Body_cav)

ggplot(sample_data(psnomito.noblnk.DESS.Cx), aes(Body_cav, Faith, fill = Body_cav))+geom_boxplot()
Faithplot<-ggplot(sample_data(psnomito.noblnk.DESS.Cx), aes(Body_cav, Faith))+geom_boxplot(aes(fill = Body_cav))+ scale_fill_manual(values=c("darkgreen","coral",'brown4','lightblue'))+
  theme_minimal()+xlab('Sample Type')+ylab("Faith's Diversity")
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/Cass_faith_plot.png", Faithplot, h=3, w=4, dpi=300, units = "in")
#############
alphaPS.Cx.s$Body_cav<-as.factor(alphaPS.Cx.s$Body_cav)
kruskal.test(alphaPS.Cx.s$Faith~alphaPS.Cx.s$Body_cav)
dunn.test::dunn.test(alphaPS.Cx.s$Faith, alphaPS.Cx.s$Body_cav, method = "Holm")

# Kruskal-Wallis rank sum test
# 
# data: x and group
# Kruskal-Wallis chi-squared = 44.1226, df = 3, p-value = 0
# 
# 
# Comparison of x by group                            
# (Holm)                                     
# Col Mean-|
#   Row Mean |       Bell        GVC        Sub
# ---------+---------------------------------
#   GVC |   5.822670
# |    0.0000*
#   |
#   Sub |  -0.484609  -4.262856
# |     0.6280    0.0001*
#   |
#   WTR |   3.470031  -0.308215   3.147192
# |    0.0010*     0.3790    0.0025*
#   
#   alpha = 0.05
# Reject Ho if p <= alpha/2

library(ggstatsplot)
vegcxtemp<-pssd2veg(psnomito.noblnk.DESS.Cx)
om<-ggbetweenstats(
  data = vegcxtemp,
  x = Body_cav,
  y = Faith,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
) + xlab("Sample Type")+ scale_colour_manual(name="Body_cav", values=c("darkgreen","coral",'brown4','lightblue'))+ ylab("Faith's Diversity Index")
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/Cass_faith_plot2.png",plot= om,height = 5,width = 6,units = "in",dpi=300)


# Calculate distance matrix ### Bray curtis dist ########################
brayDist <- phyloseq::distance(psnomito.noblnk.DESS.P, method = "bray")
# Calculate ordination
brayMDS  <- ordinate(psnomito.noblnk.DESS.P, "MDS", distance=brayDist)
## Make plot
# Create plot, store as temp variable
garf <- plot_ordination(psnomito.noblnk.DESS.P, brayMDS, color="Lat",shape="Body_cav")
garfcom<-garf + theme_linedraw()+labs(shape="Sample Type",color = "Latitude")+ggtitle("Bray-Curtis distance by sample type and latitude")+
  scale_colour_gradient(low = "coral", high = "navy")+geom_point(size=3)
garfcomseg<-garfcom+facet_wrap(~Body_cav)+ggtitle("Bray-Curtis distance split by sample type")+
  guides(color="none",shape="none")
garf_tog<-plot_grid(garfcom,garfcomseg,ncol=1)
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/BrayCur.png", plot=garf_tog, dpi=300, units="in", h=7, w=5)
######## Beta diversity ###########################
eucDist <- phyloseq::distance(psnomito.noblnk.DESS.P, method = "euclidean")
# Calculate ordination
eucMDS  <- ordinate(psnomito.noblnk.DESS.P, "PCoA", distance=eucDist)
## Make plot
# Create plot, store as temp variable, p
gare <- plot_ordination(psnomito.noblnk.DESS.P, eucMDS, color="Lat",shape="Body_cav")
garecom<-gare + theme_linedraw()+labs(shape="Sample Type",color = "Latitude")+ggtitle("Euclidean distance by sample type and latitude")+
  scale_colour_gradient(low = "coral", high = "navy")+geom_point(size=2)
garecomseg<-garecom+facet_wrap(~Body_cav)+ggtitle("Euclidean distance split by sample type")+
  guides(color="none",shape="none")
gare_tog<-plot_grid(garecom,garecomseg,ncol=1)
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/Euc_dist.png", plot=gare_tog, dpi=300, units="in", h=7, w=5)
####### Adonis2 ################################
library(tidyverse)
library(csv)
library(ape)
library("pairwiseAdonis")
library(ggpubr)

sam <- data.frame(sample_data(psnomito.noblnk.DESS.P))
adonis2(eucDist ~ Body_cav + Location + Body_cav*Location, data = sam, permutations = 999)

# adonis2(formula = eucDist ~ Body_cav + Location + Body_cav * Location, data = sam, permutations = 999)
#                   Df SumOfSqs      R2       F   Pr(>F)    
# Body_cav             3   2.3581 0.26248 11.5168  0.001 ***
#   Location           7   1.7455 0.19429  3.6535  0.001 ***
#   Body_cav:Location 21   1.8773 0.20897  1.3098  0.126    
# Residual            44   3.0030 0.33427                   
# Total               75   8.9839 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonisData <- otu_table(psnomito.noblnk.DESS.P) %>%                                      #use sequence table that has inocula removed
  as.data.frame() %>%
  rownames_to_column(var = "Sample_NAME") %>%                                          #change rownames to a column so there is a common variable to join by
  left_join(sam, by = "Sample_NAME") %>%                                            #join sample data to the sequence table
  select(-c("Sample_NAME", "Location", "Identity", "Buffer", "Species_HOST","Temp","Specimen","Lat","Long",
            "Salinity","pH","Temperature_.C.","Date","Species_O","Time","Diameter..cm.","OAL..cm.",
            "Depth..m.","Concentration..indv.m2.","Neg")) %>% #remove all metadata columns except the one to be used to compare
  select(Body_cav, everything())                                                   #rearrange so substrate column is first
dim(adonisData)
head(adonisData[,1:10])

pairwise.adonis2(adonisData[,2:34807]~Body_cav,method="euclidean",data=adonisData)
#Bell v GVC 0.001
#Bell vs Sub 0.061 
#Bell vs Water 0.012
#GVC vs Sub 0.002 
#GVC vs Water 0.001
#Sub v Water 0.003
##### Pairwise location comparisons 
psnomito.noblnk.DESS.P.GVC<-subset_samples(psnomito.noblnk.DESS.P,Body_cav=="GVC")
psnomito.noblnk.DESS.P.Bell<-subset_samples(psnomito.noblnk.DESS.P,Body_cav=="Bell")
#rerun with both GVC and bell
adonisData2 <- otu_table(psnomito.noblnk.DESS.P.GVC) %>%                                      #use sequence table that has inocula removed
  as.data.frame() %>%
  rownames_to_column(var = "Sample_NAME") %>%                                          #change rownames to a column so there is a common variable to join by
  left_join(sam, by = "Sample_NAME") %>%                                            #join sample data to the sequence table
  select(-c("Sample_NAME", "Body_cav", "Identity", "Buffer", "Species_HOST","Temp","Specimen","Lat","Long",
            "Salinity","pH","Temperature_.C.","Date","Species_O","Time","Diameter..cm.","OAL..cm.",
            "Depth..m.","Concentration..indv.m2.","Neg")) %>% #remove all metadata columns except the one to be used to compare
  select(Location, everything()) 
adonloc<-pairwise.adonis2(adonisData2[,2:34807]~Location,method="euclidean",data=adonisData2)
str(adonloc)
str(adonloc$BPK_vs_CK)
sandr<-as.data.frame(adonloc)
write.csv(sandr, file="C:/Users/kmmuf/Documents/Cass_mic_seq_test/locationadonis_GVC.csv")
#For internal only one is significant KL_vs_MK.Pr..F.	0.008
#No statistically significant differences for bells
#####Phylum barplots #################################
library(speedyseq)
ps.phyl <- psnomito.noblnk.DESS.Cx %>%
  tax_glom("Phylum") %>%
  transmute_tax_table(Kingdom, Phylum, .otu = Phylum)
ps.phyl
phylon<-plot_bar(ps.phyl,x="Location",fill="Phylum")+facet_wrap(~Body_cav)

ggsave("C:/Users/kmmuf/Downloads/Cass_shan_phylum_plot.png", plot = phylon, width = 12, height = 8, units = "in", dpi = 300)

phylon<-plot_bar(ps.phyl,x="Body_cav",fill="Phylum")+facet_grid(~Body_cav, scales = "free_x")
psphylRm = merge_samples(ps.phyl, "Body_cav")
sample_data(psphylRm)$Body_cav <- levels(sample_data(ps.phyl)$Body_cav)
psphylRm = transform_sample_counts(psphylRm, function(x) 100 * x/sum(x))

top10phy = names(sort(taxa_sums(ps.phyl), TRUE)[1:10])
b10=prune_taxa(top10phy, psphylRm)

title = "Phylum level community composition by sample type"
phyava<-plot_bar(b10, fill = "Phylum", title = title) + coord_flip() + 
  ylab("Percentage of Sequences") + scale_fill_brewer(palette ="PRGn")+theme_linedraw()
phyava
ggsave("~/phyava.png",plot = phyava, dpi = 300, h = 4, w= 5, units = 'in')
#Kingdom level
ps.kin <- psphylRm %>%
  tax_glom("Kingdom") %>%
  transmute_tax_table(Kingdom, .otu = Kingdom)
phyavaK<-plot_bar(ps.kin, fill = "Kingdom", title = "Kingdom level community composition by each sample type ") + coord_flip() + 
  ylab("Percentage of Sequences") +theme_linedraw() + scale_fill_manual(values = c("goldenrod","navy"))
phyavaK
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/phyavaK.png",plot = phyavaK, dpi = 300, h = 4, w= 5, units = 'in')
phykin<-plot_grid(phyava, phyavaK, ncol=1, rel_heights = c(1.5,1))
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/phykin.png",plot = phykin, dpi = 300, h = 6, w= 9, units = 'in')

############################################Class#########
ps.classes <- psnomito.noblnk.DESS.Cx %>%
  tax_glom("Class") %>%
  transmute_tax_table(Kingdom, Phylum, Class, .otu = Class)
psclsRm = merge_samples(ps.classes, "Body_cav")
sample_data(psclsRm)$Body_cav <- levels(sample_data(ps.classes)$Body_cav)
psclsRm = transform_sample_counts(psclsRm, function(x) 100 * x/sum(x))
top10cls = names(sort(taxa_sums(psclsRm), TRUE)[1:10])
c10=prune_taxa(top10cls, psclsRm)
summary(otu_table(c10))
##########Order and gut##############
guto<-subset_samples(psnomito.noblnk.DESS.Cx, Body_cav == "GVC")
ps.orders <- guto %>%
  tax_glom("Order") %>%
  transmute_tax_table(Kingdom, Phylum, Class, Order, .otu = Order)
ps.orders = transform_sample_counts(ps.orders, function(x) 100 * x/sum(x))
top10os = names(sort(taxa_sums(ps.orders), TRUE)[1:10])
ps.oo=prune_taxa(top10os, ps.orders)
summary(otu_table(ps.oo))
gutop<-subset_taxa(guto, Order == "Pseudomonadales")
gutop = transform_sample_counts(gutop, function(x) 100 * x/sum(x))
gutop<-filter_taxa(gutop, function(x) sum(x>0.01)>(.1*length(x)),TRUE)
summary(otu_table(gutop))

gutom<-subset_taxa(guto, Order == "Mycoplasmatales")
gutom = transform_sample_counts(gutom, function(x) 100 * x/sum(x))
summary(otu_table(gutom))

bello<-subset_samples(psnomito.noblnk.DESS.Cx, Body_cav == "Bell")
ps.orders <- bello %>%
  tax_glom("Phylum") %>%
  transmute_tax_table(Kingdom, Phylum, .otu = Phylum)
ps.orders = transform_sample_counts(ps.orders, function(x) 100 * x/sum(x))
top10os = names(sort(taxa_sums(ps.orders), TRUE)[1:10])
ps.oo=prune_taxa(top10os, ps.orders)
summary(otu_table(ps.oo))

bello = transform_sample_counts(bello, function(x) 100 * x/sum(x))
top10bos = names(sort(taxa_sums(bello), TRUE)[1:10])
ps.boo=prune_taxa(top10bos, bello)
summary(otu_table(ps.boo))
###########

wello<-subset_samples(psnomito.noblnk.DESS.Cx, Body_cav == "Sediment")
ps.orders <- wello %>%
  tax_glom("Phylum") %>%
  transmute_tax_table(Kingdom, Phylum, .otu = Phylum)
ps.orders = transform_sample_counts(ps.orders, function(x) 100 * x/sum(x))
top10os = names(sort(taxa_sums(ps.orders), TRUE)[1:10])
ps.oo=prune_taxa(top10os, ps.orders)
summary(otu_table(ps.oo))

wello<-subset_samples(psnomito.noblnk.DESS.Cx, Body_cav == "Sub")
wello = transform_sample_counts(wello, function(x) 100 * x/sum(x))
top10bos = names(sort(taxa_sums(wello), TRUE)[1:10])
ps.boo=prune_taxa(top10bos, wello)
summary(otu_table(ps.boo))




############ Top 200 OTUs in all samples together ######################
top200otus = names(sort(taxa_sums(psnomito.noblnk.DESS.P), TRUE)[1:200])
b10=prune_taxa(top200otus, psnomito.noblnk.DESS.P)
topperlabels<-plot_bar(b10,fill="Order")
topperGVC<- plot_bar(subset_samples(b10,Body_cav=="GVC"),x="Identity",fill="Order")+facet_grid(~Location,scales = "free_x")+guides(fill="none")+
  theme_linedraw()+xlab(NULL)+labs(title="A. GVC")+ylab("Relative abundance")
topperBell<- plot_bar(subset_samples(b10,Body_cav=="Bell"),x="Identity",fill="Order")+facet_grid(~Location,scales = "free_x")+guides(fill="none")+
  theme_linedraw()+xlab(NULL) +labs(title="B. Bell")+ylab("Relative abundance")
topperWater<- plot_bar(subset_samples(b10,Body_cav=="WTR"),x="Location",fill="Order")+guides(fill="none")+
  theme_linedraw()+xlab(NULL)+labs(title="C. Water")+ylab("Relative abundance")
topperSub<- plot_bar(subset_samples(b10,Body_cav=="Sub"),x="Location",fill="Order")+guides(fill="none")+
  theme_linedraw()+xlab(NULL)+labs(title="D. Benthos")
# Extract the legend. Returns a gtable
leg <- get_legend(topperlabels)
# Convert to a ggplot and print
legt<-as_ggplot(leg)
topperlow<-plot_grid(topperWater,topperSub,leg,nrow = 1,rel_widths = c(1.2,1.2,2))
topper<-plot_grid(topperGVC,topperBell,topperlow,ncol=1,rel_heights = c(3,3,2))
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/topper.png",plot = topper, dpi = 300, h = 18, w= 12, units = 'in')
############ Testing only, corallilyticus
only_endozoic_myco_vibrio<-subset_taxa(psnomito.noblnk.DESS.P,Genus=="Endozoicomonas"|Genus=="Mycoplasma"|Genus=="Vibrio")
write.csv(only_endozoic_myco_vibrio@refseq,file = "C:/Users/kmmuf/Downloads/Endo_Myco_Vib_refseq.csv")
emvGVC<-plot_bar(subset_samples(only_endozoic_myco_vibrio,Body_cav=="GVC"),fill="Species")
emvGVC
emv<-plot_bar(subset_taxa(only_endozoic_myco_vibrio,Species=="coralliilyticus"),y= "Abundance",x= "Body_cav", fill="Species")


psnomito.noblnk.DESS.P@sam_data$LocCav<-paste(psnomito.noblnk.DESS.P@sam_data$Body_cav,psnomito.noblnk.DESS.P@sam_data$Location)
corllips = merge_samples(psnomito.noblnk.DESS.P, "LocCav")
sample_data(corllips)$LocCav <- levels(sample_data(psnomito.noblnk.DESS.P)$LocCav)
corllips = transform_sample_counts(corllips, function(x) 100 * x/sum(x))
emv<-plot_bar(subset_taxa(corllips,Species=="coralliilyticus"),y= "Abundance", fill="Species")
fincorali<-
  emv +  theme_linedraw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ scale_fill_manual(values = "grey") +
  ggtitle("V. coralliilyticus by average of site and sample type") 

ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/emv.png",plot = fincorali, dpi = 300, h = 5, w= 5, units = 'in')
########### Starplot design ##############################
library(stars)
#Phylogenetic level of interest
Order.glom <- psnomito.noblnk.DESS.P %>%
  tax_glom("Order") %>%
  transmute_tax_table(Kingdom, Phylum, Class, Order, .otu = Order)
#Split by component
GVC.Order.glom<-subset_samples(Order.glom, Body_cav=="GVC")
Bell.Order.glom<-subset_samples(Order.glom, Body_cav=="Bell")
WTR.Order.glom<-subset_samples(Order.glom, Body_cav=="WTR")
Sub.Order.glom<-subset_samples(Order.glom, Body_cav=="Sub")
###Starplot production GVC
GVC.Order.glom.L = merge_samples(GVC.Order.glom, "Location")
sample_data(GVC.Order.glom.L)$Location <- levels(sample_data(GVC.Order.glom.L)$Location)
GVC.Order.glom.L = transform_sample_counts(GVC.Order.glom.L, function(x) 1 * x/sum(x))
GVC.Order.glom.L.data<-psotu2veg(GVC.Order.glom.L)
mns <- colMeans(GVC.Order.glom.L.data, na.rm=TRUE)
GVC.Order.glom.L.data.b <- as.data.frame(GVC.Order.glom.L.data[,order(mns,decreasing = TRUE)])
library(tidyverse)
#Making an "other" variable, reordering
GVC.Order.glom.L.data.b$Other<-1-(GVC.Order.glom.L.data.b$Pseudomonadales+GVC.Order.glom.L.data.b$Mycoplasmatales+GVC.Order.glom.L.data.b$Enterobacterales)
GVC.Order.glom.L.data.b<-GVC.Order.glom.L.data.b[,c(1:3,352,4:351)]
GVC.Order.glom.L.data.b<-GVC.Order.glom.L.data.b[c(3,2,1,8,7,6,4,5),]

###Starplot for bells
Bell.Order.glom.L = merge_samples(Bell.Order.glom, "Location")
sample_data(Bell.Order.glom.L)$Location <- levels(sample_data(Bell.Order.glom.L)$Location)
Bell.Order.glom.L = transform_sample_counts(Bell.Order.glom.L, function(x) 1 * x/sum(x))
Bell.Order.glom.L.data<-psotu2veg(Bell.Order.glom.L)
mns <- colMeans(Bell.Order.glom.L.data, na.rm=TRUE)
Bell.Order.glom.L.data.b <- as.data.frame(Bell.Order.glom.L.data[,order(mns,decreasing = TRUE)])
head(Bell.Order.glom.L.data.b[1:10])
Bell.Order.glom.L.data.b$Other<-1-(Bell.Order.glom.L.data.b$Pseudomonadales+Bell.Order.glom.L.data.b$Rhodobacterales+Bell.Order.glom.L.data.b$Chromatiales)
Bell.Order.glom.L.data.b<-Bell.Order.glom.L.data.b[,c(1:3,352,4:351)]
Bell.Order.glom.L.data.b<-Bell.Order.glom.L.data.b[c(3,2,1,8,7,6,4,5),]

##Starplot for WTR

WTR.Order.glom.L = transform_sample_counts(WTR.Order.glom, function(x) 1 * x/sum(x))
WTR.Order.glom.L.data<-psotu2veg(WTR.Order.glom.L)
mns <- colMeans(WTR.Order.glom.L.data, na.rm=TRUE)
WTR.Order.glom.L.data.b <- as.data.frame(WTR.Order.glom.L.data[,order(mns,decreasing = TRUE)])
head(WTR.Order.glom.L.data.b[1:10])
WTR.Order.glom.L.data.b$Other<-1-(WTR.Order.glom.L.data.b$Pseudomonadales+WTR.Order.glom.L.data.b$Rhodobacterales+WTR.Order.glom.L.data.b$Flavobacteriales)
WTR.Order.glom.L.data.b<-WTR.Order.glom.L.data.b[,c(1:3,352,4:351)]
WTR.Order.glom.L.data.b<-WTR.Order.glom.L.data.b[c(3,2,1,8,7,6,4,5),]
##Starplot for Benthos

Sub.Order.glom.L = transform_sample_counts(Sub.Order.glom, function(x) 1 * x/sum(x))
Sub.Order.glom.L.data<-psotu2veg(Sub.Order.glom.L)
mns <- colMeans(Sub.Order.glom.L.data, na.rm=TRUE)
Sub.Order.glom.L.data.b <- as.data.frame(Sub.Order.glom.L.data[,order(mns,decreasing = TRUE)])
head(Sub.Order.glom.L.data.b[1:10])
Sub.Order.glom.L.data.b$Other<-1-(Sub.Order.glom.L.data.b$Desulfobacterales+Sub.Order.glom.L.data.b$Rhodobacterales+Sub.Order.glom.L.data.b$Flavobacteriales)
Sub.Order.glom.L.data.b<-Sub.Order.glom.L.data.b[,c(1:3,352,4:351)]
Sub.Order.glom.L.data.b<-Sub.Order.glom.L.data.b[c(3,2,1,8,7,6,4,5),]


###Star plot plotting for export
#Note, key location shifts if labels are "true"

par(bg = NA,mfrow=c(4,1))

stars(GVC.Order.glom.L.data.b[, 1:4], len = 1.5, key.loc = c(2.1, 2.1), nrow = 1, ncol = 10,
      radius = TRUE, scale= FALSE, lwd = 2, flip.labels = FALSE, labels = NULL, 
      key.labels = c(substr(colnames(GVC.Order.glom.L.data.b[1:4]), 1,3)),
      cex=0.5, 
      col.lines = c("black","black","black","black","black","black","black","black"),
      col.stars = c("coral","coral","coral","coral","coral","coral","coral","coral","coral"))




stars(Bell.Order.glom.L.data.b[, 1:4], len = 1.5, key.loc = c(2.1, 2.1), nrow = 1, ncol = 10,
      radius = TRUE, scale= FALSE, flip.labels = FALSE,
      key.labels = c(substr(colnames(Bell.Order.glom.L.data.b[1:4]), 1,3)),
      cex=.5, lwd = 2,
      col.lines = c("black","black","black","black","black","black","black","black"),
      col.stars = c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen"),
      labels = NULL)

stars(WTR.Order.glom.L.data.b[, 1:4], len = 1.5, key.loc = c(2.1, 2.1), nrow = 1, ncol = 10,
      radius = TRUE, scale= FALSE, flip.labels = FALSE, 
      key.labels = c(substr(colnames(WTR.Order.glom.L.data.b[1:4]), 1,3)),
      cex=.5, lwd = 2,
      col.lines = c("black","black","black","black","black","black","black","black"),
      col.stars = c("lightblue","lightblue","lightblue","lightblue","lightblue","lightblue","lightblue","lightblue","lightblue"),
      labels = NULL)

stars(Sub.Order.glom.L.data.b[, 1:4], len = 1.5, key.loc = c(2.1, 2.1), nrow = 1, ncol = 10,
      radius = TRUE, scale= FALSE, flip.labels = FALSE, 
      key.labels = c(substr(colnames(Sub.Order.glom.L.data.b[1:4]), 1,3)),
      cex=.5, lwd = 2,
      col.lines = c("black","black","black","black","black","black","black","black"),
      col.stars = c("brown4","brown4","brown4","brown4","brown4","brown4","brown4","brown4","brown4"),
      labels = NULL)

####### Aldex2 packages##########################
library(tidyr)
library(tidyverse)
library(kableExtra)
library(ggplot2)
library(patchwork)
library(dplyr)
library(lubridate)
library(vegan)
library(ggpubr)
library(rstatix)
library(stringr)
library(emmeans)
library(multcomp)
library(multcompView)
library(dada2)
library(phyloseq)
library(ggh4x)
library(grid)
library(ALDEx2)
library(ggrepel)
################ Create body cavity comparisons (ALDEX2)#################
ps_comp <- subset_samples(psnomito.noblnk.DESS.Cx, Body_cav == "GVC"| Body_cav == "Bell")
##
pf<-ps_comp #subset phyloseq object here, you can copy and paste the next few sections to easily repeat - just change output and input names
conds <- sample_data(pf)$Body_cav
reads<-t(data.frame(otu_table(pf)))
reads<-reads[rowSums(reads) > 0,] #removes rows with only 0 in dataframe

#run clr
x <- aldex.clr(reads, conds, mc.samples=128, denom="all", verbose=F)

#run t-test because we are only comparing one variable with two levels
x.tt <- aldex.ttest(x, paired.test=TRUE, verbose=FALSE)
x.kw <- aldex.kw(x)

#run effect
x.effect <- aldex.effect(x, CI=T, verbose=FALSE, paired.test=TRUE)

#combine table
x.all <- data.frame(x.tt,x.kw,x.effect)
# Plot 
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="glm", cutoff.pval=0.05) #can use the glm model test as the sig value 
aldex.plot(x.all, type="MW", test="glm",cutoff.pval=0.05)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", cutoff.pval=0.05) #can use the glm model test as the sig value 
aldex.plot(x.all, type="MW", test="welch",cutoff.pval=0.05)

# Check for asymmetry - should be centered around 0
par(mfrow=c(1,2))
hist(x.all$diff.btw)
hist(x.all$effect)
######## Plot effect
par(mfrow=c(1,2))
plot(x.all$effect, x.all$glm.eBH, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Effect size", ylab=" BH P value", main="Effect size plot")
points(x.all$effect, x.all$we.eBH, cex=0.7, col='red', pch=19)
abline(h=0.05, lty=2, col="grey")
abline(v=2, lty=2, col="grey")
abline(v=-2, lty=2, col="grey")

plot(x.all$diff.btw, x.all$glm.eBH, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Difference", ylab="P value", main="Volcano plot")
points(x.all$effect, x.all$we.eBH, cex=0.7, col='red', pch=19)
abline(h=0.05, lty=2, col="grey")

#subset ASVs with an effect size > 2
#pvalue for glm < 0.05 and ttest pvalue < 0.05
#nr. sigtaxa

x.all.comp1<-x.all
x.all.sig.comp1<-x.all.comp1 %>% filter(effect > 2 | effect < -2) %>% filter(glm.eBH < 0.05) %>% filter(we.eBH < 0.05)
x.all.sig.comp1
sig_features_comp1<-rownames(x.all.sig.comp1)
length(sig_features_comp1)

#volcano plot for Differential abundance

aldplot<-ggplot(x.all.comp1, aes(y=-log10(we.eBH), x=effect)) + 
  geom_point() + 
  geom_point(data=x.all.sig.comp1, aes(y=-log10(we.eBH), x=effect),color='blue') + 
  geom_hline(yintercept = -log10(0.05),color='gray', linetype='dashed') + 
  geom_vline(xintercept = 2, color='gray', linetype='dashed') +  geom_vline(xintercept = -2,color='gray', linetype='dashed')+
  geom_text_repel(data=x.all.sig.comp1,aes(y=-log10(we.eBH), x=effect,label=rownames(x.all.sig.comp1)), min.segment.length=0, max.overlaps = 15) + 
  theme_bw() + ggtitle(label = 'Aldex2 plot of ASV differential 
  abundance between Bell vs. GVC') + xlab("Effect size")

ggsave(file="C:/Users/kmmuf/Documents/Cass_Mic_seq_test/aldex_dif_plot.png",plot= aldplot, h=3,w=4,dpi = 300, units = "in")

########## Core microbiome discovery ####################################
library(vegan)
library(phyloseq)
library(ggrepel)
library(patchwork)
library(microbiome)
library(RColorBrewer)
# subset for control group
ps.GVC<-subset_samples(psnomito.noblnk.DESS.Cx, Body_cav=='GVC')
ps.bell<-subset_samples(psnomito.noblnk.DESS.Cx, Body_cav=='Bell')

#removes taxa with 0 counts
ps.GVC <- prune_taxa(taxa_sums(ps.GVC) > 30, ps.GVC)
ps.GVC
#2507 taxa
ps.bell <- prune_taxa(taxa_sums(ps.bell) > 30, ps.bell)
ps.bell
#8555 taxa

pseq.GVC.rel <- microbiome::transform(ps.GVC, "compositional")
pseq.bell.rel <- microbiome::transform(ps.bell, "compositional")


# Core with compositionals (check popu:
head(prevalence(pseq.GVC.rel, detection = 0.1/100, sort = TRUE),20)
#Thresh 70%+ ASV2, 1
#50%+ ASV 4, 11
#ASV1      ASV4      ASV2     ASV11      ASV6 
#0.9000000 0.8000000 0.8000000 0.6000000 0.5000000
head(prevalence(pseq.bell.rel, detection = 0.1/100, sort = TRUE),20)
#ASV1    ASV104     ASV86     ASV95      ASV6 
#0.7333333 0.7000000 0.7000000 0.6666667 0.6666667 
# ASV135    ASV105      ASV4
# 0.6333333 0.6333333 0.6000000
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
# Only include core taxa

pseq.core.GVC <- core(pseq.GVC.rel, detection =0, prevalence = 70/100)
pseq.core.bell <- core(pseq.bell.rel, detection =0, prevalence = 70/100)

colnames(pseq.core.GVC@tax_table)
colnames(pseq.core.GVC@tax_table)<-c("Domain", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")
colnames(pseq.core.bell@tax_table)<-c("Domain", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")

pseq.core.GVC <- microbiome::add_besthit(pseq.core.GVC)

pseq.core.bell <- microbiome::add_besthit(pseq.core.bell)

p <- plot_core(pseq.core.GVC, plot.type = "heatmap",
               colours = rev(brewer.pal(5, "Spectral")),           
               prevalences = prevalences,
               detections = detections
) +  theme_linedraw()+
  labs(x = "Detection Threshold",y= "ASV") + ggtitle("Gastrovascular Cavity",subtitle = "Sample/Abundance minimum matrix")
pbell <- plot_core(pseq.core.bell, plot.type = "heatmap",
               colours = rev(brewer.pal(5, "Spectral")),           
               prevalences = prevalences,
               detections = detections
) +
  labs(x = "Detection Threshold",y= "ASV") + ggtitle("Bell",subtitle = "Sample/Abundance minimum matrix") + theme_linedraw()
pbell
gvcandbell<-plot_grid(p,pbell,ncol=1,rel_heights = c(1,1.2))
ggsave(filename = "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/gvc_bell_core.png", plot= gvcandbell, h = 7, w = 7, dpi = 300, units = "in")
####
core.taxa <- core_members(pseq.bell.rel, detection = 0, prevalence = 70/100)
core.taxa2 <- core_members(pseq.GVC.rel, detection = 0, prevalence = 70/100)
core.complete<-c("ASV1", "ASV2", "ASV4", "ASV6", "ASV11",  "ASV37",  "ASV86","ASV95",  "ASV104", "ASV105", "ASV132", "ASV135")
ps.core<-subset_taxa(psnomito.noblnk.DESS.P, rownames(tax_table(psnomito.noblnk.DESS.P)) %in% core.taxa2) #subsets relative abundance table from the full dataset 
ps.core@tax_table
# Taxonomy Table:     [ 5 taxa by 7 taxonomic ranks ]:
#   Kingdom  Phylum         Class               Order           Family Genus Species
# <chr>    <chr>          <chr>               <chr>           <chr>  <chr> <chr>  
#   1 ASV1   Bacteria Proteobacteria Gammaproteobacteria Pseudomonadales Endoz… Endo… atrinae
# 2 ASV2   Bacteria Firmicutes     Bacilli             Mycoplasmatales Mycop… Myco… <NA>   
#   3 ASV4   Bacteria Proteobacteria Gammaproteobacteria Pseudomonadales Endoz… Endo… atrinae
# 4 ASV11  Bacteria Proteobacteria Gammaproteobacteria Enterobacteral… Vibri… Vibr… corall…
# 5 ASV132 Bacteria Proteobacteria Gammaproteobacteria Enterobacteral… Vibri… Vibr… hepata…
colnames(ps.core@tax_table)<-c("Domain", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")
ps.core <- microbiome::add_besthit(ps.core)

GVCcore<-plot_bar(subset_samples(ps.core,Body_cav=="GVC"|Body_cav=="Bell"), x= "Identity", fill="OTU")+facet_grid(~Body_cav+reorder(Location, Lat), scales = "free", space = "free")+
  #scale_fill_manual(values = c("deepskyblue1","palegreen1","palegreen3","salmon", "deepskyblue3"))+
  guides(fill="none")+
  theme_linedraw()+xlab("Bell and GVC samples organized by site and body cavity South to North")+ ylab("Abundance (0 - 100%)")
GVCcoreinsubs<-plot_bar(subset_samples(ps.core,Body_cav=="Sub"|Body_cav=="WTR"), x= "Location", fill="OTU")+facet_grid(~Body_cav+reorder(Location, Lat), scales = "free", space = "free")+
  #scale_fill_manual(values = c("deepskyblue1","palegreen1","palegreen3","salmon", "deepskyblue3"))+
  theme_linedraw()+xlab("Environmental samples organized by site South to North") + ylab("Abundance (0 - 0.5%)")
GVCplotgrid<-plot_grid(GVCcore,GVCcoreinsubs, rel_heights = c(1,.5), ncol = 1)
GVCplotgrid
ggsave(filename = "C:/Users/kmmuf/Documents/Cass_Mic_seq_test/GVCplotgrid.png", plot= GVCplotgrid, h = 6, w = 10, dpi = 300, units = "in")

################ individualized linear models ###################
#(no significance)
library(lmtest)
library(ggfortify)

core.60_01_GVC<-c("ASV1", "ASV2", "ASV4", "ASV11")
ps.core2<-subset_taxa(subset_samples(psnomito.noblnk.DESS.P, Body_cav == "GVC"), rownames(tax_table(psnomito.noblnk.DESS.P)) %in% core.60_01_GVC) #subsets relative abundance table from the full dataset 
core4_melted<-psmelt(ps.core2)
core4_melted$Time<- gsub("[:]", "" , core4_melted$Time, perl=TRUE)
core4_melted$Time<-as.numeric(core4_melted$Time)
core4_melted_ASV1<-subset(core4_melted,OTU == "ASV1")
core4_melted_ASV2<-subset(core4_melted,OTU == "ASV2")
core4_melted_ASV4<-subset(core4_melted,OTU == "ASV4")
core4_melted_ASV11<-subset(core4_melted,OTU == "ASV11")
colnames(core4_melted)
# [1] "OTU"                     "Sample"                 
# [3] "Abundance"               "Sample_NAME"            
# [5] "Location"                "Identity"               
# [7] "Body_cav"                "Buffer"                 
# [9] "Species_HOST"            "Temp"                   
# [11] "Specimen"                "Lat"                    
# [13] "Long"                    "Salinity"               
# [15] "pH"                      "Temperature_.C."        
# [17] "Date"                    "Species_O"              
# [19] "Time"                    "Diameter..cm."          
# [21] "OAL..cm."                "Depth..m."              
# [23] "Concentration..indv.m2." "Neg"                    
# [25] "LocCav"                  "Kingdom"                
# [27] "Phylum"                  "Class"                  
# [29] "Order"                   "Family"                 
# [31] "Genus"                   "Species" 
###
library(corrplot)
corrplot(cor(core4_melted_ASV1[,c(10,12,14:16,19:22)]),method = "number")
core4_melted_ASV1<-core4_melted_ASV1[-24,]
#Depth and salinity and temp correlated (.7)
histogram(core4_melted_ASV1$Abundance)
histogram(sqrt(core4_melted_ASV1$Abundance))
histogram(log(core4_melted_ASV1$Abundance))
#Data cannot be made normal easily, but is reasonably sized n


ASVall1<-lm(Abundance~Salinity+Temp+pH+Diameter..cm.+Lat+Depth..m., data=core4_melted_ASV1)
ASV1_pd<-lm(Abundance~pH*Depth..m., data=core4_melted_ASV1)
ASV1_d<-lm(Abundance~Depth..m., data=core4_melted_ASV1)
summary(ASVall1)
summary(ASV1_pd)
summary(ASV1_d)
anova(ASVall1,ASV1_pd,ASV1_d)
par(mfrow = c(2,2))
plot(ASV1_pd)
#overall, not terrible
ggplot(core4_melted_ASV1,aes(x=Depth..m., y=Abundance, fill = pH))+geom_point(aes(color = pH))
#Not acceptable
########ASV2
core4_melted_ASV2<-core4_melted_ASV2[-24,]
histogram(core4_melted_ASV2$Abundance)
histogram(sqrt(core4_melted_ASV2$Abundance))
histogram(log(core4_melted_ASV2$Abundance))
#Data is reasonable when sqrt transformed


ASVall2<-lm(sqrt(Abundance)~Salinity+Temp+pH+Diameter..cm.+Lat+Depth..m., data=core4_melted_ASV2)
ASV2_<-lm(Abundance~Lat, data=core4_melted_ASV2)
summary(ASVall2)
summary(ASV2_)
#No explanatory value

###### ASV4
core4_melted_ASV4<-core4_melted_ASV1[-24,]
histogram(core4_melted_ASV1$Abundance)
histogram(sqrt(core4_melted_ASV1$Abundance))
histogram(log(core4_melted_ASV1$Abundance))
ASVall4<-lm(Abundance~Salinity+Temp+pH+Diameter..cm.+Lat+Depth..m., data=core4_melted_ASV4)
ASV4_tpd<-lm(Abundance~Temp*pH*Depth..m., data = core4_melted_ASV4)
ASV4_td<-lm(Abundance~Temp*Depth..m., data = core4_melted_ASV4)

#### Temp pH and depth covary. Selecting Depth
ASV4_d<-lm(Abundance~Depth..m., data= core4_melted_ASV4)
summary(ASVall4)
summary(ASV4_tpd)
summary(ASV4_td)
summary(ASV4_d)
anova(ASVall4,ASV4_tpd,ASV4_td,ASV4_d)
par(mfrow = c(2,2))
plot(ASV4_tpd)
#Residuals look ok
ggplot(core4_melted_ASV4,aes(x=Depth..m., y=Abundance, fill = Temp))+geom_point(aes(color = pH))
# Poor explanatory value

####ASV11
core4_melted_ASV11<-core4_melted_ASV11[-24,]
histogram(core4_melted_ASV11$Abundance)
histogram(sqrt(core4_melted_ASV11$Abundance))
histogram(log(core4_melted_ASV11$Abundance))
ASVall11<-lm(log(Abundance+0.00001)~Salinity+Temp+pH+Diameter..cm.+Lat+Depth..m., data=core4_melted_ASV11)
summary(ASVall11)
#No explanatory value

#############Cassiopea_andromeda_and_xam
psnomito.noblnk.DESS.GVCBELL.CACX<-subset_samples(psnomito.noblnk.DESS, Body_cav=="GVC"|Body_cav=="Bell")
psnomito.noblnk.DESS.GVCBELL.CACX<-subset_samples(psnomito.noblnk.DESS.GVCBELL.CACX, Location=="CK"|Location=="KL"|Location=="GB")
psnomito.noblnk.DESS.GVCBELL.CACX.1<-subset_taxa(psnomito.noblnk.DESS.GVCBELL.CACX, taxa_sums(psnomito.noblnk.DESS.GVCBELL.CACX)>1)
psnomito.noblnk.DESS.GVCBELL.CACX.1= transform_sample_counts(psnomito.noblnk.DESS.GVCBELL.CACX.1, function(x) 1 * x/sum(x))
anga<-plot_bar(psnomito.noblnk.DESS.GVCBELL.CACX.1,x="Identity", y="Abundance", fill="Phylum") + facet_grid(Species_HOST~Body_cav, scales = "free")
ggsave("C:/Users/kmmuf/Downloads/anga.png", anga, h=7, w= 20, units="in")

eucDist.CA <- phyloseq::distance(psnomito.noblnk.DESS.GVCBELL.CACX.1, method = "euclidean")
# Calculate ordination
eucMDS.CA  <- ordinate(psnomito.noblnk.DESS.GVCBELL.CACX.1, "PCoA", distance=eucDist.CA)
## Make plot
# Create plot, store as temp variable, p
gare.CA <- plot_ordination(psnomito.noblnk.DESS.GVCBELL.CACX.1, eucMDS.CA, color="Species_HOST",shape="Body_cav")
garecom.CA<-gare.CA + theme_linedraw()+labs(shape="Sample Type",color = "Species")+ggtitle("Euclidean distance by sample type and species")+
  scale_colour_manual(values = c("grey", "black"))+geom_point(size=2)
garecom.CA
garecomseg<-garecom+facet_wrap(~Body_cav)+ggtitle("Euclidean distance split by sample type")+
  guides(color="none",shape="none")
gare_tog<-plot_grid(garecom,garecomseg,ncol=1)
