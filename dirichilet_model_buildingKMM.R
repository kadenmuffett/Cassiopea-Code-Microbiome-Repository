## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
## --- Install-----

BiocManager::install("DirichletMultinomial")

install.packages("lattice")
install.packages("parallel")


## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----library, message = FALSE-------------------------------------------------
library(DirichletMultinomial)
library(lattice)
library(parallel)
library(microbiome)

## ----colors---------------------------------------------------------
options(width=70, digits=2)
full <- TRUE
.qualitative <- DirichletMultinomial:::.qualitative

## ----data-input-----------------------------------------------------
psnomito.noblnk.DESS.Cx_for_dmm<-psnomito.noblnk.DESS.Cx
colnames(psnomito.noblnk.DESS.Cx_for_dmm@tax_table)<-c("Domain", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")

psnomito.noblnk.DESS.Cx_for_dmm <- microbiome::add_besthit(psnomito.noblnk.DESS.Cx_for_dmm)
count <- psotu2veg(prune_taxa(taxa_sums(psnomito.noblnk.DESS.Cx_for_dmm) > 300, psnomito.noblnk.DESS.Cx_for_dmm))
count[1:5, 1:3]

## ----taxon-counts---------------------------------------------------
cnts <- log10(colSums(count))
range(cnts)
dev.off()

densityplot(
    cnts, xlim=range(cnts),
    xlab="Taxon representation (log 10 count)"
)

## ----fit------------------------------------------------------------
if (full) {
    fit <- mclapply(1:7, dmn, count=count, verbose=TRUE)
    save(fit, file=file.path(tempdir(), "fit.rda"))
} else data(fit)
fit[[4]]

## ----min-laplace, figure=TRUE---------------------------------------
lplc <- sapply(fit, laplace)
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
(best <- fit[[which.min(lplc)]])

## ----mix-weight-----------------------------------------------------
mixturewt(best)
head(mixture(best), 3)

## ----fitted---------------------------------------------------------
splom(log(fitted(best)))

## ----posterior-mean-diff--------------------------------------------
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p4 <- fitted(best, scale=TRUE)
colnames(p4) <- paste("m", 1:2, sep="")
(meandiff <- colSums(abs(p4 - as.vector(p0))))
sum(meandiff)

## ----table-1--------------------------------------------------------
diff <- rowSums(abs(p4 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- cbind(Mean=p0[o], p4[o,], diff=diff[o], cdiff)
DT::datatable(df) |>
    DT::formatRound(colnames(df), digits = 4)

## ----heatmap-similarity---------------------------------------------
heatmapdmn(count, fit[[1]], best, 30)

## ----twin-pheno-----------------------------------------------------
pheno0 <- sample_data(psnomito.noblnk.DESS.Cx)$Body_cav
names(pheno0) <- rownames(count)
table(pheno0)
pheno<-pheno0
table(pheno)

## ----subsets--------------------------------------------------------
counts <- lapply(levels(pheno), csubset, count, pheno)
sapply(counts, dim)
keep <- c("Bell", "GVC")
count <- count[pheno %in% keep,]
pheno <- factor(pheno[pheno %in% keep], levels=keep)

## ----fit-several----------------------------------------------------
if (full) {
    bestgrp <- dmngroup(
        count, pheno, k=1:5, verbose=TRUE, mc.preschedule=FALSE
    )
    save(bestgrp, file=file.path(tempdir(), "bestgrp.rda"))
} else data(bestgrp)

## ----best-several---------------------------------------------------
bestgrp
lapply(bestgrp, mixturewt)
c(
    sapply(bestgrp, laplace),
    'Bell+GVC' = sum(sapply(bestgrp, laplace)),
    Single = laplace(best)
)
########
heatmapdmn(count, fit[[1]], best, 30)

##########
library(reshape2)
library(magrittr)
library(dplyr)
for (k in seq(ncol(fitted(bestgrp$Bell)))) {
  d <- melt(fitted(bestgrp$Bell))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.995))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() + ylab("Relative Community Assignment Value")+xlab("ASV")+
    labs(title = paste("Top drivers of DMM: Bell community"))
  print(j)
}

DDMdrivers<-cowplot::plot_grid(j, p,nrow = 2)
DDMdrivers
ggsave("C:/Users/kmmuf/Documents/Cass_Mic_seq_test/DMM_loadings.png", plot = DDMdrivers, h = 4.5, w=6, dpi=300, units="in")
## ----confusion------------------------------------------------------
xtabs(~pheno + predict(bestgrp, count, assign=TRUE))
# pheno  Bell GVC
# Bell   27   3
# GVC     2  28
## ----cross-validate-------------------------------------------------
if (full) {
    ## full leave-one-out; expensive!
    xval <- cvdmngroup(
        nrow(count), count, c(Bell=1, GVC=3), pheno,
        verbose=TRUE, mc.preschedule=FALSE
    )
    save(xval, file=file.path(tempdir(), "xval.rda"))
} else data(xval)

## ----ROC-dmngroup---------------------------------------------------
bst <- roc(pheno[rownames(count)] == "GVC",
predict(bestgrp, count)[,"GVC"])
bst$Label <- "Single"
two <- roc(pheno[rownames(xval)] == "GVC", xval[,"GVC"])
two$Label <- "Two group"
both <- rbind(bst, two)
pars <- list(superpose.line=list(col=.qualitative[1:2], lwd=2))
xyplot(
    TruePostive ~ FalsePositive, group=Label, both,
    type="l", par.settings=pars,
    auto.key=list(lines=TRUE, points=FALSE, x=.6, y=.1),
    xlab="False Positive", ylab="True Positive"
)

## ----sessionInfo----------------------------------------------------
sessionInfo()

