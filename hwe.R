load("./data/pheno.RData")
library('data.table')
library('HardyWeinberg')
source("funs.R")
set.seed(1234)

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=FALSE, drop_related=TRUE)

# load covariates & snp
cov <- get_covariates()
pc <- get_genetic_principal_components()
snp <- extract_variant_from_bgen("2", 136608646, "G", "A")

# merge data
dat <- merge(linker, pheno, by.x="app16729", by.y="eid")
dat <- merge(dat, cov, "appieu")
dat <- merge(dat, pc, "appieu")
dat <- merge(dat, snp, "appieu")

# select fields and drop missing data
t <- dat[,c("chr2_136608646_G_A", "ethnic_background.21000.0.0")]
t <- na.omit(t)

# HWE by ethnicity
run_hwe_pval <- function(x){
    result = tryCatch({
        t<-table(round(x))
        l<-c(AA=as.integer(t[names(t)=="0"]), BB=as.integer(t[names(t)=="2"]), AB=as.integer(t[names(t)=="1"]))
        HWChisq(l, cc=0)
        return(fit$pval)
    }, warning = function(warn_condition) {
        return(NA)
    }, error = function(error_condition) {
        return(NA)
    })
}
run_hwe_chi <- function(x){
    result = tryCatch({
        t<-table(round(x))
        l<-c(AA=as.integer(t[names(t)=="0"]), BB=as.integer(t[names(t)=="2"]), AB=as.integer(t[names(t)=="1"]))
        fit <- tryCatch(HWChisq(l, cc=0),error=function(e) NA, warning=function(w) NA)
        if (is.na(fit)){
            return(NA)
        } else {
            return(fit$chisq)
        }
    }, error = function(error_condition) {
        return(NA)
    })
}
run_hwe_n <- function(x){
    result = tryCatch({
        return(length(x))
    }, error = function(error_condition) {
        return(NA)
    })
}
run_gt_n <- function(x){
    return(c(CC=sum(x == 0), CT=sum(x == 1), TT=sum(x == 2)))
}

# test
hwe.p <- aggregate(t$chr2_136608646_G_A, by=list(t$ethnic_background.21000.0.0), run_hwe_pval)
names(hwe.p)[2]<-'pval'
hwe.chi <- aggregate(t$chr2_136608646_G_A, by=list(t$ethnic_background.21000.0.0), run_hwe_chi)
names(hwe.chi)[2]<-'chi'
hwe.n <- aggregate(t$chr2_136608646_G_A, by=list(t$ethnic_background.21000.0.0), run_hwe_n)
names(hwe.n)[2]<-'n'
gt.counts <- aggregate(round(t$chr2_136608646_G_A), by=list(t$ethnic_background.21000.0.0), run_gt_n)

# tidy table
res <- merge(hwe.p, hwe.chi, "Group.1")
res <- merge(res, hwe.n, "Group.1")
res <- merge(res,gt.counts, "Group.1")
coding <- fread("http://biobank.ctsu.ox.ac.uk/showcase/codown.cgi?id=1001&btn_glow=Download")
coding <- coding[,c("coding", "meaning")]
names(coding)[2]<- "Population"
res <- merge(res, coding, by.x="Group.1", by.y="coding")
res <- cbind(res, as.data.frame(res$x))
res$x <- NULL
res$af <- ((res$CT + (res$TT * 2)) /2 ) / res$n
# round values
res$af <- sapply(res$af, tbl_round)
res$pval <- sapply(res$pval, tbl_round)
res$chi <- sapply(res$chi, tbl_round)
res <- res[,c("Population", "af", "n", "TT", "CT", "CC", "pval", "chi")]
names(res) <- c("Population", "AF", "N", "AA", "GA", "GG", "P", "x^2")
res <- res[order(res$N, decreasing=T),]
res$Population[which(res$Population == "British")] <- "White, British"
res$Population[which(res$Population == "Irish")] <- "White, Irish"
res <- res[res$Population != "White",]
res <- res[res$Population != "Mixed",]
res <- res[res$Population != "Asian or Asian British",]
res <- res[res$Population != "Black or Black British",]

# gnomAD HWE
gnomad <- fread("gnomad.txt", fill=T)
gnomad$AA <- gnomad$Homozygotes
gnomad$AA_AC <- gnomad$AA * 2
gnomad$GA <- gnomad[["Allele count"]] - gnomad$AA_AC
gnomad$GG <- gnomad[["Allele number"]] / 2 - (gnomad$AA + gnomad$GA)
(gnomad$GA + gnomad$AA * 2) / (gnomad$GA *2 + gnomad$AA *2 + gnomad$GG *2)
all(gnomad$GA + gnomad$AA * 2 == gnomad[["Allele count"]])
all(gnomad$GA*2 + gnomad$AA*2 +gnomad$GG*2  == gnomad[["Allele number"]])
gnomad$AA_AC <- NULL

gnomad_hwe <- data.frame()
for (i in 1:nrow(gnomad)){
    l<-c(AA=gnomad$GG[i], BB=gnomad$AA[i], AB=gnomad$GA[i])
    fit <- tryCatch(HWChisq(l, cc=0),error=function(e) NA, warning=function(w) NA)
    if (is.na(fit)){
        gnomad_hwe <- rbind(gnomad_hwe, data.frame(population=gnomad$Population[i], chisq=NA, pval=NA))
    } else {
        gnomad_hwe <- rbind(gnomad_hwe, data.frame(population=gnomad$Population[i], chisq=fit$chisq, pval=fit$pval))
    }
}
