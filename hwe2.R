load("./data/pheno.RData")

library('data.table')
library('dplyr')
library('HardyWeinberg')
source("funs.R")
set.seed(1234)

hwe_wrapper <- function(x){
    counts <- c("AA"=sum(x==0), "BB"=sum(x==2), "AB"=sum(x==1))
    result <- tryCatch({
        fit <- HWChisq(counts, cc=0)
        res <- data.frame(chisq=fit$chisq, pval=fit$pval)
        res <- cbind(res, t(counts))
        return(res)
    }, warning = function(warn_condition) {
        res <- data.frame(chisq=NA, pval=NA)
        res <- cbind(res, t(counts))
        return(res)
    }, error = function(error_condition) {
        res <- data.frame(chisq=NA, pval=NA)
        res <- cbind(res, t(counts))
        return(res)
    })

    return(result)
}

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

# select fields, drop missing data and round genotypes
dat <- dat %>% 
    dplyr::select(chr2_136608646_G_A, ethnic_background.21000.0.0) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(chr2_136608646_G_A = round(chr2_136608646_G_A))

# estimate HWE test statistics without Yates correction
results <- dat %>% 
    dplyr::group_by(ethnic_background.21000.0.0) %>%
    dplyr::summarize(hwe_wrapper(chr2_136608646_G_A))

# append coding
coding <- fread("http://biobank.ctsu.ox.ac.uk/showcase/codown.cgi?id=1001&btn_glow=Download")
coding <- coding %>% dplyr::select("coding", "meaning") %>% 
    dplyr::rename(Population = "meaning")
results <- merge(results, coding, by.x="ethnic_background.21000.0.0", by.y="coding")

# add N and AF
results <- results %>% 
    dplyr::group_by(ethnic_background.21000.0.0) %>%
    dplyr::mutate(N=sum(AA, AB, BB), AF=((AB + (BB * 2)) /2 ) / (AA+AB+BB))

# update populations
results$Population[which(results$Population == "British")] <- "White, British"
results$Population[which(results$Population == "Irish")] <- "White, Irish"
results <- results[results$Population != "White",]
results <- results[results$Population != "Mixed",]
results <- results[results$Population != "Asian or Asian British",]
results <- results[results$Population != "Black or Black British",]

write.csv(results, file="data/ukbb-hwe.csv", row.names =F)

# gnomad
gnomad <- fread("data/gnomad.txt", fill=T)
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
    l <- c(AA=gnomad$GG[i], BB=gnomad$AA[i], AB=gnomad$GA[i])
    N <- sum(l)
    AF <- (gnomad$GA[i] + gnomad$AA[i] * 2) * 0.5 / N
    fit <- tryCatch(HWChisq(l, cc=0),error=function(e) NA, warning=function(w) NA)
    if (is.na(fit)){
        gnomad_hwe <- rbind(gnomad_hwe, data.frame(population=gnomad$Population[i], chisq=NA, pval=NA, N, AF, GA=gnomad$GA[i], AA=gnomad$AA[i], GG=gnomad$GG[i]))
    } else {
        gnomad_hwe <- rbind(gnomad_hwe, data.frame(population=gnomad$Population[i], chisq=fit$chisq, pval=fit$pval, N, AF, GA=gnomad$GA[i], AA=gnomad$AA[i], GG=gnomad$GG[i]))
    }
}