load("./data/pheno.RData")

library('data.table')
library('broom')
library('optparse')
library('jsonlite')
library('httr')
library('TwoSampleMR')
library('rbgen')
library('stringr')
source("funs.R")
set.seed(1234)

# load data
cov <- get_covariates()
pc <- get_genetic_principal_components()
snp <- extract_variant_from_bgen("2", 136608646, "G", "A")
linker <- get_filtered_linker(drop_standard_excl=FALSE, drop_non_white_british=FALSE, drop_related=FALSE)
dat <- merge(linker, pheno, by.x="app16729", by.y="eid")
dat <- merge(dat, cov, "appieu")
dat <- merge(dat, pc, "appieu")
dat <- merge(dat, snp, "appieu")

# load spouse list and restrict to couples
spouses <- fread("./Spouse-pairs-IEU.txt")
names(spouses) <- c("appieu", "couple")
spouses <- merge(spouses, dat, "appieu")
spouses <- spouses %>% 
    dplyr::select("couple", "chr2_136608646_G_A", "sex.31.0.0", "age_at_recruitment.21022.0.0", "milk_type_used.1418.0.0", paste0("PC",seq(1,40)))
spouses <- na.omit(spouses)

# select one partner and compare with another partner
female <- spouses %>% dplyr::filter(sex.31.0.0 == 0)
male <- spouses %>% dplyr::filter(sex.31.0.0 == 1)
names(female) <- paste0("female.", names(female))
names(male) <- paste0("male.", names(male))
spouses <- merge(female, male, by.x="female.couple", by.y="male.couple")

# additive test for association of LP genotype between mates
f <- as.formula(paste0("female.chr2_136608646_G_A ~ male.chr2_136608646_G_A + female.age_at_recruitment.21022.0.0 + male.age_at_recruitment.21022.0.0 + ", paste0("female.PC",seq(1,40), collapse=" + "), " + ", paste0("male.PC",seq(1,40),  collapse=" + ")))
fit <- lm(f, data=spouses)

# dominant test for association of LP genotype between mates
spouses$male.chr2_136608646_G_A.bin <- round(spouses$male.chr2_136608646_G_A) > 0
spouses$female.chr2_136608646_G_A.bin <- round(spouses$female.chr2_136608646_G_A) > 0
f <- as.formula(paste0("female.chr2_136608646_G_A.bin ~ male.chr2_136608646_G_A.bin + female.age_at_recruitment.21022.0.0 + male.age_at_recruitment.21022.0.0 + ", paste0("female.PC",seq(1,40), collapse=" + "), " + ", paste0("male.PC",seq(1,40),  collapse=" + ")))
fit <- glm(f, data=spouses, family="binomial")