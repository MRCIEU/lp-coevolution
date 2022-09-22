load("./data/pheno.RData")

library("broom")
library("ggplot2")
source("funs.R")
set.seed(1234)

# unrelated

# load covariates & snp
cov <- get_covariates()
pc <- get_genetic_principal_components()

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE)

# unrelated
unrelated <- merge(linker, pheno, by.x="app16729", by.y="eid")
unrelated <- merge(unrelated, cov, "appieu")
unrelated <- merge(unrelated, pc, "appieu")

# load BMI IVs
giant <- fread("data/giant.csv")
giant <- giant %>% dplyr::arrange(p) %>% head(n=20)

# load BMI SNPs
for (i in 1:nrow(giant)){
    snp <- extract_variant_from_bgen(as.character(giant$chr[i]), as.double(giant$position[i]), giant$nea[i], giant$ea[i])
    unrelated <- merge(unrelated, snp, "appieu")
}

# load LP SNP
snp <- extract_variant_from_bgen("2", 136608646, "G", "A")
unrelated <- merge(unrelated, snp, "appieu")

# SD norm BMI
unrelated$body_mass_index.21001.0.0 <- unrelated$body_mass_index.21001.0.0 / sd(unrelated$body_mass_index.21001.0.0, na.rm=T)

# test for asscoiation with LP
results <- data.frame()
snp <- "chr2_136608646_G_A"
f <- paste0("body_mass_index.21001.0.0 ~ ", snp, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC",seq(1, 40), collapse="+"))
fit_all <- lm(as.formula(f), data=unrelated) %>% tidy %>% dplyr::filter(term == snp) %>% dplyr::mutate(Model="All")
fit_milk <- lm(as.formula(f), data=unrelated %>% dplyr::filter(milk_type_used.1418.0.0 == 1)) %>% tidy %>% dplyr::filter(term == snp) %>% dplyr::mutate(Model="Mainly use cows' milk")
fit_nomilk <- lm(as.formula(f), data=unrelated %>% dplyr::filter(milk_type_used.1418.0.0 == 0)) %>% tidy %>% dplyr::filter(term == snp) %>% dplyr::mutate(Model="Soya or never/rarely have cows' milk")
results <- rbind(results, fit_all)
results <- rbind(results, fit_milk)
results <- rbind(results, fit_nomilk)

# test association of SNPs with BMI
for (i in 1:nrow(giant)){
    snp <- paste0("chr", giant$chr[i], "_", giant$position[i], "_", giant$nea[i], "_", giant$ea[i])
    f <- paste0("body_mass_index.21001.0.0 ~ ", snp, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC",seq(1, 40), collapse="+"))
    fit_all <- lm(as.formula(f), data=unrelated) %>% tidy %>% dplyr::filter(term == snp) %>% dplyr::mutate(Model="All")
    fit_milk <- lm(as.formula(f), data=unrelated %>% dplyr::filter(milk_type_used.1418.0.0 == 1)) %>% tidy %>% dplyr::filter(term == snp) %>% dplyr::mutate(Model="Mainly use cows' milk")
    fit_nomilk <- lm(as.formula(f), data=unrelated %>% dplyr::filter(milk_type_used.1418.0.0 == 0)) %>% tidy %>% dplyr::filter(term == snp) %>% dplyr::mutate(Model="Soya or never/rarely have cows' milk")
    results <- rbind(results, fit_all)
    results <- rbind(results, fit_milk)
    results <- rbind(results, fit_nomilk)
}

# rename LP and FTO
results <- results %>% dplyr::mutate(term=dplyr::recode(term, "chr16_53803574_T_A"="chr16_53803574_T_A (FTO)"))
results <- results %>% dplyr::mutate(term=dplyr::recode(term, "chr2_136608646_G_A"="chr2_136608646_G_A (LP)"))

# plot
results$SNP <- factor(results$SNP, levels=rev(unique(results$SNP)))
results$conf.low <- results$estimate - (1.96 * results$std.error)
results$conf.high <- results$estimate + (1.96 * results$std.error)
names(results)[1] <- "SNP"
names(results)[2] <- "Effect"
pdf("bmi.plot.pdf")
ggplot(results, aes(x=SNP, y=Effect, ymin=conf.low, ymax=conf.high, group=Model, color=Model)) +
    geom_point(position=position_dodge(width=1)) +
    geom_errorbar(position=position_dodge(width=1)) +
    coord_flip() +
    theme_classic() +
    labs(y="Effect of SNP on BMI (SD)") +
    theme(legend.position="bottom") +
    geom_hline(data=NULL, mapping=NULL, yintercept=0, linetype="dashed", colour="grey")
dev.off()