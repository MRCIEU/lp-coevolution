library('data.table')
library('dplyr')
library('survival')
library('sandwich')
library('lmtest')
library('plm')
library("broom")
library('optparse')
library('forestplot')
library('data.table')
library('stringr')
library("gridExtra")
source("funs.R")
load("./data/pheno.RData")
setEPS()
set.seed(1234)

# Adapted from: https://github.com/nmdavies/within_family_mr/blob/master/HUNT/reg_2_option_c_individual_snps-v5.R
related_plm <- function(df, out, sd=T){
    if (sd){
        df[[out]] <- df[[out]] / sd(df[[out]], na.rm=T)
    }
    # extract variables
    t <- df[,c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "famid", "chr2_136608646_G_A", out),with=F]
    t$famid <- as.factor(t$famid)

    # drop missing values
    t <- na.omit(t)

    # drop singletons
    t <- t %>%
        group_by(famid) %>%
        filter(n() > 1)
    
    #Run plm model with family fixed effects + robust standard errors
    f <- as.formula(paste(out, "~ chr2_136608646_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    estimates<-plm(f, data = t, index = "famid", model = "within", inst.method = "bvk")
    estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
    
    beta<-estimates_robust[1,1]
    se<-estimates_robust[1,2]
    pvalue<-estimates_robust[1,4]
    sample_size<-nrow(model.frame(estimates))
    EAF<-mean(model.frame(estimates)[,2])/2

    return(data.frame(bin=F, beta, se, pvalue, sample_size, out, model="within-family"))
}

unrelated_lm <- function(df, out, bin=F, sd=T){
    if (sd){
        df[[out]] <- df[[out]] / sd(df[[out]], na.rm=T)
    }
    if (bin){
        fit <- glm(as.formula(paste(out, "~ chr2_136608646_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=df, family="binomial")
    } else {
        fit <- lm(as.formula(paste(out, "~ chr2_136608646_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=df)
    }
    n=length(resid(fit))
    fit <- tidy(fit)
    return(data.frame(bin,model="unrelated", sample_size=n, beta=fit$estimate[2], pvalue=fit$p.value[2], se=fit$std.error[2], out))
}

# Adapted from: https://github.com/LaurenceHowe/SiblingGWAS/blob/master/resources/regression/unified_regression.R
related_lm <- function(df, out, sd=T, bin=F){
    if (sd){
        df[[out]] <- df[[out]] / sd(df[[out]], na.rm=T)
    }
    # extract variables
    t <- df[,c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "famid", "chr2_136608646_G_A", out),with=F]
    t$famid <- as.factor(t$famid)

    # drop missing values
    t <- na.omit(t)

    # drop singletons
    t <- t %>%
        group_by(famid) %>%
        filter(n() > 1)

    # set params
    t$FAM_MEAN <- ave(as.numeric(unlist(t[,"chr2_136608646_G_A"])), t$famid, FUN=mean)
    t$CENTREDGENOTYPE <- t$chr2_136608646_G_A-t$FAM_MEAN

    # linear model
    f <- as.formula(paste(out, "~ FAM_MEAN + CENTREDGENOTYPE + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    if (bin){
        model2 <- glm(f, data=t, family="binomial")
    } else {
        model2 <- lm(f, data=t)
    }

    # effective sample size
    n <- length(resid(model2))

    # mean effect of X on Y WF
    BETA_WF <- as.numeric(model2$coefficients[3])

    # save the variance covariance matrix
    vcv_matrix <- vcovCL(model2, cluster=t$famid)

    # save the clustered SEs and corresponding P-values for WF/BF
    test_matrix <- coeftest(model2, vcov.=vcv_matrix)

    # SE taking account for family design
    SE_BETA_WF <- test_matrix[3,2]
    P_BETA_WF <- test_matrix[3,4]

    return(data.frame(bin, model="within-family", sample_size=n, beta=BETA_WF, pvalue=P_BETA_WF, se=SE_BETA_WF, out))
}

# load covariates & snp
cov <- get_covariates()
pc <- get_genetic_principal_components()
snp <- extract_variant_from_bgen("2", 136608646, "G", "A")

# unrelated

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE)

# unrelated
unrelated <- merge(linker, pheno, by.x="app16729", by.y="eid")
unrelated <- merge(unrelated, cov, "appieu")
unrelated <- merge(unrelated, pc, "appieu")
unrelated <- merge(unrelated, snp, "appieu")

# sibs
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=FALSE, drop_related=FALSE)
sibs <- get_siblings()
fam <- merge(linker, sibs, "appieu")
related <- merge(fam, pheno, by.x="app16729", by.y="eid")
related <- merge(related, cov, "appieu")
related <- merge(related, pc, "appieu")
related <- merge(related, snp, "appieu")

# unrelated SNP-estimates
main <- data.frame()
main <- rbind(main, unrelated_lm(unrelated, "vitamin_d.30890.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "calcium.30680.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "igf.30770.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "standing_height.50.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "heel_bone_mineral_density.78.0.0", sd=F))
main <- rbind(main, unrelated_lm(unrelated, "body_mass_index.21001.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "mothers_age_at_death.3526.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "fathers_age_at_death.1807.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "number_of_live_births.2734.0.0"))
main <- rbind(main, unrelated_lm(unrelated, "number_of_children_fathered.2405.0.0"))
main$lci <- main$beta - (main$se * 1.96)
main$uci <- main$beta + (main$se * 1.96)

# within-fam SNP-estimates
wf <- data.frame()
wf <- rbind(wf, related_plm(related, "vitamin_d.30890.0.0"))
wf <- rbind(wf, related_plm(related, "calcium.30680.0.0"))
wf <- rbind(wf, related_plm(related, "igf.30770.0.0"))
wf <- rbind(wf, related_plm(related, "standing_height.50.0.0"))
wf <- rbind(wf, related_plm(related, "heel_bone_mineral_density.78.0.0", sd=F))
wf <- rbind(wf, related_plm(related, "body_mass_index.21001.0.0"))
wf <- rbind(wf, related_plm(related, "mothers_age_at_death.3526.0.0"))
wf <- rbind(wf, related_plm(related, "fathers_age_at_death.1807.0.0"))
wf <- rbind(wf, related_plm(related, "number_of_live_births.2734.0.0"))
wf <- rbind(wf, related_plm(related, "number_of_children_fathered.2405.0.0"))
wf$lci <- wf$beta - (wf$se * 1.96)
wf$uci <- wf$beta + (wf$se * 1.96)

results1 <- merge(main, wf, "out", suffixes = c(".main",".wf"))
results1$out <- as.character(results1$out)
results1$out <- str_split(results1$out, "\\.", simplify = TRUE)[,1]
results1$out <- gsub("_", " ", results1$out)
results1$out <- stringr::str_to_title(results1$out)
results1$out <- str_replace(results1$out, "Igf", "IGF-1")
results1$out <- str_replace(results1$out, "Mothers Age At Death", "Mother's Age At Death")
results1$out <- str_replace(results1$out, "Fathers Age At Death", "Father's Age At Death")

# OR/HR

# Cox Proportional-Hazards Model
fit <- coxph(Surv(time, dead) ~ chr2_136608646_G_A + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = unrelated)

main <- data.frame()
main <- rbind(main, unrelated_lm(unrelated, "milk_type_used.1418.0.0", bin=T, sd=F))
main <- rbind(main, data.frame(bin=T, out="mortality", model="unrelated", sample_size=length(resid(fit)), beta=tidy(fit)$estimate[1], pvalue=tidy(fit)$p.value[1], se=tidy(fit)$std.error[1]))
main$lci <- exp(main$beta - (main$se * 1.96))
main$uci <- exp(main$beta + (main$se * 1.96))
main$beta <- exp(main$beta)

wf <- data.frame()
wf <- rbind(wf, related_lm(related, "milk_type_used.1418.0.0", bin=T, sd=F))
wf <- rbind(wf, data.frame(bin=T, out="mortality", model="within-family", sample_size=NA, beta=NA, pvalue=NA, se=NA))
wf$lci <- exp(wf$beta - (wf$se * 1.96))
wf$uci <- exp(wf$beta + (wf$se * 1.96))
wf$beta <- exp(wf$beta)

results2 <- merge(main, wf, "out", suffixes = c(".main",".wf"))
results2$out <- as.character(results2$out)
results2$out <- str_split(results2$out, "\\.", simplify = TRUE)[,1]
results2$out <- gsub("_", " ", results2$out)
results2$out <- stringr::str_to_title(results2$out)
results2$out <- str_replace(results2$out, "Milk Type Used", "Cows' Milk Consumer (OR)")
results2$out <- str_replace(results2$out, "Mortality", "Mortality (HR)")

# Print two plots side by side using the grid
# package's layout option for viewports
postscript("forest.eps", width=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(
    results1$out,
    boxsize = 0.05,
    xticks = c(-0.25, -0.1, 0, 0.1, 0.25),
    legend = c("Unrelated", "Within-family"),
    mean = results1[,c("beta.main", "beta.wf")],
    lower = results1[,c("lci.main", "lci.wf")],
    upper = results1[,c("uci.main", "uci.wf")],
    col=fpColors(box=c("darkblue", "darkred"), lines=c("darkblue", "darkred")),
    xlab="Lactase persistence allele effect on outcome (SD)",
    new_page=F,
    title="A"
)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
forestplot(
    results2$out,
    boxsize = 0.05,
    xticks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75),
    legend = c("Unrelated", "Within-family"),
    mean = results2[,c("beta.main", "beta.wf")],
    lower = results2[,c("lci.main", "lci.wf")],
    upper = results2[,c("uci.main", "uci.wf")],
    col=fpColors(box=c("darkblue", "darkred"), lines=c("darkblue", "darkred")),
    zero = 1,
    xlab="Lactase persistence allele effect on outcome (OR/HR)",
    new_page=F,
    title="B"
)
popViewport(2)
dev.off()