library('data.table')
library('broom')
library('optparse')
library('jsonlite')
library('httr')
library('TwoSampleMR')
library('rbgen')
library('stringr')
set.seed(1244)

# lactase persistence 
snp <- 'rs4988235'

# load spouse list and restrict to couples
spouses <- fread("./Spouse-pairs-IEU.txt")
names(spouses) <- c("appieu", "couple")

# obtain SNPs from UKBB
for (chr in c(2)){
    message(paste("Obtaining SNPs from chr:", chr))

    if (chr < 10){
        path <- paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr0", chr, ".bgen")
    } else {
        path <- paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr", chr, ".bgen")
    }

    # get variant
    variant <- bgen.load(
        path,
        rsids = snp
    )

    # convert dosage for each genotype to copies of alt allele
    dosage <- as.data.frame(
        apply(variant$data, 1, function(data) { return(data[,1]*0 + data[,2]*1 + data[,3]*2) })
    )
    names(dosage) <- paste(names(dosage), variant$variants$allele0, variant$variants$allele1, sep="_")
    dosage$appieu <- row.names(dosage)

    # store
    spouses <- merge(spouses, dosage, "appieu")
}

# load phenotypes
load("../../pheno.RData")

# link & merge genotypes & phenotypes
linker <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/16729/2019-04-29/data/ieu_id_linker.csv")
spouses <- merge(spouses, linker, "appieu")
spouses <- merge(spouses, pheno[c("f.eid", "Milk.lactose", "Sex")], by.y="f.eid", by.x="app16729")
names(spouses)[names(spouses) == "Milk.lactose"] <- "index.Milk.lactose"
names(spouses)[names(spouses) == "Sex"] <- "index.Sex"
spouses <- na.omit(spouses)

# drop participants where one partner has missings data
counts <- as.data.frame(table(spouses$couple))
counts <- counts[which(counts$Freq == 2),]
spouses <- spouses[spouses$couple %in% counts$Var1]

# merge partner variables
spouses$couple.Sex <- paste(spouses$couple, spouses$index.Sex, sep=".")
partner <- spouses[,c("index.Sex", "index.Milk.lactose", "couple", "rs4988235_G_A")]
partner$couple.Sex <- paste(partner$couple, as.numeric(partner$index.Sex != 1), sep=".")
names(partner) <- c("partner.Sex", "partner.Milk.lactose", "couple", "partner.rs4988235_G_A", "couple.Sex")
partner$couple <- NULL
spouses <- merge(spouses, partner, "couple.Sex")

# Main

# snp-exposure estimate (snp-index-Milk.lactose assoc)
index.snp.assoc <- data.frame()
for (snp in names(spouses)[grep("^rs", names(spouses))]){
    f <- as.formula(paste0("spouses$index.Milk.lactose ~ spouses$index.Sex + spouses[[snp]] + spouses$Age.at.recruitment.x + ", paste0("spouses$PC",seq(1,10), ".x", collapse=" + ")))
    fit <- tidy(glm(f, family="binomial"))
    j <- str_split(snp, "_", simplify = TRUE)
    index.snp.assoc <- rbind(data.frame(SNP=j[1], other_allele=j[2], effect_allele=j[3], beta=fit$estimate[3], se=fit$std.error[3], pval=fit$p.value[3]), index.snp.assoc, stringsAsFactors=F)
}

# snp-outcome estimate (snp-partner-Milk.lactose assoc)
partner.snp.assoc <- data.frame()
for (snp in names(spouses)[grep("^rs", names(spouses))]){
    f <- as.formula(paste0("spouses$partner.Milk.lactose ~ spouses$partner.Sex + spouses[[snp]] + spouses$Age.at.recruitment.x + spouses$Age.at.recruitment.y + ", paste0("spouses$PC",seq(1,10), ".x", collapse=" + "), " + ", paste0("spouses$PC",seq(1,10), ".y", collapse=" + ")))
    fit <- tidy(glm(f, family="binomial"))
    j <- str_split(snp, "_", simplify = TRUE)
    partner.snp.assoc <- rbind(data.frame(SNP=j[1], other_allele=j[2], effect_allele=j[3], beta=fit$estimate[3], se=fit$std.error[3], pval=fit$p.value[3]), partner.snp.assoc, stringsAsFactors=F)
}

# TwoSampleMR
exp_dat <- format_data(index.snp.assoc, type="exposure")
out_dat <- format_data(partner.snp.assoc, type="outcome")

# harmonise
dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat
)

# MR
res <- mr(dat)
res$lci <- res$b - (res$se * 1.96)
res$uci <- res$b + (res$se * 1.96)
res$b <- exp(res$b)
res$lci <- exp(res$lci)
res$uci <- exp(res$uci)

# Sensitivity

# snp-exposure estimate (snp-index-Milk.lactose assoc)
index.snp.assoc <- data.frame()
for (snp in names(spouses)[grep("^rs", names(spouses))]){
    f <- as.formula(paste0("spouses$index.Milk.lactose ~ spouses$index.Sex + spouses[[snp]] + spouses$Age.at.recruitment.x + ", paste0("spouses$PC",seq(1,40), ".x", collapse=" + ")))
    fit <- tidy(glm(f, family="binomial"))
    j <- str_split(snp, "_", simplify = TRUE)
    index.snp.assoc <- rbind(data.frame(SNP=j[1], other_allele=j[2], effect_allele=j[3], beta=fit$estimate[3], se=fit$std.error[3], pval=fit$p.value[3]), index.snp.assoc, stringsAsFactors=F)
}

# snp-outcome estimate (snp-partner-Milk.lactose assoc)
partner.snp.assoc <- data.frame()
for (snp in names(spouses)[grep("^rs", names(spouses))]){
    f <- as.formula(paste0("spouses$partner.Milk.lactose ~ spouses$partner.Sex + spouses[[snp]] + spouses$Age.at.recruitment.x + spouses$Age.at.recruitment.y + ", paste0("spouses$PC",seq(1,40), ".x", collapse=" + "), " + ", paste0("spouses$PC",seq(1,40), ".y", collapse=" + ")))
    fit <- tidy(glm(f, family="binomial"))
    j <- str_split(snp, "_", simplify = TRUE)
    partner.snp.assoc <- rbind(data.frame(SNP=j[1], other_allele=j[2], effect_allele=j[3], beta=fit$estimate[3], se=fit$std.error[3], pval=fit$p.value[3]), partner.snp.assoc, stringsAsFactors=F)
}

# TwoSampleMR
exp_dat <- format_data(index.snp.assoc, type="exposure")
out_dat <- format_data(partner.snp.assoc, type="outcome")

# harmonise
dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat
)

# MR
res <- mr(dat)
res$lci <- res$b - (res$se * 1.96)
res$uci <- res$b + (res$se * 1.96)
res$b <- exp(res$b)
res$lci <- exp(res$lci)
res$uci <- exp(res$uci)


# genetic correlation

## additive; main
f <- as.formula(paste0("spousesgc$rs4988235_G_A.y ~ spousesgc$rs4988235_G_A.x + spousesgc$Age.at.recruitment.x + spousesgc$Age.at.recruitment.y + ", paste0("spousesgc$PC",seq(1,10), ".x", collapse=" + "), " + ", paste0("spousesgc$PC",seq(1,10), ".y", collapse=" + ")))
fit <- lm(f)

## additive; sens
f <- as.formula(paste0("spousesgc$rs4988235_G_A.y ~ spousesgc$rs4988235_G_A.x + spousesgc$Age.at.recruitment.x + spousesgc$Age.at.recruitment.y + ", paste0("spousesgc$PC",seq(1,40), ".x", collapse=" + "), " + ", paste0("spousesgc$PC",seq(1,40), ".y", collapse=" + ")))
fit <- lm(f)

## dominant; main
f <- as.formula(paste0("spousesgc$rs4988235_G_A_dom_rd.y ~ spousesgc$rs4988235_G_A_dom_rd.x + spousesgc$Age.at.recruitment.x + spousesgc$Age.at.recruitment.y + ", paste0("spousesgc$PC",seq(1,10), ".x", collapse=" + "), " + ", paste0("spousesgc$PC",seq(1,10), ".y", collapse=" + ")))
fit <- glm(f, family="binomial")

## dominant; sens
f <- as.formula(paste0("spousesgc$rs4988235_G_A_dom_rd.y ~ spousesgc$rs4988235_G_A_dom_rd.x + spousesgc$Age.at.recruitment.x + spousesgc$Age.at.recruitment.y + ", paste0("spousesgc$PC",seq(1,40), ".x", collapse=" + "), " + ", paste0("spousesgc$PC",seq(1,40), ".y", collapse=" + ")))
fit <- glm(f, family="binomial")