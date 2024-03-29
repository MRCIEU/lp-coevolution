library('data.table')
library('dplyr')
source("funs.R")
set.seed(1234)

sdiet_id<-c()
sdiet_name<-c()
for (i in 0:4){
    sdiet_id <- c(sdiet_id, paste0("x20086_", i, "_", seq(0, 5)))
    sdiet_name <- c(sdiet_name, paste0("type_of_special_diet_followed_yesterday.20086.", i, ".", seq(0, 5)))
}
tdiet_id<-c()
tdiet_name<-c()
for (i in 0:4){
    tdiet_id <- c(tdiet_id, paste0("x100020_", i, "_", 0))
    tdiet_name <- c(tdiet_name, paste0("typical_diet_yesterday.100020.", i, ".0"))
}

# load phenotypes
f <- "/tmp/tmp.5c45ZEm0Nt/data.43017.phesant.tab"
pheno <- fread(f, select=c(
        "eid",
        "x2734_0_0", 
        "x2405_0_0",
        "x21000_0_0",
        "x40007_0_0",
        "x34_0_0",
        "x21022_0_0",
        "x31_0_0",
        "x30890_0_0",
        "x30680_0_0",
        "x30770_0_0",
        "x50_0_0",
        "x51_0_0",
        "x1687_0_0",
        "x1697_0_0",
        "x2714_0_0",
        "x2375_0_0",
        "x2385_0_0",
        "x78_0_0",
        "x1418_0_0",
        "x21001_0_0",
        "x3526_0_0",
        "x1807_0_0",
        "x100530_0_0",
        "x100530_1_0",
        "x100530_2_0",
        "x100530_3_0",
        "x100530_4_0",
        "x100520_0_0",
        "x100520_1_0",
        "x100520_2_0",
        "x100520_3_0",
        "x100520_4_0",
        sdiet_id,
        "x20722_0_0",
        "x20726_0_0",
        "x20747_0_0",
        "x1677_0_0",
        "x20151_0_0",
        "x20150_0_0",
        "x20116_0_0",
        tdiet_id
    ),
    col.names=c(
        "eid",
        "number_of_live_births.2734.0.0",
        "number_of_children_fathered.2405.0.0",
        "ethnic_background.21000.0.0",
        "age_at_death.40007.0.0",
        "year_of_birth.34.0.0",
        "age_at_recruitment.21022.0.0",
        "sex.31.0.0",
        "vitamin_d.30890.0.0",
        "calcium.30680.0.0",
        "igf.30770.0.0",
        "standing_height.50.0.0",
        "seated_height.51.0.0",
        "comparative_body_size_at_age_10.1687.0.0",
        "comparative_height_size_at_age_10.1697.0.0",
        "age_when_periods_started_menarche.2714.0.0",
        "relative_age_of_first_facial_hair.2375.0.0",
        "relative_age_voice_broke.2385.0.0",
        "heel_bone_mineral_density.78.0.0",
        "milk_type_used.1418.0.0",
        "body_mass_index.21001.0.0",
        "mothers_age_at_death.3526.0.0",
        "fathers_age_at_death.1807.0.0",
        "flavoured_milk_intake_yesterday.100530.0.0",
        "flavoured_milk_intake_yesterday.100530.1.0",
        "flavoured_milk_intake_yesterday.100530.2.0",
        "flavoured_milk_intake_yesterday.100530.3.0",
        "flavoured_milk_intake_yesterday.100530.4.0",
        "milk_intake_yesterday.100520.0.0",
        "milk_intake_yesterday.100520.1.0",
        "milk_intake_yesterday.100520.2.0",
        "milk_intake_yesterday.100520.3.0",
        "milk_intake_yesterday.100520.4.0",
        sdiet_name,
        "liking_for_skimmed_milk.20722.0.0",
        "liking_for_soya_milk.20726.0.0",
        "liking_for_whole_milk.20747.0.0",
        "breastfed_as_a_baby.1677.0.0",
        "forced_vital_capacity.20151.0.0",
        "forced_expiratory_volume.20150.0.0",
        "smoking_status.20116.0.0",
        tdiet_name
    )
)
unlink(f)

# process phenotypes
pheno <- pheno %>% mutate_at(c('mothers_age_at_death.3526.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('mothers_age_at_death.3526.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('fathers_age_at_death.1807.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('fathers_age_at_death.1807.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('number_of_live_births.2734.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('number_of_live_births.2734.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('number_of_children_fathered.2405.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('number_of_children_fathered.2405.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('ethnic_background.21000.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('ethnic_background.21000.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('comparative_body_size_at_age_10.1687.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('comparative_body_size_at_age_10.1687.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('comparative_height_size_at_age_10.1697.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('comparative_height_size_at_age_10.1697.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('age_when_periods_started_menarche.2714.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('age_when_periods_started_menarche.2714.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('relative_age_of_first_facial_hair.2375.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('relative_age_of_first_facial_hair.2375.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('relative_age_voice_broke.2385.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('relative_age_voice_broke.2385.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('milk_type_used.1418.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('milk_type_used.1418.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('milk_type_used.1418.0.0'), na_if, 5)
pheno <- pheno %>% mutate_at(c('liking_for_skimmed_milk.20722.0.0'), na_if, -121)
pheno <- pheno %>% mutate_at(c('liking_for_skimmed_milk.20722.0.0'), na_if, -818)
pheno <- pheno %>% mutate_at(c('liking_for_soya_milk.20726.0.0'), na_if, -121)
pheno <- pheno %>% mutate_at(c('liking_for_soya_milk.20726.0.0'), na_if, -818)
pheno <- pheno %>% mutate_at(c('liking_for_whole_milk.20747.0.0'), na_if, -121)
pheno <- pheno %>% mutate_at(c('liking_for_whole_milk.20747.0.0'), na_if, -818)
pheno <- pheno %>% mutate_at(c('breastfed_as_a_baby.1677.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('breastfed_as_a_baby.1677.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('smoking_status.20116.0.0'), na_if, -3)

# recode
pheno <- pheno %>% mutate(comparative_body_size_at_age_10.1687.0.0 = replace(comparative_body_size_at_age_10.1687.0.0, comparative_body_size_at_age_10.1687.0.0 == 1, -1))
pheno <- pheno %>% mutate(comparative_body_size_at_age_10.1687.0.0 = replace(comparative_body_size_at_age_10.1687.0.0, comparative_body_size_at_age_10.1687.0.0 == 3, 0))
pheno <- pheno %>% mutate(comparative_body_size_at_age_10.1687.0.0 = replace(comparative_body_size_at_age_10.1687.0.0, comparative_body_size_at_age_10.1687.0.0 == 2, 1))

pheno <- pheno %>% mutate(comparative_height_size_at_age_10.1697.0.0 = replace(comparative_height_size_at_age_10.1697.0.0, comparative_height_size_at_age_10.1697.0.0 == 1, -1))
pheno <- pheno %>% mutate(comparative_height_size_at_age_10.1697.0.0 = replace(comparative_height_size_at_age_10.1697.0.0, comparative_height_size_at_age_10.1697.0.0 == 3, 0))
pheno <- pheno %>% mutate(comparative_height_size_at_age_10.1697.0.0 = replace(comparative_height_size_at_age_10.1697.0.0, comparative_height_size_at_age_10.1697.0.0 == 2, 1))

pheno <- pheno %>% mutate(relative_age_of_first_facial_hair.2375.0.0 = replace(relative_age_of_first_facial_hair.2375.0.0, relative_age_of_first_facial_hair.2375.0.0 == 2, 0))
pheno <- pheno %>% mutate(relative_age_of_first_facial_hair.2375.0.0 = replace(relative_age_of_first_facial_hair.2375.0.0, relative_age_of_first_facial_hair.2375.0.0 == 1, -1))
pheno <- pheno %>% mutate(relative_age_of_first_facial_hair.2375.0.0 = replace(relative_age_of_first_facial_hair.2375.0.0, relative_age_of_first_facial_hair.2375.0.0 == 3, 1))

pheno <- pheno %>% mutate(relative_age_voice_broke.2385.0.0 = replace(relative_age_voice_broke.2385.0.0, relative_age_voice_broke.2385.0.0 == 2, 0))
pheno <- pheno %>% mutate(relative_age_voice_broke.2385.0.0 = replace(relative_age_voice_broke.2385.0.0, relative_age_voice_broke.2385.0.0 == 1, -1))
pheno <- pheno %>% mutate(relative_age_voice_broke.2385.0.0 = replace(relative_age_voice_broke.2385.0.0, relative_age_voice_broke.2385.0.0 == 3, 1))

pheno <- pheno %>% mutate(milk_type_used.1418.0.0 = replace(milk_type_used.1418.0.0, milk_type_used.1418.0.0 == 2, 1))
pheno <- pheno %>% mutate(milk_type_used.1418.0.0 = replace(milk_type_used.1418.0.0, milk_type_used.1418.0.0 == 3, 1))
pheno <- pheno %>% mutate(milk_type_used.1418.0.0 = replace(milk_type_used.1418.0.0, milk_type_used.1418.0.0 == 4, 0))
pheno <- pheno %>% mutate(milk_type_used.1418.0.0 = replace(milk_type_used.1418.0.0, milk_type_used.1418.0.0 == 6, 0))

pheno <- pheno %>% mutate(milk_intake_yesterday.100520.0.0 = replace(milk_intake_yesterday.100520.0.0, milk_intake_yesterday.100520.0.0 == 555, 0.5))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.0.0 = replace(milk_intake_yesterday.100520.0.0, milk_intake_yesterday.100520.0.0 == 600, 6))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.1.0 = replace(milk_intake_yesterday.100520.1.0, milk_intake_yesterday.100520.1.0 == 555, 0.5))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.1.0 = replace(milk_intake_yesterday.100520.1.0, milk_intake_yesterday.100520.1.0 == 600, 6))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.2.0 = replace(milk_intake_yesterday.100520.2.0, milk_intake_yesterday.100520.2.0 == 555, 0.5))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.2.0 = replace(milk_intake_yesterday.100520.2.0, milk_intake_yesterday.100520.2.0 == 600, 6))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.3.0 = replace(milk_intake_yesterday.100520.3.0, milk_intake_yesterday.100520.3.0 == 555, 0.5))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.3.0 = replace(milk_intake_yesterday.100520.3.0, milk_intake_yesterday.100520.3.0 == 600, 6))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.4.0 = replace(milk_intake_yesterday.100520.4.0, milk_intake_yesterday.100520.4.0 == 555, 0.5))
pheno <- pheno %>% mutate(milk_intake_yesterday.100520.4.0 = replace(milk_intake_yesterday.100520.4.0, milk_intake_yesterday.100520.4.0 == 600, 6))

# mean daily milk intake
milk_intake_yesterday.100520 <- apply(pheno[,paste0("milk_intake_yesterday.100520.", seq(0, 4), ".0"), with=F], 1, function(x) {m <- mean(x, na.rm=T); if (is.nan(m)) {NA} else {m}})
milk_intake_yesterday.100520 <- as.data.table(milk_intake_yesterday.100520)
pheno <- cbind(pheno, milk_intake_yesterday.100520)

pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.0.0 = replace(flavoured_milk_intake_yesterday.100530.0.0, flavoured_milk_intake_yesterday.100530.0.0 == 555, 0.5))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.0.0 = replace(flavoured_milk_intake_yesterday.100530.0.0, flavoured_milk_intake_yesterday.100530.0.0 == 600, 6))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.1.0 = replace(flavoured_milk_intake_yesterday.100530.1.0, flavoured_milk_intake_yesterday.100530.1.0 == 555, 0.5))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.1.0 = replace(flavoured_milk_intake_yesterday.100530.1.0, flavoured_milk_intake_yesterday.100530.1.0 == 600, 6))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.2.0 = replace(flavoured_milk_intake_yesterday.100530.2.0, flavoured_milk_intake_yesterday.100530.2.0 == 555, 0.5))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.2.0 = replace(flavoured_milk_intake_yesterday.100530.2.0, flavoured_milk_intake_yesterday.100530.2.0 == 600, 6))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.3.0 = replace(flavoured_milk_intake_yesterday.100530.3.0, flavoured_milk_intake_yesterday.100530.3.0 == 555, 0.5))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.3.0 = replace(flavoured_milk_intake_yesterday.100530.3.0, flavoured_milk_intake_yesterday.100530.3.0 == 600, 6))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.4.0 = replace(flavoured_milk_intake_yesterday.100530.4.0, flavoured_milk_intake_yesterday.100530.4.0 == 555, 0.5))
pheno <- pheno %>% mutate(flavoured_milk_intake_yesterday.100530.4.0 = replace(flavoured_milk_intake_yesterday.100530.4.0, flavoured_milk_intake_yesterday.100530.4.0 == 600, 6))

# mean daily flavoured milk intake
flavoured_milk_intake_yesterday.100530 <- apply(pheno[,paste0("flavoured_milk_intake_yesterday.100530.", seq(0, 4), ".0"), with=F], 1, function(x) {m <- mean(x, na.rm=T); if (is.nan(m)) {NA} else {m}})
flavoured_milk_intake_yesterday.100530 <- as.data.table(flavoured_milk_intake_yesterday.100530)
pheno <- cbind(pheno, flavoured_milk_intake_yesterday.100530)

# derive lactose-free diet variable
lactose_free_diet <- apply(pheno[,sdiet_name,with=F], 1, function(x) {sum(x==9, na.rm=T)>0})
lactose_free_diet <- as.data.table(lactose_free_diet)
pheno <- cbind(pheno, lactose_free_diet)
lactose_free_na <- apply(pheno[,tdiet_name,with=F], 1, function(x) {all(is.na(x))})
pheno[lactose_free_na]$lactose_free_diet <- NA

# derive lactose-free diet variable2
pheno$lactose_free_diet.0 <- apply(pheno[,paste0("type_of_special_diet_followed_yesterday.20086.", 0, ".", seq(0, 5)),with=F], 1, function(x) {sum(x==9, na.rm=T)>0})
pheno$lactose_free_diet.1 <- apply(pheno[,paste0("type_of_special_diet_followed_yesterday.20086.", 1, ".", seq(0, 5)),with=F], 1, function(x) {sum(x==9, na.rm=T)>0})
pheno$lactose_free_diet.2 <- apply(pheno[,paste0("type_of_special_diet_followed_yesterday.20086.", 2, ".", seq(0, 5)),with=F], 1, function(x) {sum(x==9, na.rm=T)>0})
pheno$lactose_free_diet.3 <- apply(pheno[,paste0("type_of_special_diet_followed_yesterday.20086.", 3, ".", seq(0, 5)),with=F], 1, function(x) {sum(x==9, na.rm=T)>0})
pheno$lactose_free_diet.4 <- apply(pheno[,paste0("type_of_special_diet_followed_yesterday.20086.", 4, ".", seq(0, 5)),with=F], 1, function(x) {sum(x==9, na.rm=T)>0})
pheno$lactose_free_na.0 <- is.na(pheno$typical_diet_yesterday.100020.0.0)
pheno$lactose_free_na.1 <- is.na(pheno$typical_diet_yesterday.100020.1.0)
pheno$lactose_free_na.2 <- is.na(pheno$typical_diet_yesterday.100020.2.0)
pheno$lactose_free_na.3 <- is.na(pheno$typical_diet_yesterday.100020.3.0)
pheno$lactose_free_na.4 <- is.na(pheno$typical_diet_yesterday.100020.4.0)
pheno$lactose_free_diet.0[pheno$lactose_free_na.0] <- NA
pheno$lactose_free_diet.1[pheno$lactose_free_na.1] <- NA
pheno$lactose_free_diet.2[pheno$lactose_free_na.2] <- NA
pheno$lactose_free_diet.3[pheno$lactose_free_na.3] <- NA
pheno$lactose_free_diet.4[pheno$lactose_free_na.4] <- NA

# derived phenotypes
pheno$leg_length <- pheno$standing_height.50.0.0 - pheno$seated_height.51.0.0
pheno$leg_to_trunk_length_ratio <- pheno$leg_length / pheno$seated_height.51.0.0
pheno$dead <- as.integer(!is.na(pheno$age_at_death.40007.0.0))
pheno <- pheno %>% mutate(time = case_when(
    is.na(age_at_death.40007.0.0) ~ 2020 - year_of_birth.34.0.0,
    !is.na(age_at_death.40007.0.0)  ~ age_at_death.40007.0.0
    )
)

# save data
save.image(file = "data/pheno.RData")