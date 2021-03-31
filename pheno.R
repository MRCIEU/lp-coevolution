library('data.table')
library('dplyr')
source("funs.R")
set.seed(1234)

# load phenotypes
f <- "/tmp/tmp.krHob0dvnc/data.29244.tab"
pheno <- fread(f, select=c(
        "f.eid",
        "f.2734.0.0", 
        "f.2405.0.0",
        "f.21000.0.0",
        "f.40007.0.0",
        "f.34.0.0",
        "f.21022.0.0",
        "f.31.0.0",
        "f.30890.0.0",
        "f.30680.0.0",
        "f.30770.0.0",
        "f.50.0.0",
        "f.51.0.0",
        "f.1687.0.0",
        "f.1697.0.0",
        "f.2714.0.0",
        "f.2375.0.0",
        "f.2385.0.0",
        "f.78.0.0",
        "f.1418.0.0",
        "f.21001.0.0",
        "f.3526.0.0",
        "f.1807.0.0",
        "f.20086.0.0",
        "f.20086.0.1",
        "f.20086.0.2",
        "f.20086.0.3",
        "f.20086.0.4",
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
        "type_of_special_diet_followed_yesterday.20086.0.0",
        "type_of_special_diet_followed_yesterday.20086.0.1",
        "type_of_special_diet_followed_yesterday.20086.0.2",
        "type_of_special_diet_followed_yesterday.20086.0.3",
        "type_of_special_diet_followed_yesterday.20086.0.4"
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