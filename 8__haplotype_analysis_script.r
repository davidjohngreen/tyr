library(dplyr)
library(tidyverse)
library(logistf)
library(janitor)

haplotypes = read.csv(file = 'haplotypes.csv', sep='\t',stringsAsFactors = T, header = T)

sample_stats = read.csv(file = 'sample_stats.txt', header = T, stringsAsFactors = T, sep = '\t')
sample_stats = sample_stats %>%
    select(Participant.Id, Platekey, Participant.Phenotypic.Sex,
            Pred.African.Ancestries, Pred.South.Asian.Ancestries,
            Pred.East.Asian.Ancestries, Pred.European.Ancestries,
            Pred.American.Ancestries) %>%
    mutate(ancestry = case_when(Pred.European.Ancestries >= 0.8 ~ "eur",
            Pred.African.Ancestries >= 0.8 ~ "afr",
            Pred.South.Asian.Ancestries >= 0.8 ~ "sas",
            Pred.East.Asian.Ancestries >= 0.8 ~ "eas",
            Pred.American.Ancestries >= 0.8 ~ "amr")) %>%
    select(Participant.Id, Platekey, Participant.Phenotypic.Sex, ancestry) %>%
    rename(LPid = Platekey, participant_id=Participant.Id,sex=Participant.Phenotypic.Sex)


combined = inner_join(haplotypes, sample_stats, by='LPid') %>%
          rename(plate_key = LPid)
additional = read.csv('additional_alleles.csv', header=T, stringsAsFactors=T, sep='\t')
combined = left_join(combined, additional, by='plate_key') %>%
          replace_na(list(additional = 0))

albinism = read.csv(file = 'albinism_IDs.csv', stringsAsFactors = T, header = T, sep = '\t')
albinism = albinism %>% rename(participant_id = partid)
combined = left_join(combined, albinism, by="participant_id")
combined <- mutate(combined, albinism = ifelse(is.na(albinism), 0, albinism))
combined <- mutate(combined, cohort = 'GeL')

common_count <- read.csv(file = 'recoded-status-by-ID.txt', stringsAsFactors = T, header = T, sep = '\t')
common_count <- common_count %>% select('plate_key', 'common_TYR')
testing <- left_join(combined, common_count, by='plate_key')

french_cohort = read.csv(file = 'french_cohort_haplotypes.csv', header=T, stringsAsFactors = T, sep='\t')

full_dataset = bind_rows(testing, french_cohort)
full_dataset$additional <- as.factor(full_dataset$additional)

full_dataset <- full_dataset %>% mutate(CAA = case_when(haplotype == 'CAA' ~ 1)) %>% replace_na(list(CAA = 0))
full_dataset <- full_dataset %>% mutate(CAG = case_when(haplotype == 'CAG' ~ 1)) %>% replace_na(list(CAG = 0))
full_dataset <- full_dataset %>% mutate(CCA = case_when(haplotype == 'CCA' ~ 1)) %>% replace_na(list(CCA = 0))
full_dataset <- full_dataset %>% mutate(CCG = case_when(haplotype == 'CCG' ~ 1)) %>% replace_na(list(CCG = 0))
full_dataset <- full_dataset %>% mutate(TCA = case_when(haplotype == 'TCA' ~ 1)) %>% replace_na(list(TCA = 0))
full_dataset <- full_dataset %>% mutate(TCG = case_when(haplotype == 'TCG' ~ 1)) %>% replace_na(list(TCG = 0))

"""subset for different analyses"""
eur_dataset <- full_dataset %>% filter(ancestry == 'eur')
eur_gel_dataset <- full_dataset %>% filter(ancestry == 'eur', cohort == 'GeL')
gel_dataset <- full_dataset %>% filter(cohort == 'GeL')

"""full_dataset ORs"""
fit1 <- logistf(data=full_dataset, albinism~CAA+additional+ancestry+sex, firth=TRUE, pl=TRUE)
ORs_full_dataset <- data.frame(exp(cbind(OR = coef(fit1),confint(fit1)))) %>% mutate_if(is.numeric, round, digits = 2)

"""eur_dataset ORs"""
fit2 <- logistf(data=eur_dataset, albinism~CAA+additional+sex, firth=TRUE, pl=TRUE)
ORs_eur_dataset <- data.frame(exp(cbind(OR = coef(fit2),confint(fit2)))) %>% mutate_if(is.numeric, round, digits = 2)

"""eur_gel_dataset ORs"""
fit3 <- logistf(data=eur_gel_dataset, albinism~CAA+additional+sex, firth=TRUE, pl=TRUE)
ORs_eur_gel_dataset <- data.frame(exp(cbind(OR = coef(fit3),confint(fit3)))) %>% mutate_if(is.numeric, round, digits = 2)

"""gel_dataset ORs"""
fit4 <- logistf(data=gel_dataset, albinism~CAA+additional+ancestry+sex, firth=TRUE, pl=TRUE)
ORs_gel_dataset <- data.frame(exp(cbind(OR = coef(fit4),confint(fit4)))) %>% mutate_if(is.numeric, round, digits = 2)
