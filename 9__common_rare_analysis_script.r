library(dplyr)
library(tidyverse)
library(logistf)
library(janitor)

haplotypes = read.csv(file = 'haplotype_groupings_full_genotypes.csv', sep='\t',stringsAsFactors = T, header = T)

sample_stats = read.csv(file = 'sample_stats.txt', header = T, stringsAsFactors = T, sep = '\t')
sample_stats = sample_stats %>% select(Participant.Id, Platekey, Participant.Phenotypic.Sex, Pred.African.Ancestries, Pred.South.Asian.Ancestries, Pred.East.Asian.Ancestries, Pred.European.Ancestries, Pred.American.Ancestries)
sample_stats = sample_stats %>% mutate(ancestry = case_when(Pred.European.Ancestries >= 0.8 ~ "eur",
                            Pred.African.Ancestries >= 0.8 ~ "afr",
                            Pred.South.Asian.Ancestries >= 0.8 ~ "sas",
                            Pred.East.Asian.Ancestries >= 0.8 ~ "eas",
                            Pred.American.Ancestries >= 0.8 ~ "amr"))

sample_stats = sample_stats %>% select(Participant.Id, Platekey, Participant.Phenotypic.Sex, ancestry)
sample_stats = sample_stats %>% rename(LPid = Platekey, participant_id=Participant.Id,sex=Participant.Phenotypic.Sex)

combined = inner_join(haplotypes, sample_stats, by='LPid')
combined = combined %>% rename(plate_key = LPid)

additional = read.csv('additional.csv', header=T, stringsAsFactors=T, sep='\t')
combined = left_join(combined, additional, by='plate_key')
combined <- combined %>% replace_na(list(additional = 0))


albinism = read.csv(file = 'albinism_IDs.csv', stringsAsFactors = T, header = T, sep = '\t')
albinism = albinism %>% rename(participant_id = partid)

combined = left_join(combined, albinism, by="participant_id")
combined  <- mutate(combined, albinism = ifelse(is.na(albinism), 0, albinism))
combined <- mutate(combined, cohort = 'GeL')

common_count <- read.csv(file = 'recoded-status-by-ID-full-genotypes.txt', stringsAsFactors = T, header = T, sep = '\t')
common_count <- common_count %>% select('plate_key', 'common_TYR')
testing <- inner_join(combined, common_count, by='plate_key')

french = read.csv(file = 'french_cohort_common_rare.csv', header=T, stringsAsFactors = T, sep='\t')

final = bind_rows(testing, french)

"""get those with zero common and zero rare"""
test <- final %>% mutate(zerozero = case_when(common_TYR == 0 & additional == 0 ~ 1)) %>% replace_na(list(zerozero = 0))


"""run logistf to get the OR for those with zero common and zero rare variants"""
fit1<-logistf(data=test, albinism~zerozero+ancestry+sex, firth=TRUE, pl=TRUE)
ratios <- data.frame(exp(cbind(OR = coef(fit1),confint(fit1))))
ratios_zero <- data.frame(ratios)



"""get those with different combinations of common and rare, setting 0 common, 0 rare as the reference group"""
test = final %>% mutate(common_rare = case_when(
                            common_TYR == 6 & additional == 0 ~ "6_common_0_rare",
                            common_TYR == 5 & additional == 0 ~ "5_common_0_rare",
                            common_TYR == 4 & additional == 0 ~ "4_common_0_rare",
                            common_TYR == 3 & additional == 0 ~ "3_common_0_rare",
                            common_TYR == 2 & additional == 0 ~ "2_common_0_rare",
                            common_TYR == 1 & additional == 0 ~ "1_common_0_rare",
                            common_TYR == 0 & additional == 0 ~ "0_common_0_rare",

                            common_TYR == 6 & additional == 1 ~ "6_common_1_rare",
                            common_TYR == 5 & additional == 1 ~ "5_common_1_rare",
                            common_TYR == 4 & additional == 1 ~ "4_common_1_rare",
                            common_TYR == 3 & additional == 1 ~ "3_common_1_rare",
                            common_TYR == 2 & additional == 1 ~ "2_common_1_rare",
                            common_TYR == 1 & additional == 1 ~ "1_common_1_rare",
                            common_TYR == 0 & additional == 1 ~ "0_common_1_rare",

                            common_TYR == 0 & additional == 2 ~ "0_common_2_rare"))


test$common_rare <- as.factor(test$common_rare)


"""run the regression on the total dataset"""
fit1<-logistf(data=test, albinism~common_rare+ancestry+sex, firth=TRUE, pl=TRUE)
ratios <- data.frame(exp(cbind(OR = coef(fit1),confint(fit1)))) %>% mutate_if(is.numeric, round, digits = 2)
ratios_common_rare <- data.frame(ratios)



"""get those with different combinations of common and rare TYR variants, setting 0 common, 0 rare as the reference group"""
test = test %>% mutate(common_rare_TYR = case_when(
                            common_TYR == 6 & additional_tyr == 1 ~ "6_common_1_rare_tyr",
                            common_TYR == 5 & additional_tyr == 1 ~ "5_common_1_rare_tyr",
                            common_TYR == 4 & additional_tyr == 1 ~ "4_common_1_rare_tyr",
                            common_TYR == 3 & additional_tyr == 1 ~ "3_common_1_rare_tyr",
                            common_TYR == 2 & additional_tyr == 1 ~ "2_common_1_rare_tyr",
                            common_TYR == 1 & additional_tyr == 1 ~ "1_common_1_rare_tyr",
                            common_TYR == 0 & additional_tyr == 1 ~ "0_common_1_rare_tyr",
                            common_TYR == 0 & additional_tyr == 2 ~ "0_common_2_rare_tyr",
                            common_TYR == 0 & additional     == 0 ~ "0_common_0_rare"))

test$common_rare_TYR <- as.factor(test$common_rare_TYR)


fit1<-logistf(data=test, albinism~common_rare_TYR+ancestry+sex, firth=TRUE, pl=TRUE)
ratios <- data.frame(exp(cbind(OR = coef(fit1),confint(fit1))))
tyr <- data.frame(ratios)
