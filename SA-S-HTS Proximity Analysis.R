#South African SARS-CoV-2 Household Transmission Study
    #Data used for Kleynhans et al 2023 
      #Household transmission of SARS-CoV-2 from adult index cases living with 
        #and without HIV in South Africa, 2020-2021: a case-ascertained, 
        #prospective observational household transmission study

#Purpose: Describe close proximity events and assess association between close proximity events and SARS-CoV-2 transmission
#Created by: Jackie Kleynhans, National Institute for Communicable Diseases
#Contact: jackiel@nicd.ac.za

#Updates made March 2023 (post peer-review)
  #Dataset: includes variable (members_excluded) to flag households where some members were excluded from analysis due to baseline sero-positivity
  #Script
      #Additional contact parameters added: maximum contact frequency and duration, average frequency, cumulative time in contact

#Clear workspace
rm(list=ls())

#Load Libraries
library(readr)
library(tidyverse)
library(readxl)
library(dplyr)
library(lme4)
library(broom.mixed)
library(tableone)
library(ggpubr)
library(fields)
library(gtsummary)

#Set working directory
setwd("~/COVID-19 HTS/Proximity/Data/Kleynhans 2023 SA-S-HTS Proximity")

#Import data
hts_cnet_filt_sym <- read_csv("sashts_contact_network.csv")
meta <- read_csv("sashts_metadata.csv") %>% 
  mutate(agegrp9 = if_else(agegrp9=="12-May", "5-12", agegrp9)) %>% 
  mutate(agegrp9 = if_else(agegrp9=="60", "=60", agegrp9)) %>% 
  mutate(ixagegrp9 = if_else(ixagegrp9=="60", "=60", ixagegrp9))

#Set reference categories
meta$site <- factor(meta$site, 
                    levels = c("Klerksdorp", "Soweto"), 
                    labels=c("Klerksdorp", "Soweto"))
meta$site <- relevel(meta$site, 
                     ref="Klerksdorp")
meta$ixagegrp9 <- factor(meta$ixagegrp9, 
                         levels = c("18-34","35-59", "=60"), 
                         labels = c("18-34","35-59", "\u226560"))
meta$ixagegrp9 <- relevel(meta$ixagegrp9, 
                          ref="18-34")
meta$ixminctcat1 <- factor(meta$ixminctcat1, 
                           levels = c("<25", "25-35", ">35", "Unknown"), 
                           labels = c("<25", "25-35", ">35", "Unknown"))
meta$ixminctcat1 <- relevel(meta$ixminctcat1, 
                            ref=">35")
meta$ixminctcat1 <- factor(meta$ixminctcat1, 
                           levels = c("<25", "25-35", ">35", "Unknown"), 
                           labels = c("<25", "25-35", ">35", "Unknown"))
meta$ixminctcat1 <- relevel(meta$ixminctcat1, ref=">35")
meta$ixminctcat2 <- factor(meta$ixminctcat2, 
                           levels = c(">35", "30-35","<30","Unknown"))
meta$ixminctcat2 <- relevel(meta$ixminctcat2, 
                            ref=">35")
meta$ixesarsvarf1 <- factor(meta$ixesarsvarf1, 
                            levels = c("Variant Unknown", "non-Alpha/Beta/Delta", "Alpha", "Beta", "Delta"))
meta$ixesarsvarf1 <- relevel(meta$ixesarsvarf1, ref="Beta")
meta$agegrp9 <- factor(meta$agegrp9, 
                       levels = c("<5", "5-12", "13-17", "18-34","35-59", "=60"), 
                       labels = c("<5", "5-12", "13-17", "18-34","35-59","\u226560"))
meta$agegrp9 <- relevel(meta$agegrp9, ref="<5")
meta$smokecignow1 <- factor(meta$smokecignow1)
meta$smokecignow1 <- relevel(meta$smokecignow1, 
                             ref="Yes")
meta$bmicat <- factor(meta$bmicat, 
                      levels = c("Underweight", "Normal weight", "Overweight", "Obese", "Unknown"))
meta$bmicat <- relevel(meta$bmicat, 
                       ref="Normal weight")
meta$sex <- factor(meta$sex, 
                   levels = c("Male", "Female"))
meta$sex <- relevel(meta$sex, 
                    ref="Male")
meta$ixsex <- factor(meta$ixsex, 
                   levels = c("Male", "Female"))
meta$ixsex <- relevel(meta$ixsex, 
                    ref="Male")
meta$sars_indid1 <- factor(meta$sars, 
                          levels = c("Negative", "Positive"))


hts_cnet_filt_sym$age_indid1 <- factor(hts_cnet_filt_sym$age_indid1, 
                       levels = c("<5", "5-12", "13-17", "18-34","35-59", "60"), 
                       labels = c("<5", "5-12", "13-17", "18-34","35-59","\u226560"))
hts_cnet_filt_sym$age_indid2 <- factor(hts_cnet_filt_sym$age_indid2, 
                                       levels = c("<5", "5-12", "13-17", "18-34","35-59", "60"), 
                                       labels = c("<5", "5-12", "13-17", "18-34","35-59","\u226560"))
hts_cnet_filt_sym$sars_indid1 <- factor(hts_cnet_filt_sym$sars_indid1, 
                                       levels = c("Negative", "Positive"))
hts_cnet_filt_sym$sars_indid2 <- factor(hts_cnet_filt_sym$sars_indid2, 
                                        levels = c("Negative", "Positive"))

#Combined contact features for household ---

#Number of nodes
hts_cnet_filt_sym <- hts_cnet_filt_sym %>% 
  group_by(hh) %>% 
  mutate(nodes = n_distinct(indid1))

#Full deployment period
contact_summary_hh<- hts_cnet_filt_sym %>% 
  group_by(hh, hcir, nodes) %>% 
  summarize(contact_freq_full = n(),
            contact_dur_full = sum(duration_sec),
            contact_days_full = n_distinct(date))

#Normalized duration based on nodes
contact_summary_hh$contact_dur_norm_full <- contact_summary_hh$contact_dur_full/(contact_summary_hh$nodes*(contact_summary_hh$nodes-1))

#Overall contact summary, includes index ----
ind_daily_contacts_inc_index <- hts_cnet_filt_sym %>%
  group_by(hh, indid1, date) %>%
  summarize(ind_contact_freq_p_day = n(),
            ind_contact_dur_p_day = sum(duration_sec),
            ind_contact_freq_full = n(),
            ind_ave_dur_p_day = ind_contact_dur_p_day/ind_contact_freq_p_day,
            ind_contact_freq_median_p_day = median(ind_contact_freq_p_day),
            members_excluded = max(members_excluded)) %>% 
  ungroup() %>% 
  complete(indid1,
           fill = list(ind_contact_freq_p_day = 0, ind_contact_dur_p_day = 0, 
                       ind_contact_freq_full = 0, 
                       ind_ave_dur_p_day = 0, ind_contact_freq_median_p_day = 0)) %>% 
  mutate(hh = substr(indid1, 1, 4)) %>% 
  group_by(hh) %>% 
  mutate(members_excluded = max(members_excluded, na.rm = TRUE)) %>% 
  ungroup()

#Obtain number of days 
ind_days_nodes <- hts_cnet_filt_sym %>% 
  group_by(indid1) %>% 
  summarise(ind_contact_days_full = n_distinct(date))

#Merge into daily summary
ind_daily_contacts_inc_index <- left_join(ind_daily_contacts_inc_index, ind_days_nodes)

#Calculate daily values
ind_contact_summary_inc_index <- ind_daily_contacts_inc_index %>%
  group_by(indid1) %>%
  summarize(ind_contact_days_full = sum(ind_contact_days_full),
            ind_contact_freq_full = sum(ind_contact_freq_full),
            ind_contact_dur_full = sum(ind_contact_dur_p_day/60),
            ind_contact_ave_dur_full = if_else(ind_contact_dur_full!=0, (ind_contact_dur_full/ind_contact_freq_full), 0),
            ind_contact_freq_mean_p_day = mean(ind_contact_freq_p_day),
            ind_contact_freq_median_p_day = median(ind_contact_freq_p_day),
            ind_contact_dur_mean_p_day = mean(ind_contact_dur_p_day/60),
            ind_contact_dur_median_p_day = median(ind_contact_dur_p_day/60),
            ind_ave_dur_mean_p_day = mean(ind_ave_dur_p_day/60),
            ind_ave_dur_median_p_day = median(ind_ave_dur_p_day/60),
            ind_contact_freq_max = max(ind_contact_freq_p_day, na.rm = TRUE),
            ind_contact_dur_max = max(ind_contact_dur_p_day/60, na.rm = TRUE),
            ind_contact_ave_freq_p_day = if_else(ind_contact_freq_full!=0, ind_contact_freq_full/ind_contact_days_full, 0),
            ind_contact_cum_p_day = if_else(ind_contact_dur_full!=0, ind_contact_dur_full/ind_contact_days_full, 0),
            members_excluded = max(members_excluded, na.rm = TRUE))
rm(ind_daily_contacts_inc_index)

#Individual analysis - contacts with index ----
#Obtain contact features on individual level 
#Summarize individual contact patterns 
ind_daily_contacts <- hts_cnet_filt_sym %>% 
  mutate(index_contact = if_else(grepl("-001", indid2), 1, 0)) %>% 
  group_by(hh, indid1,index_contact, date) %>% 
  summarize(ind_contact_freq_p_day = n(),
            ind_contact_dur_p_day = sum(duration_sec),
            ind_contact_freq_full = n(),
            ind_ave_dur_p_day = ind_contact_dur_p_day/ind_contact_freq_p_day,
            ind_contact_freq_median_p_day = median(ind_contact_freq_p_day),
            members_excluded = max(members_excluded)) %>%
  ungroup() %>% 
  complete(indid1,index_contact,
           fill = list(ind_contact_freq_p_day = 0, ind_contact_dur_p_day = 0, 
                       ind_contact_freq_full = 0, 
                       ind_ave_dur_p_day = 0, ind_contact_freq_median_p_day = 0)) %>% 
  filter(index_contact==1 & grepl("-001", indid1)==FALSE) %>% 
  mutate(hh = substr(indid1, 1, 4)) %>% 
  group_by(hh) %>% 
  mutate(members_excluded = max(members_excluded, na.rm = TRUE)) %>% 
  ungroup()

#Merge into daily summary
ind_daily_contacts <- left_join(ind_daily_contacts, ind_days_nodes)

#Calculate daily values
ind_contact_summary <- ind_daily_contacts %>% 
  group_by(indid1) %>% 
  summarize(ind_contact_days_full = sum(ind_contact_days_full),
            ind_contact_freq_full = sum(ind_contact_freq_full),
            ind_contact_dur_full = sum(ind_contact_dur_p_day/60),
            ind_contact_ave_dur_full = if_else(ind_contact_dur_full!=0, (ind_contact_dur_full/ind_contact_freq_full), 0),
            ind_contact_freq_mean_p_day = mean(ind_contact_freq_p_day),
            ind_contact_freq_median_p_day = median(ind_contact_freq_p_day),
            ind_contact_dur_mean_p_day = mean(ind_contact_dur_p_day/60),
            ind_contact_dur_median_p_day = median(ind_contact_dur_p_day/60),
            ind_ave_dur_mean_p_day = mean(ind_ave_dur_p_day/60),
            ind_ave_dur_median_p_day = median(ind_ave_dur_p_day/60),
            ind_contact_freq_max = max(ind_contact_freq_p_day, na.rm = TRUE),
            ind_contact_dur_max = max(ind_contact_dur_p_day/60, na.rm = TRUE),
            ind_contact_ave_freq_p_day = if_else(ind_contact_freq_full!=0, ind_contact_freq_full/ind_contact_days_full, 0),
            ind_contact_cum_p_day = if_else(ind_contact_dur_full!=0, ind_contact_dur_full/ind_contact_days_full, 0),
            members_excluded = max(members_excluded, na.rm = TRUE))

rm(ind_daily_contacts, ind_days_nodes)

#Grouped analysis - Contact parameters with pos individuals ----
grp_daily_contacts <- hts_cnet_filt_sym %>% 
  group_by(hh, indid1, sars_indid2, pair_sars, hcir, date) %>%
  summarize(grp_contact_freq_p_day = n(),
            grp_contact_dur_p_day = sum(duration_sec),
            grp_contact_freq_full = n(),
            grp_ave_dur_p_day = grp_contact_dur_p_day/grp_contact_freq_p_day,
            grp_contact_freq_median_p_day = median(grp_contact_freq_p_day),
            members_excluded = max(members_excluded)) %>% 
  ungroup() %>% 
  complete(indid1, sars_indid2,
           fill = list(grp_contact_freq_p_day = 0, grp_contact_dur_p_day = 0, 
                       grp_contact_freq_full = 0, 
                       grp_ave_dur_p_day = 0, grp_contact_freq_median_p_day = 0)) %>% 
  filter(sars_indid2=="Positive") %>% 
  mutate(hh = substr(indid1, 1, 4)) %>% 
  group_by(hh) %>% 
  mutate(members_excluded = max(members_excluded, na.rm = TRUE)) %>% 
  ungroup()

#Obtain number of days and nodes
grp_days_nodes <- hts_cnet_filt_sym %>% 
  group_by(indid1) %>% 
  summarise(grp_contact_days_full = n_distinct(date),
            nodes_in_contact_d_grp = n_distinct(indid2))

#Merge into daily summary
grp_daily_contacts <- left_join(grp_daily_contacts, grp_days_nodes)

#Calculate daily values
grp_contact_summary <- grp_daily_contacts %>% 
  group_by(indid1) %>% 
  summarize(grp_contact_days_full = sum(grp_contact_days_full),
            grp_contact_freq_full = sum(grp_contact_freq_full),
            grp_contact_dur_full = sum(grp_contact_dur_p_day),
            grp_contact_ave_dur_full = if_else(grp_contact_dur_full!=0, grp_contact_dur_full/grp_contact_freq_full, 0),
            grp_contact_freq_mean_p_day = mean(grp_contact_freq_p_day),
            grp_contact_freq_median_p_day = median(grp_contact_freq_p_day),
            grp_contact_dur_mean_p_day = mean(grp_contact_dur_p_day),
            grp_contact_dur_median_p_day = median(grp_contact_dur_p_day),
            grp_ave_dur_mean_p_day = mean(grp_ave_dur_p_day),
            grp_ave_dur_median_p_day = median(grp_ave_dur_p_day),grp_contact_freq_max = max(grp_contact_freq_p_day, na.rm = TRUE),
            grp_contact_dur_max = max(grp_contact_dur_p_day/60, na.rm = TRUE),
            grp_contact_freq_max = max(grp_contact_freq_p_day, na.rm = TRUE),
            grp_contact_ave_freq_p_day = if_else(grp_contact_freq_full!=0, grp_contact_freq_full/grp_contact_days_full, 0),
            grp_contact_cum_p_day = if_else(grp_contact_dur_full!=0, grp_contact_dur_full/grp_contact_days_full, 0),
            nodes_in_contact_d_grp = max(nodes_in_contact_d_grp, na.rm = TRUE),
            members_excluded = max(members_excluded))
rm(grp_daily_contacts, grp_days_nodes)

#Prepare for analysis

#Combine with metadata
#Excluding index
ind_contact_summary<-left_join(ind_contact_summary, meta, by=c("indid1"="indid"))

#Including index
ind_contact_summary_inc_index <-full_join(ind_contact_summary_inc_index, meta, by=c("indid1"="indid"))

#Grouped analysis
grp_contact_summary <-left_join(grp_contact_summary, meta, by=c("indid1"="indid"))

#Fill variant across hh (grouped analysis)
grp_contact_summary <- grp_contact_summary %>% 
  group_by(hhid) %>% 
  mutate(ixesarsvarf1 = first(na.omit(ixesarsvarf1))) %>% 
  ungroup()

#---- Analysis 

#Summary of overall contact parameters ----
cont_overall_sum <- ind_contact_summary_inc_index %>% 
rename("Site" = site,
       "Index" = index,
       "Age (years)" = agegrp9,
       "Sex" = sex,
       "SARS-CoV-2 infection" = sars,
       "Median daily duration" = ind_contact_dur_median_p_day,
       "Maximum duration" = ind_contact_dur_max,
       "Median daily frequency" = ind_contact_freq_median_p_day,
       "Maximum frequency" = ind_contact_freq_max,
       "Median average daily duration" = ind_ave_dur_median_p_day,
       "Daily average frequency" = ind_contact_ave_freq_p_day,
       "Cumulative time in contact" = ind_contact_cum_p_day)

tableTwo <- CreateTableOne(vars = c("Median daily duration", "Maximum duration", "Median average daily duration", "Cumulative time in contact",
                                    "Median daily frequency", "Maximum frequency",  "Daily average frequency"),
                           strata = c("Age (years)"), 
                           data = cont_overall_sum,
                           addOverall = TRUE)

print(tableTwo, addOverall = TRUE, 
      nonnormal  = c("Median daily duration", "Maximum duration", 
                     "Median daily frequency", "Maximum frequency", 
                     "Median average daily duration", "Daily average frequency", 
                     "Cumulative time in contact"),
      quote = TRUE, noSpaces = TRUE, contDigits=1)

#By site

cont_overall_sum_k <- cont_overall_sum %>% 
    filter(Site == "Klerksdorp") 

tableTwo <- CreateTableOne(vars = c("Median daily duration", "Maximum duration", "Median average daily duration", "Cumulative time in contact",
                                    "Median daily frequency", "Maximum frequency",  "Daily average frequency"),
                           strata = c("Age (years)"), 
                           data = cont_overall_sum_k,
                           addOverall = TRUE)

print(tableTwo, addOverall = TRUE, 
      nonnormal  = c("Median daily duration", "Maximum duration", 
                     "Median daily frequency", "Maximum frequency", 
                     "Median average daily duration", "Daily average frequency", 
                     "Cumulative time in contact"),
      quote = TRUE, noSpaces = TRUE, contDigits=1)

cont_overall_sum_s <- cont_overall_sum %>% 
         filter(Site == "Soweto") 

tableTwo <- CreateTableOne(vars = c("Median daily duration", "Maximum duration", "Median average daily duration", "Cumulative time in contact",
                                    "Median daily frequency", "Maximum frequency",  "Daily average frequency"),
                           strata = c("Age (years)"), 
                           data = cont_overall_sum_s,
                           addOverall = TRUE)

print(tableTwo, addOverall = TRUE, 
      nonnormal  = c("Median daily duration", "Maximum duration", 
                     "Median daily frequency", "Maximum frequency", 
                     "Median average daily duration", "Daily average frequency", 
                     "Cumulative time in contact"),
      quote = TRUE, noSpaces = TRUE, contDigits=1)

rm(tableTwo)

#As figures

colfunc <-colorRampPalette(c("#006666", "#99e6e6","#8e02b1"))
col_scheme1 <- colfunc(7)

#Combine one dataset with values by site and overall, and for age groups and all ages
cont_overall_sum_oa <-  cont_overall_sum %>% 
  mutate(`Age (years)` = "All")

cont_overall_sum <- rbind(cont_overall_sum, cont_overall_sum_oa)

cont_overall_sum_oa <-  cont_overall_sum %>% 
  mutate(Site = "Overall")

cont_overall_sum <- rbind(cont_overall_sum, cont_overall_sum_oa)

#Plots
#Median daily duration
figdurmed <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=(`Median daily duration`), fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Median daily duration \n(minutes, log scale)\n") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Maximum duration
figdurmax <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=(`Maximum duration`), fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Maximum duration \n(minutes, log scale)\n") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Median average daily duration
figdurave <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=(`Median average daily duration`), fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Median average daily duration \n(minutes, log scale)\n") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Cumulative time in contact
figdurcum <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=(`Cumulative time in contact`), fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Cumulative time in contact \n(minutes, log scale)") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(legend.position="bottom",
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8)) +
  scale_y_continuous(trans='log10') +
  guides(fill = guide_legend(nrow = 1))

#Median daily frequency
figfreqmed <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=`Median daily frequency`, fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Median daily frequency \n(log scale)") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Maximum daily frequency
figfreqmax <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=`Maximum frequency`, fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Maximum frequency \n(log scale)") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Daily average frequency
figfreqave <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=`Daily average frequency`, fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Daily average frequency \n(log scale)") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10')  +
  theme(legend.position="bottom",
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8)) +
  scale_y_continuous(trans='log10') +
  guides(fill = guide_legend(nrow = 1))

ggarrange(figdurmed, figdurmax, figdurave, figdurcum, 
          ncol = 1, nrow = 4,
          heights = c(1,1,1,1.45))
ggsave("Contact Parameters Overall Duration.tiff", width = 16, height = 16, units = "cm", dpi=300)

ggarrange(figfreqmed, figfreqmax, figfreqave, 
          ncol = 1, nrow = 3,
          heights = c(1,1,1.45))
ggsave("Contact Parameters Overall Frequency.tiff", width = 16, height = 16, units = "cm", dpi=300)

rm(cont_overall_sum, cont_overall_sum_oa, cont_overall_sum_k, cont_overall_sum_s, 
   figdurmed, figdurmax, figdurave, figdurcum,
   figfreqmed, figfreqmax, figfreqave)

#Age-based contact matrices ----

#List of ages 
ind_inc_ana_prox_ages <- grp_contact_summary %>% 
  select(indid1, agegrp9) %>% 
  mutate(agegrp9 = relevel(agegrp9, ref="<5")) %>% 
  rename("indid" = "indid1","Age (years)" = "agegrp9")


#Use symmetrical network, renaming age fir graph
hts_cnet_filt_sym_temp <- hts_cnet_filt_sym %>% 
  rename(`Age (years) 2` = age_indid2,
         `Age (years) 1` = age_indid1)

#Normalize based on number of individuals per age group
#Also this will need manual coding if matrices needed in other age groups

#List of individuals with contact data
#Check IDs under ind1
ind1_list <- hts_cnet_filt_sym_temp %>%
  rename(indid = indid1,
         `Age (years)` = `Age (years) 1`) %>% 
  distinct(indid, `Age (years)`)
#Check IDs under id2
ind2_list <- hts_cnet_filt_sym_temp %>%
  rename(indid = indid2,
         `Age (years)` = `Age (years) 2`) %>% 
  distinct(indid, `Age (years)`)

#Merge these two list
ind_list <- rbind(ind1_list, ind2_list)
rm(ind1_list, ind2_list)
#Keep unique values
ind_list <- ind_list %>% 
  distinct() 

#Population count per age group - Klerksdorp
age_counts <- data.frame(ind_list %>% 
                           filter(grepl("J", indid)) %>%   
                           group_by(`Age (years)`) %>% 
                           summarise(count = n()) %>% 
                           rename(age =  `Age (years)`))
age_list <- age_counts$count
names(age_list) <- age_counts$age

for (a in 1:length(age_list)) {
  nam <- paste("age_", names(age_list[a]), sep = "")
  assign(nam, age_list + age_list[a])
  print(names(age_list[a]))
}

pop_matrix_k <- matrix(c(`age_<5`, `age_5-12`, `age_13-17`, `age_18-34`, `age_35-59`, `age_=60`), ncol = 6)
for (a in 1:ncol(pop_matrix_k)) {
  pop_matrix_k[a,a] = pop_matrix_k[a,a]/2
}

#Population count per age group - Soweto
age_counts <- data.frame(ind_list %>% 
                           filter(grepl("S", indid)) %>%   
                           group_by(`Age (years)`) %>% 
                           summarise(count = n()) %>% 
                           rename(age =  `Age (years)`))
age_list <- age_counts$count
names(age_list) <- age_counts$age

for (a in 1:length(age_list)) {
  nam <- paste("age_", names(age_list[a]), sep = "")
  assign(nam, age_list + age_list[a])
  print(names(age_list[a]))
}

pop_matrix_s <- matrix(c(`age_<5`, `age_5-12`, `age_13-17`, `age_18-34`, `age_35-59`, `age_=60`), ncol = 6)
for (a in 1:ncol(pop_matrix_s)) {
  pop_matrix_s[a,a] = pop_matrix_s[a,a]/2
}

#Combined
pop_matrix_comb <- pop_matrix_k + pop_matrix_s

#Normalized matrices (based on combined count in that age group)
#Overall
#Duration
dur_matrix <- as.matrix(hts_cnet_filt_sym_temp %>% 
                          group_by(`Age (years) 1`, `Age (years) 2`) %>% 
                          summarise(norm_dur = sum(duration_sec)/60) %>% 
                          spread(`Age (years) 1`, norm_dur) %>% 
                          select(`<5`, `5-12`, `13-17`, `18-34`, `35-59`, `=60`)) / pop_matrix_comb
row.names(dur_matrix) <- c("<5", "5-12", "13-17", "18-34", "35-59", "\u226560")

#Frequency
freq_matrix <- as.matrix(hts_cnet_filt_sym_temp %>% 
                           group_by(`Age (years) 1`, `Age (years) 2`) %>% 
                           summarise(norm_freq = n()) %>% 
                           spread(`Age (years) 1`, norm_freq) %>% 
                           select(`<5`, `5-12`, `13-17`, `18-34`, `35-59`, `=60`)) / pop_matrix_comb
row.names(freq_matrix) <- c("<5", "5-12", "13-17", "18-34", "35-59", "\u226560")

#By Site
#Klerksdorp
#Duration
dur_matrix_klerk <- as.matrix(hts_cnet_filt_sym_temp %>% 
                                filter(grepl("J", indid1)) %>% 
                                group_by(`Age (years) 1`, `Age (years) 2`) %>% 
                                summarise(norm_dur = (sum(duration_sec)/60)) %>% 
                                spread(`Age (years) 1`, norm_dur) %>% 
                                select(`<5`, `5-12`, `13-17`, `18-34`, `35-59`, `=60`)) / pop_matrix_k
row.names(dur_matrix_klerk) <- c("<5", "5-12", "13-17", "18-34", "35-59", "\u226560")

#Frequency
freq_matrix_klerk <- as.matrix(hts_cnet_filt_sym_temp %>% 
                                 filter(grepl("J", indid1)) %>%  
                                 group_by(`Age (years) 1`, `Age (years) 2`) %>% 
                                 summarise(norm_freq = n()) %>% 
                                 spread(`Age (years) 1`, norm_freq) %>% 
                                 select(`<5`, `5-12`, `13-17`, `18-34`, `35-59`, `=60`)) / pop_matrix_k
row.names(freq_matrix_klerk) <- c("<5", "5-12", "13-17", "18-34", "35-59", "\u226560")

#Soweto
#Duration
dur_matrix_sow <- as.matrix(hts_cnet_filt_sym_temp %>% 
                              filter(grepl("S", indid1)) %>% 
                              group_by(`Age (years) 1`, `Age (years) 2`) %>% 
                              summarise(norm_dur = (sum(duration_sec)/60)) %>% 
                              spread(`Age (years) 1`, norm_dur) %>% 
                              select(`<5`, `5-12`, `13-17`, `18-34`, `35-59`, `=60`)) / pop_matrix_s
row.names(dur_matrix_klerk) <- c("<5", "5-12", "13-17", "18-34", "35-59", "\u226560")

#Frequency
freq_matrix_sow <- as.matrix(hts_cnet_filt_sym_temp %>% 
                               filter(grepl("S", indid1)) %>%  
                               group_by(`Age (years) 1`, `Age (years) 2`) %>%
                               summarise(norm_freq = n()) %>% 
                               spread(`Age (years) 1`, norm_freq) %>% 
                               select(`<5`, `5-12`, `13-17`, `18-34`, `35-59`, `=60`)) / pop_matrix_s
row.names(freq_matrix_sow) <- c("<5", "5-12", "13-17", "18-34", "35-59", "\u226560")

#Define color scheme used
colfunc <-colorRampPalette(c("#006666", "#99e6e6","#440154"))
col_scheme <- colfunc(1000)

#Output duration matrix - Overall
tiff(file=paste0("Contact Matrices.tiff"),
     width=16, height=24, unit = "cm", res=600)
layout(matrix(c(1,2,3,4,5,6), nrow = 3))

par(mar = c(2.5, 2.5, 2.5, 3)) 
image(dur_matrix,  xaxt= "n", yaxt= "n", frame.plot=F,
      col = col_scheme)
mtext( "A", font=2, side=3, at=0)
title(xlab = "Age (years)", line=1, cex.lab=1)
title(ylab = "Age (years)", line=1, cex.lab=1)
axis(1, at=seq(0,1,0.2), labels=rownames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
axis(2, at=seq(0,1,0.2), labels=colnames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
image.plot(dur_matrix,  xaxt= "n", yaxt= "n", frame.plot=F,
           col = col_scheme, 
           axis.args = list(cex.axis=0.75, tck=-0.025, mgp=c(0,0.5,0)),
           horizontal = FALSE,
           legend.width = 1,
           legend.mar = 3,
           legend.only = TRUE) 

image(freq_matrix,  xaxt= "n", yaxt= "n", frame.plot=F,
      col = col_scheme)
mtext( "D", font=2, side=3, at=0)
title(xlab = "Age (years)", line=1, cex.lab=1)
title(ylab = "Age (years)", line=1, cex.lab=1)
axis(1, at=seq(0,1,0.2), labels=rownames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
axis(2, at=seq(0,1,0.2), labels=colnames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
image.plot(freq_matrix,  xaxt= "n", yaxt= "n", frame.plot=F,
           col = col_scheme, 
           axis.args = list(cex.axis=0.75, tck=-0.025, mgp=c(0,0.5,0)),
           horizontal = FALSE,
           legend.width = 1,
           legend.mar = 3,
           legend.only = TRUE) 

image(dur_matrix_klerk,  xaxt= "n", yaxt= "n", frame.plot=F,
      col = col_scheme)
mtext( "B", font=2, side=3, at=0)
title(xlab = "Age (years)", line=1, cex.lab=1)
title(ylab = "Age (years)", line=1, cex.lab=1)
axis(1, at=seq(0,1,0.2), labels=rownames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
axis(2, at=seq(0,1,0.2), labels=colnames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
image.plot(dur_matrix_klerk,  xaxt= "n", yaxt= "n", frame.plot=F,
           col = col_scheme, 
           axis.args = list(cex.axis=0.75, tck=-0.025, mgp=c(0,0.5,0)),
           horizontal = FALSE,
           legend.width = 1,
           legend.mar = 3,
           legend.only = TRUE) 

image(freq_matrix_klerk,  xaxt= "n", yaxt= "n", frame.plot=F,
      col = col_scheme)
mtext( "E", font=2, side=3, at=0)
title(xlab = "Age (years)", line=1, cex.lab=1)
title(ylab = "Age (years)", line=1, cex.lab=1)
axis(1, at=seq(0,1,0.2), labels=rownames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
axis(2, at=seq(0,1,0.2), labels=colnames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
image.plot(freq_matrix_klerk,  xaxt= "n", yaxt= "n", frame.plot=F,
           col = col_scheme, 
           axis.args = list(cex.axis=0.75, tck=-0.025, mgp=c(0,0.5,0)),
           horizontal = FALSE,
           legend.width = 1,
           legend.mar = 3,
           legend.only = TRUE) 

image(dur_matrix_sow,  xaxt= "n", yaxt= "n", frame.plot=F,
      col = col_scheme)
mtext( "C", font=2, side=3, at=0)
title(xlab = "Age (years)", line=1, cex.lab=1)
title(ylab = "Age (years)", line=1, cex.lab=1)
axis(1, at=seq(0,1,0.2), labels=rownames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
axis(2, at=seq(0,1,0.2), labels=colnames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
image.plot(dur_matrix_sow,  xaxt= "n", yaxt= "n", frame.plot=F,
           col = col_scheme, 
           axis.args = list(cex.axis=0.75, tck=-0.025, mgp=c(0,0.5,0)),
           horizontal = FALSE,
           legend.width = 1,
           legend.mar = 3,
           legend.only = TRUE) 

image(freq_matrix_sow,  xaxt= "n", yaxt= "n", frame.plot=F,
      col = col_scheme)
mtext( "F", font=2, side=3, at=0)
title(xlab = "Age (years)", line=1, cex.lab=1)
title(ylab = "Age (years)", line=1, cex.lab=1)
axis(1, at=seq(0,1,0.2), labels=rownames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
axis(2, at=seq(0,1,0.2), labels=colnames(dur_matrix), lty = 0, cex.axis=0.75, tck=-0.0025, mgp=c(-0.75,-0.05,0))
image.plot(freq_matrix_sow,  xaxt= "n", yaxt= "n", frame.plot=F,
           col = col_scheme, 
           axis.args = list(cex.axis=0.75, tck=-0.025, mgp=c(0,0.5,0)),
           horizontal = FALSE,
           legend.width = 1,
           legend.mar = 3,
           legend.only = TRUE) 

dev.off()

rm("hts_cnet_filt_sym_temp", "a", "age_<5", "age_5-12", "age_13-17", "age_18-34", "age_35-59", 
   "age_=60", "age_counts" ,"age_list", "col_scheme", "col_scheme1","colfunc",
   "dur_matrix", "dur_matrix_klerk", "dur_matrix_sow", 
   "ind_inc_ana_prox_ages", "ind_list", "nam","pop_matrix_comb", 
   "pop_matrix_k" , "pop_matrix_s", freq_matrix, freq_matrix_klerk, freq_matrix_sow
   )

#Regression ----

ind_contact_summary$sars <- factor(ind_contact_summary$sars)
ind_contact_summary$agegrp9 <- relevel(ind_contact_summary$agegrp9, ref="18-34")

#Due to privacy concerns not all meta data included. Only final model shown
  
#Univariate analysis
#Generate vector with list of variables to investigate
ind_var_list <- c("site", "ixagegrp9", "ixminctcat2", "ixesarsvarf1", 
                  "agegrp9", "sex", "sleep_room_ix", "cared_by_ix", 
                  "ind_contact_dur_median_p_day", "ind_contact_dur_max", "ind_ave_dur_median_p_day", "ind_contact_cum_p_day",
                  "ind_contact_freq_median_p_day", "ind_contact_freq_max","ind_contact_ave_freq_p_day")

#Frequency table
ind_contact_summary %>% 
  select(c(ind_var_list, "sars")) %>% # keep only columns of interest
  tbl_summary(     
    by = sars,
    type = list(c(sleep_room_ix,cared_by_ix) ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{median} ({p25}-{p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"),
    percent="row",
    label  = list(                                              # display labels for column names
      site   ~ "Site",                           
      ixagegrp9 ~ "Index Age (years)",
      ixminctcat2 ~ "Index Minimum Ct value",
      agegrp9 ~  "Contact Age (years)",
      sex ~  "Contact Sex",
      ixesarsvarf1 ~ "SARS-CoV-2 variant",
      sleep_room_ix ~ "Sleep in same room as index",
      cared_by_ix ~ "Cared for by index",
      ind_contact_dur_median_p_day ~  "Median daily duration",
      ind_contact_dur_max ~ "Maximum duration",
      ind_ave_dur_median_p_day ~ "Median average daily duration",
      ind_contact_cum_p_day ~ "Cumulative time in contact",
      ind_contact_freq_median_p_day ~  "Median daily frequency",
      ind_contact_freq_max ~ "Maximum frequency",
      ind_contact_ave_freq_p_day ~ "Daily average frequency"),) 

#To complete full univariate analysis, define list of variables and function to perform model
univ_res <- lapply(ind_var_list,
                   
                   function(var) {
                     
                     formula    <- as.formula(paste("sars ~", var, "+ (1 | site) + (1 | hhid)"))
                     res.logist <- glmer(formula, data = ind_contact_summary,
                                         family = binomial("logit"))
                     data.frame(tidy(res.logist,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
                       mutate_if(is.numeric, round, 4)
                   })

print(univ_res, quote=TRUE, noSpaces = TRUE)
rm(univ_res)

#Build model: Individual-level ----

#Complete model excluding contact parameters
ind_fin_mod <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                       agegrp9 + sex + 
                       (1 | site) + (1 | hhid), 
                     data = ind_contact_summary,
                     family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median daily duration of contacts per day with index (ind_contact_dur_median_p_day)
ind_fin_mod_mdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contact_dur_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Maximum duration of contacts with index (ind_contact_dur_max)
ind_fin_mod_mxdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contact_dur_max + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mxdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median average daily duration of contacts with index (ind_ave_dur_median_p_day)
ind_fin_mod_avdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                             agegrp9 + sex + 
                             ind_ave_dur_median_p_day + 
                             (1 | site) + (1 | hhid), 
                           data = ind_contact_summary,
                           family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_avdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Cumulative time in contact with index (ind_contact_cum_p_day)
ind_fin_mod_cdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                             agegrp9 + sex + 
                             ind_contact_cum_p_day + 
                             (1 | site) + (1 | hhid), 
                           data = ind_contact_summary,
                           family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_cdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median frequency of contacts per day with index (ind_contact_freq_median_p_day)
ind_fin_mod_mfeq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contact_freq_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Maximum frequency of contacts per day with index (ind_contact_freq_max)
ind_fin_mod_mxfeq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                             ind_contact_freq_max + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mxfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Daily average duration of contacts overall with index (ind_contact_ave_freq_p_day)
ind_fin_mod_afreq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                             ind_contact_ave_freq_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_afreq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)


#Build model: Grouped analysis ----

grp_contact_summary$sars <- factor(grp_contact_summary$sars)
grp_contact_summary$agegrp9 <- relevel(grp_contact_summary$agegrp9, ref="18-34")

#Univariate analysis
#Generate vector with list of variables to investigate
grp_var_list <- c("site", "agegrp9", "sex", "bmicat", "smokecignow1", "ixesarsvarf1",
                  "grp_contact_dur_median_p_day", "grp_contact_dur_max", "grp_ave_dur_median_p_day", "grp_contact_cum_p_day",
                  "grp_contact_freq_median_p_day", "grp_contact_freq_max","grp_contact_ave_freq_p_day")

#Frequency table
grp_contact_summary %>% 
  select(c(grp_var_list, "sars")) %>% # keep only columns of interest
  tbl_summary(     
    by = sars,
    percent="row",
    statistic = list(
      all_continuous() ~ "{median} ({p25}-{p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    label  = list(                                              # display labels for column names
      site   ~ "Site",                           
      agegrp9 ~  "Contact Age (years)",
      sex ~  "Contact Sex",
      bmicat ~ "Body mass index",
      smokecignow1 ~ "Current smoking",
      ixesarsvarf1 ~ "SARS-CoV-2 variant",
      grp_contact_dur_median_p_day ~  "Median daily duration",
      grp_contact_dur_max ~ "Maximum duration",
      grp_ave_dur_median_p_day ~ "Median average daily duration",
      grp_contact_cum_p_day ~ "Cumulative time in contact",
      grp_contact_freq_median_p_day ~  "Median daily frequency",
      grp_contact_freq_max ~ "Maximum frequency",
      grp_contact_ave_freq_p_day ~ "Daily average frequency"),) 

#To complete full univariate analysis, define list of variables and function to perform model

univ_res <- lapply(grp_var_list,
                   
                   function(var) {
                     
                     formula    <- as.formula(paste("sars ~", var, "+ (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp)"))
                     res.logist <- glmer(formula, data = grp_contact_summary,
                                         family = binomial("logit"))
                     data.frame(tidy(res.logist,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
                       mutate_if(is.numeric, round, 4)
                   })
print(univ_res, quote=TRUE, noSpaces = TRUE)
rm(univ_res)

#Build model: group-level ----

#Complete model excluding contact parameters
grp_fin_mod <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                       (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                     data = grp_contact_summary,
                     family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median daily duration of contacts per day with cases (grp_contact_dur_median_p_day)
grp_fin_mod_mdur <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                       grp_contact_dur_median_p_day + 
                       (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                     data = grp_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_mdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Maximum duration of contacts with cases (grp_contact_dur_max)
grp_fin_mod_mxdur <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                             grp_contact_dur_max + 
                             (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                           data = grp_contact_summary,
                           family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_mxdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median average daily duration of contacts with cases (grp_ave_dur_median_p_day)
grp_fin_mod_avdur <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                             grp_ave_dur_median_p_day + 
                             (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                           data = grp_contact_summary,
                           family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_avdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Cumulative time in contact with cases (grp_contact_cum_p_day)
grp_fin_mod_cdur <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 +
                            grp_contact_cum_p_day + 
                            (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                          data = grp_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_cdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median frequency of contacts per day with cases (grp_contact_freq_median_p_day)
grp_fin_mod_mfeq <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                            grp_contact_freq_median_p_day + 
                            (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                          data = grp_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_mfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Maximum frequency of contacts per day with cases (grp_contact_freq_max)
grp_fin_mod_mxfeq <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                             grp_contact_freq_max + 
                             (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                           data = grp_contact_summary,
                           family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_mxfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Daily average frequency of contacts overall with cases (grp_contact_ave_freq_p_day)
grp_fin_mod_afreq <- glmer(sars ~ agegrp9 + bmicat + smokecignow1 +  ixesarsvarf1 + 
                             grp_contact_ave_freq_p_day + 
                             (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                           data = grp_contact_summary,
                           family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_afreq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Assess association with Wilcoxon rank-sum test ----

#Median daily duration
ind_contact_dur_median_p_day <- wilcox.test( ind_contact_dur_median_p_day ~ sars, data = ind_contact_summary,
                                             paired = FALSE, 
                                             exact = FALSE)
#Maximum duration
ind_contact_dur_max <- wilcox.test( ind_contact_dur_max ~ sars, data = ind_contact_summary,
                                    paired = FALSE, 
                                    exact = FALSE)
#Median average daily duration 
ind_ave_dur_median_p_day <- wilcox.test( ind_ave_dur_median_p_day ~ sars, data = ind_contact_summary,
                                         paired = FALSE, 
                                         exact = FALSE)
#Cumulative time in contact
ind_contact_cum_p_day <- wilcox.test( ind_contact_cum_p_day ~ sars, data = ind_contact_summary,
                                      paired = FALSE, 
                                      exact = FALSE)
#Median daily frequency
ind_contact_freq_median_p_day <- wilcox.test( ind_contact_freq_median_p_day ~ sars, data = ind_contact_summary,
                                              paired = FALSE, 
                                              exact = FALSE)
#Maximum frequency
ind_contact_freq_max <- wilcox.test( ind_contact_freq_max ~ sars, data = ind_contact_summary,
                                     paired = FALSE, 
                                     exact = FALSE)
#Daily average frequency
ind_contact_ave_freq_p_day <- wilcox.test( ind_contact_ave_freq_p_day ~ sars, data = ind_contact_summary,
                                           paired = FALSE, 
                                           exact = FALSE)

#Median daily duration
grp_contact_dur_median_p_day <- wilcox.test(  grp_contact_dur_median_p_day ~ sars, data = grp_contact_summary,
                                              paired = FALSE, 
                                              exact = FALSE)
#Maximum duration
grp_contact_dur_max <- wilcox.test( grp_contact_dur_max ~ sars, data = grp_contact_summary,
                                    paired = FALSE, 
                                    exact = FALSE)
#Median average daily duration
grp_ave_dur_median_p_day  <- wilcox.test(grp_ave_dur_median_p_day ~ sars, data = grp_contact_summary,
                                         paired = FALSE, 
                                         exact = FALSE)
#Cumulative time in contact
grp_contact_cum_p_day <- wilcox.test( grp_contact_cum_p_day ~ sars, data = grp_contact_summary,
                                      paired = FALSE, 
                                      exact = FALSE)
#Median daily frequency
grp_contact_freq_median_p_day <- wilcox.test( grp_contact_freq_median_p_day ~ sars, data = grp_contact_summary,
                                              paired = FALSE, 
                                              exact = FALSE)
#Maximum frequency
grp_contact_freq_max <- wilcox.test( grp_contact_freq_max ~ sars, data = grp_contact_summary,
                                     paired = FALSE, 
                                     exact = FALSE)
#Daily average frequency
grp_contact_ave_freq_p_day <- wilcox.test( grp_contact_ave_freq_p_day ~ sars, data = grp_contact_summary,
                                           paired = FALSE, 
                                           exact = FALSE)


wilcox_results_complete <- data.frame(rbind(ind_contact_dur_median_p_day, ind_contact_dur_max, ind_ave_dur_median_p_day, 
                        ind_contact_cum_p_day, ind_contact_freq_median_p_day, ind_contact_freq_max, 
                        ind_contact_ave_freq_p_day, 
                        grp_contact_dur_median_p_day, grp_contact_dur_max, grp_ave_dur_median_p_day, 
                        grp_contact_cum_p_day, grp_contact_freq_median_p_day, grp_contact_freq_max, 
                        grp_contact_ave_freq_p_day)) %>% 
  select(p.value) %>% 
  rename("P value complete analysis"="p.value")

print(wilcox_results_complete, quote = TRUE)

#Sensitivity analysis: repeat results with only households without excluded members ----

ind_contact_summary_excl <- ind_contact_summary %>% 
  filter(members_excluded == 0) 

grp_contact_summary_excl <- grp_contact_summary %>% 
  filter(members_excluded == 0) 

#Sensitivity analysis: Logistic regression ----
#Individual analysis (sensitivity) ----

#Frequency table
ind_contact_summary_excl %>% 
  select(c(ind_var_list, "sars")) %>% # keep only columns of interest
  tbl_summary(     
    by = sars,
    type = list(c(sleep_room_ix,cared_by_ix) ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{median} ({p25}-{p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"),
    percent="row",
    label  = list(                                              # display labels for column names
      site   ~ "Site",                           
      ixagegrp9 ~ "Index Age (years)",
      ixminctcat2 ~ "Index Minimum Ct value",
      agegrp9 ~  "Contact Age (years)",
      sex ~  "Contact Sex",
      ixesarsvarf1 ~ "SARS-CoV-2 variant",
      sleep_room_ix ~ "Sleep in same room as index",
      cared_by_ix ~ "Cared for by index",
      ind_contact_dur_median_p_day ~  "Median daily duration",
      ind_contact_dur_max ~ "Maximum duration",
      ind_ave_dur_median_p_day ~ "Median average daily duration",
      ind_contact_cum_p_day ~ "Cumulative time in contact",
      ind_contact_freq_median_p_day ~  "Median daily frequency",
      ind_contact_freq_max ~ "Maximum frequency",
      ind_contact_ave_freq_p_day ~ "Daily average frequency"),) 

univ_res <- lapply(ind_var_list,
                   
                   function(var) {
                     
                     formula    <- as.formula(paste("sars ~", var, "+ (1 | site) + (1 | hhid)"))
                     res.logist <- glmer(formula, data = ind_contact_summary_excl,
                                         family = binomial("logit"))
                     data.frame(tidy(res.logist,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
                       mutate_if(is.numeric, round, 4)
                   })

print(univ_res, quote=TRUE, noSpaces = TRUE)

#Complete model excluding contact parameters
ind_fin_mod <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                       agegrp9 + sex + 
                       (1 | site) + (1 | hhid), 
                     data = ind_contact_summary_excl,
                     family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median daily duration of contacts per day with index (ind_contact_dur_median_p_day)
ind_fin_mod_mdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contact_dur_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary_excl,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Maximum duration of contacts with index (ind_contact_dur_max)
ind_fin_mod_mxdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                             agegrp9 + sex + 
                             ind_contact_dur_max + 
                             (1 | site) + (1 | hhid), 
                           data = ind_contact_summary_excl,
                           family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mxdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median average daily duration of contacts with index (ind_ave_dur_median_p_day)
ind_fin_mod_avdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                             agegrp9 + sex + 
                             ind_ave_dur_median_p_day + 
                             (1 | site) + (1 | hhid), 
                           data = ind_contact_summary_excl,
                           family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_avdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Cumulative time in contact with index (ind_contact_cum_p_day)
ind_fin_mod_cdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contact_cum_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary_excl,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_cdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median frequency of contacts per day with index (ind_contact_freq_median_p_day)
ind_fin_mod_mfeq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contact_freq_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary_excl,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Maximum frequency of contacts per day with index (ind_contact_freq_max)
ind_fin_mod_mxfeq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                             agegrp9 + sex + 
                             ind_contact_freq_max + 
                             (1 | site) + (1 | hhid), 
                           data = ind_contact_summary_excl,
                           family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mxfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Daily average duration of contacts overall with index (ind_contact_ave_freq_p_day)
ind_fin_mod_afreq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                             agegrp9 + sex + 
                             ind_contact_ave_freq_p_day + 
                             (1 | site) + (1 | hhid), 
                           data = ind_contact_summary_excl,
                           family = binomial("logit"))

print(data.frame(tidy(ind_fin_mod_afreq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)


#Sensitivity analysis: wilcoxon rank sum ----

#Median daily duration
ind_contact_dur_median_p_day <- wilcox.test( ind_contact_dur_median_p_day ~ sars, data = ind_contact_summary_excl,
                                             paired = FALSE, 
                                             exact = FALSE)
#Maximum duration
ind_contact_dur_max <- wilcox.test( ind_contact_dur_max ~ sars, data = ind_contact_summary_excl,
                                    paired = FALSE, 
                                    exact = FALSE)
#Median average daily duration 
ind_ave_dur_median_p_day <- wilcox.test( ind_ave_dur_median_p_day ~ sars, data = ind_contact_summary_excl,
                                         paired = FALSE, 
                                         exact = FALSE)
#Cumulative time in contact
ind_contact_cum_p_day <- wilcox.test( ind_contact_cum_p_day ~ sars, data = ind_contact_summary_excl,
                                      paired = FALSE, 
                                      exact = FALSE)
#Median daily frequency
ind_contact_freq_median_p_day <- wilcox.test( ind_contact_freq_median_p_day ~ sars, data = ind_contact_summary_excl,
                                              paired = FALSE, 
                                              exact = FALSE)
#Maximum frequency
ind_contact_freq_max <- wilcox.test( ind_contact_freq_max ~ sars, data = ind_contact_summary_excl,
                                     paired = FALSE, 
                                     exact = FALSE)
#Daily average frequency
ind_contact_ave_freq_p_day <- wilcox.test( ind_contact_ave_freq_p_day ~ sars, data = ind_contact_summary_excl,
                                           paired = FALSE, 
                                           exact = FALSE)

#Median daily duration
grp_contact_dur_median_p_day <- wilcox.test(  grp_contact_dur_median_p_day ~ sars, data = grp_contact_summary_excl,
                                              paired = FALSE, 
                                              exact = FALSE)
#Maximum duration
grp_contact_dur_max <- wilcox.test( grp_contact_dur_max ~ sars, data = grp_contact_summary_excl,
                                    paired = FALSE, 
                                    exact = FALSE)
#Median average daily duration
grp_ave_dur_median_p_day  <- wilcox.test(grp_ave_dur_median_p_day ~ sars, data = grp_contact_summary_excl,
                                         paired = FALSE, 
                                         exact = FALSE)
#Cumulative time in contact
grp_contact_cum_p_day <- wilcox.test( grp_contact_cum_p_day ~ sars, data = grp_contact_summary_excl,
                                      paired = FALSE, 
                                      exact = FALSE)
#Median daily frequency
grp_contact_freq_median_p_day <- wilcox.test( grp_contact_freq_median_p_day ~ sars, data = grp_contact_summary_excl,
                                              paired = FALSE, 
                                              exact = FALSE)
#Maximum frequency
grp_contact_freq_max <- wilcox.test( grp_contact_freq_max ~ sars, data = grp_contact_summary_excl,
                                     paired = FALSE, 
                                     exact = FALSE)
#Daily average frequency
grp_contact_ave_freq_p_day <- wilcox.test( grp_contact_ave_freq_p_day ~ sars, data = grp_contact_summary_excl,
                                           paired = FALSE, 
                                           exact = FALSE)


wilcox_results_rest <- data.frame(rbind(ind_contact_dur_median_p_day, ind_contact_dur_max, ind_ave_dur_median_p_day, 
                                            ind_contact_cum_p_day, ind_contact_freq_median_p_day, ind_contact_freq_max, 
                                            ind_contact_ave_freq_p_day, 
                                            grp_contact_dur_median_p_day, grp_contact_dur_max, grp_ave_dur_median_p_day, 
                                            grp_contact_cum_p_day, grp_contact_freq_median_p_day, grp_contact_freq_max, 
                                            grp_contact_ave_freq_p_day)) %>% 
  select(p.value) %>% 
  rename("P value restricted analysis"="p.value")

print(wilcox_results_rest, quote = TRUE)

print(cbind(wilcox_results_complete, wilcox_results_rest), quote = FALSE)
