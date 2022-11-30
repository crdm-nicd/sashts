#South African SARS-CoV-2 Household Transmission Study
    #Data used for Kleynhans et al 2023 
      #Household transmission of SARS-CoV-2 from adult index cases living with 
        #and without HIV in South Africa, 2020-2021: a case-ascertained, 
        #prospective observational household transmission study

#Purpose: Describe close proximity evenrs and assess association between close proximity events and SARS-CoV-2 transmission
#Created by: Jackie Kleynhans, National Institute for Communicable Diseases
#Contact: jackiel@nicd.ac.za

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

#Set working directory
setwd("~/COVID-19 HTS/Proximity/Data/Kleynhans 2023 SA-S-HTS Proximity")

#Import data
hts_cnet_filt_sym <- read_csv("sashts_contact_network.csv")
meta <- read_csv("sashts_metadata.csv")

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
meta$sars_indid1 <- factor(meta$sars, 
                          levels = c("Negative", "Positive"))


hts_cnet_filt_sym$age_indid1 <- factor(hts_cnet_filt_sym$age_indid1, 
                       levels = c("<5", "5-12", "13-17", "18-34","35-59", "=60"), 
                       labels = c("<5", "5-12", "13-17", "18-34","35-59","\u226560"))
hts_cnet_filt_sym$age_indid2 <- factor(hts_cnet_filt_sym$age_indid2, 
                                       levels = c("<5", "5-12", "13-17", "18-34","35-59", "=60"), 
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
  group_by(indid1, date) %>%
  summarize(ind_contact_freq_p_day = n(),
            ind_contact_dur_p_day = sum(duration_sec),
            ind_contact_freq_full = n(),
            ind_ave_dur_p_day = ind_contact_dur_p_day/ind_contact_freq_p_day,
            ind_contatc_freq_median_p_day = median(ind_contact_freq_p_day)) %>% 
  ungroup() %>% 
  complete(indid1,
           fill = list(ind_contact_freq_p_day = 0, ind_contact_dur_p_day = 0, 
                       ind_contact_freq_full = 0, 
                       ind_ave_dur_p_day = 0, ind_contatc_freq_median_p_day = 0))

#Obtain number of days 
ind_days_nodes <- hts_cnet_filt_sym %>% 
  group_by(indid1) %>% 
  summarise(ind_contact_days_full = n_distinct(date))

#Merge into daily summary
ind_daily_contacts_inc_index <- left_join(ind_daily_contacts_inc_index, ind_days_nodes)

#Generate summary variables
ind_contact_summary_inc_index <- ind_daily_contacts_inc_index %>%
  group_by(indid1) %>%
  summarize(ind_contact_days_full = sum(ind_contact_days_full),
            ind_contact_freq_full = sum(ind_contact_freq_full),
            ind_contact_dur_full = sum(ind_contact_dur_p_day),
            ind_contact_ave_dur_full = if_else(ind_contact_dur_full!=0, ind_contact_dur_full/ind_contact_freq_full, 0),
            ind_contatc_freq_mean_p_day = mean(ind_contact_freq_p_day),
            ind_contatc_freq_median_p_day = median(ind_contact_freq_p_day),
            ind_contatc_dur_mean_p_day = mean(ind_contact_dur_p_day),
            ind_contatc_dur_median_p_day = median(ind_contact_dur_p_day),
            ind_ave_dur_mean_p_day = mean(ind_ave_dur_p_day),
            ind_ave_dur_median_p_day = median(ind_ave_dur_p_day))
rm(ind_daily_contacts_inc_index)

#Individual analysis - contacts with index 
#Obtain contact features on individual level 
#Summarize individual contact patterns 
ind_daily_contacts <- hts_cnet_filt_sym %>% 
  mutate(index_contact = if_else(grepl("-001", indid2), 1, 0)) %>% 
  group_by(indid1,index_contact, date) %>% 
  summarize(ind_contact_freq_p_day = n(),
            ind_contact_dur_p_day = sum(duration_sec),
            ind_contact_freq_full = n(),
            ind_ave_dur_p_day = ind_contact_dur_p_day/ind_contact_freq_p_day,
            ind_contatc_freq_median_p_day = median(ind_contact_freq_p_day)) %>% 
  ungroup() %>% 
  complete(indid1,index_contact,
           fill = list(ind_contact_freq_p_day = 0, ind_contact_dur_p_day = 0, 
                       ind_contact_freq_full = 0, 
                       ind_ave_dur_p_day = 0, ind_contatc_freq_median_p_day = 0)) %>% 
  filter(index_contact==1 & grepl("-001", indid1)==FALSE)

#Merge into daily summary
ind_daily_contacts <- left_join(ind_daily_contacts, ind_days_nodes)

#Generate summary variables
ind_contact_summary <- ind_daily_contacts %>% 
  group_by(indid1) %>% 
  summarize(ind_contact_days_full = sum(ind_contact_days_full),
            ind_contact_freq_full = sum(ind_contact_freq_full),
            ind_contact_dur_full = sum(ind_contact_dur_p_day),
            ind_contact_ave_dur_full = if_else(ind_contact_dur_full!=0, ind_contact_dur_full/ind_contact_freq_full, 0),
            ind_contatc_freq_mean_p_day = mean(ind_contact_freq_p_day),
            ind_contatc_freq_median_p_day = median(ind_contact_freq_p_day),
            ind_contatc_dur_mean_p_day = mean(ind_contact_dur_p_day),
            ind_contatc_dur_median_p_day = median(ind_contact_dur_p_day),
            ind_ave_dur_mean_p_day = mean(ind_ave_dur_p_day),
            ind_ave_dur_median_p_day = median(ind_ave_dur_p_day))

rm(ind_daily_contacts, ind_days_nodes)

#Grouped analysis
#Contact parameters with pos individuals (grp_contact_summary, generated from hts_cnet_filt_sym)
grp_daily_contacts <- hts_cnet_filt_sym %>% 
  group_by(hh, indid1, sars_indid2, pair_sars, hcir, date) %>%
  summarize(grp_contact_freq_p_day = n(),
            grp_contact_dur_p_day = sum(duration_sec),
            grp_contact_freq_full = n(),
            grp_ave_dur_p_day = grp_contact_dur_p_day/grp_contact_freq_p_day,
            grp_contatc_freq_median_p_day = median(grp_contact_freq_p_day)) %>% 
  ungroup() %>% 
  complete(indid1, sars_indid2, 
           fill = list(grp_contact_freq_p_day = 0, grp_contact_dur_p_day = 0, 
                       grp_contact_freq_full = 0, 
                       grp_ave_dur_p_day = 0, grp_contatc_freq_median_p_day = 0)) %>% 
  filter(sars_indid2=="Positive")

#Obtain number of days and nodes
grp_days_nodes <- hts_cnet_filt_sym %>% 
  group_by(indid1) %>% 
  summarise(grp_contact_days_full = n_distinct(date),
            nodes_in_contact_d_grp = n_distinct(indid2))
#Merge into daily summary
grp_daily_contacts <- left_join(grp_daily_contacts, grp_days_nodes)
#Generate summary variables
grp_contact_summary <- grp_daily_contacts %>% 
  group_by(indid1) %>% 
  summarize(grp_contact_days_full = sum(grp_contact_days_full),
            grp_contact_freq_full = sum(grp_contact_freq_full),
            grp_contact_dur_full = sum(grp_contact_dur_p_day),
            grp_contact_ave_dur_full = if_else(grp_contact_dur_full!=0, grp_contact_dur_full/grp_contact_freq_full, 0),
            grp_contatc_freq_mean_p_day = mean(grp_contact_freq_p_day),
            grp_contatc_freq_median_p_day = median(grp_contact_freq_p_day),
            grp_contatc_dur_mean_p_day = mean(grp_contact_dur_p_day),
            grp_contatc_dur_median_p_day = median(grp_contact_dur_p_day),
            grp_ave_dur_mean_p_day = mean(grp_ave_dur_p_day),
            grp_ave_dur_median_p_day = median(grp_ave_dur_p_day),
            nodes_in_contact_d_grp = max(nodes_in_contact_d_grp, na.rm = TRUE))
rm(grp_daily_contacts, grp_days_nodes)

#Prepare for analysis

#Combine with metadata
#Excluding index
ind_contact_summary<-left_join(ind_contact_summary, meta, by=c("indid1"="indid"))

#Including index
ind_contact_summary_inc_index <-full_join(ind_contact_summary_inc_index, meta, by=c("indid1"="indid"))

#Grouped analysis
grp_contact_summary <-left_join(grp_contact_summary, meta, by=c("indid1"="indid"))

#---- Analysis 

#Summary of overall contact parameters
cont_overall_sum <- ind_contact_summary_inc_index %>% 
rename("Site" = site,
       "Index" = index,
       "Age (years)" = agegrp9,
       "Sex" = sex,
       "SARS-CoV-2 infection" = sars,
       "Median daily contact duration" = ind_contatc_dur_median_p_day,
       "Median daily contact frequency" = ind_contatc_freq_median_p_day,
       "Median daily average contact duration" = ind_ave_dur_median_p_day)

tableTwo <- CreateTableOne(vars = c("Median daily contact duration", "Median daily contact frequency", "Median daily average contact duration"),
                           strata = c("Age (years)"), 
                           data = cont_overall_sum,
                           addOverall = TRUE)

print(tableTwo, addOverall = TRUE, 
      nonnormal  = c("Median daily contact duration", "Median daily contact frequency", "Median daily average contact duration"),
      quote = TRUE, noSpaces = TRUE, contDigits=0)

#By site

cont_overall_sum_k <- cont_overall_sum %>% 
    filter(Site == "Klerksdorp") 

tableTwo <- CreateTableOne(vars = c("Median daily contact duration", "Median daily contact frequency", "Median daily average contact duration"),
                           strata = c("Age (years)"), 
                           data = cont_overall_sum_k,
                           addOverall = TRUE)

print(tableTwo, addOverall = TRUE,
      nonnormal  = c("Median daily contact duration", "Median daily contact frequency", "Median daily average contact duration"),
      quote = TRUE, noSpaces = TRUE, contDigits=0)

cont_overall_sum_s <- cont_overall_sum %>% 
         filter(Site == "Soweto") 

tableTwo <- CreateTableOne(vars = c("Median daily contact duration", "Median daily contact frequency", "Median daily average contact duration"),
                           strata = c("Age (years)"), 
                           data = cont_overall_sum_s,
                           addOverall = TRUE)

print(tableTwo, addOverall = TRUE, 
      nonnormal  = c("Median daily contact duration", "Median daily contact frequency", "Median daily average contact duration"),
      quote = TRUE, noSpaces = TRUE, contDigits=0)
rm(tableTwo)

#As figure

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
#Frequency
figfreq <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=`Median daily contact frequency`, fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Median daily \nclose-proximity event frequency \n(log scale)\n") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Duration
figdur <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=(`Median daily contact duration`/60), fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Median daily cumulative \nclose-proximity events \n(minutes, log scale)") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1,
                    guide = FALSE) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 8))

#Average duration
figavedur <- ggplot(cont_overall_sum, aes(x=`Age (years)`, y=`Median daily average contact duration`, fill=`Age (years)`)) +
  geom_boxplot(notch=FALSE) +
  facet_wrap(~Site) +
  labs(y= "Median daily \naverage close-proximity event \nduration (seconds, log scale)\n") +
  theme_classic()+
  scale_fill_manual(values = col_scheme1) +
  theme(legend.position="bottom",
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8)) +
  scale_y_continuous(trans='log10') +
  guides(fill = guide_legend(nrow = 1))


ggarrange(figdur, figfreq, figavedur, 
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3,
          heights = c(1,1,1.45))
#ggsave("Contact Parameters Overall.tiff", width = 16, height = 16, units = "cm", dpi=300)

rm(cont_overall_sum, cont_overall_sum_oa, cont_overall_sum_k, cont_overall_sum_s, figavedur, figdur, figfreq)

#Age-based contact matrices

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
   "pop_matrix_k" , "pop_matrix_s", freq_matrix, freq_matrix_klerk, freq_matrix_sow,
   )

#Additional figures
#Distribution of contact frequency, duration, average duration

#Individual
#Duration
fig_dur_ind <- ggplot(ind_contact_summary, aes(ind_contatc_dur_median_p_day, fill=as.character(sars))) + 
  geom_histogram(bins = 20) +
  labs(x = "Close-proximiy event \nduration (seconds, log scale)", y = "Number of \nindividuals") +
  scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c( "#1FA088FF", "#440154FF"),
                    guide = "none") +
  theme_classic()  +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  facet_wrap(~site)
#Frequency
fig_freq_ind  <- ggplot(ind_contact_summary, aes(ind_contatc_freq_median_p_day, fill=as.character(sars))) + 
  geom_histogram(bins = 20) +
  labs(x = "Close-proximiy event \nfrequency (log scale)", y = "Number of \nindividuals") +
  scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c( "#1FA088FF", "#440154FF"),
                    guide = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  facet_wrap(~site)
#Ave duration
fig_avedur_ind  <-ggplot(ind_contact_summary, aes(ind_ave_dur_median_p_day, fill=as.character(sars))) + 
  geom_histogram(bins = 20) +
  labs(x = "Average close-proximiy event \n duration (seconds, log scale)", y = "Number of \nindividuals") + 
  guides(fill=guide_legend(title="SARS-CoV-2")) +
  scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c( "#1FA088FF", "#440154FF"), labels = c("Negative", "Positive")) +
  theme_classic() +
  theme(legend.position="bottom",
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  facet_wrap(~site)

#Grouped

#Duration
fig_dur_grp <- ggplot(grp_contact_summary, aes(grp_contatc_dur_median_p_day, fill=as.character(sars))) + 
  geom_histogram(bins = 20) +
  labs(x = "Close-proximiy event \nduration (seconds, log scale)", y = "Number of \nindividuals") +
  scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c( "#1FA088FF", "#440154FF"),
                    guide = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  facet_wrap(~site)
#Frequency
fig_freq_grp  <-ggplot(grp_contact_summary, aes(grp_contatc_freq_median_p_day, fill=as.character(sars))) + 
  geom_histogram(bins = 20) +
  labs(x = "Close-proximiy event \nfrequency (log scale)", y = "Number of \nindividuals") +
  scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c( "#1FA088FF", "#440154FF"),
                    guide = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  facet_wrap(~site)
#Ave duration
fig_avedur_grp <-ggplot(grp_contact_summary, aes(grp_ave_dur_median_p_day, fill=as.character(sars))) + 
  geom_histogram(bins = 20) +
  labs(x = "Average close-proximiy event \n duration (seconds, log scale)", y = "Number of \nindividuals") + 
  guides(fill=guide_legend(title="SARS-CoV-2")) +
  scale_x_continuous(trans='log10') +
  scale_fill_manual(values = c( "#1FA088FF", "#440154FF"), labels = c("Negative", "Positive")) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  facet_wrap(~site)

ggarrange(fig_dur_ind, fig_dur_grp, fig_freq_ind, fig_freq_grp, fig_avedur_ind,   fig_avedur_grp, 
          labels = c("A", "D", "B", "E", "C", "F"),
          ncol = 2, nrow = 3,
          heights = c(1,1,1.45))
#ggsave("Contact Distribution Overall.tiff", width = 16, height = 16, units = "cm", dpi=600)
rm(fig_dur_ind, fig_dur_grp, fig_freq_ind, fig_freq_grp, fig_avedur_ind,   fig_avedur_grp)

ind_contact_summary$sars <- factor(ind_contact_summary$sars)
#Build model: Individual-level 
#Due to privacy concerns not all meta data included. Only final model shown
ind_fin_mod <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                       agegrp9 + sex + 
                       (1 | site) + (1 | hhid), 
                     data = ind_contact_summary,
                     family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median duration of contacts per day with index (ind_contatc_dur_median_p_day)
ind_fin_mod_mdur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contatc_dur_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Median frequency of contacts per day with index (ind_contatc_freq_median_p_day)
ind_fin_mod_mfeq <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_contatc_freq_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(ind_fin_mod_mfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Average duration of contacts overall with index (ind_contact_ave_dur_full)
ind_fin_mod_adur <- glmer(sars ~ ixagegrp9 + ixminctcat2 + ixesarsvarf1 + 
                            agegrp9 + sex + 
                            ind_ave_dur_median_p_day + 
                            (1 | site) + (1 | hhid), 
                          data = ind_contact_summary,
                          family = binomial("logit"))
#Model does not converge!
print(data.frame(tidy(ind_fin_mod_adur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote = TRUE, noSpaces = TRUE)

#Build model: Grouped analysis ----
grp_contact_summary$sars <- factor(grp_contact_summary$sars)
grp_fin_mod <- glmer(sars ~ agegrp9 + smokecignow1 + bmicat + ixesarsvarf1 +
                       (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                     data = grp_contact_summary,
                     family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote=TRUE, noSpaces = TRUE)

#Median duration of contacts per day with infected (grp_contatc_dur_median_p_day)
grp_fin_mod_mdur <- glmer(sars ~ agegrp9 + smokecignow1 + bmicat + ixesarsvarf1 +
                            grp_contatc_dur_median_p_day + 
                            (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                          data = grp_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_mdur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote=TRUE, noSpaces = TRUE)

#Median frequency of contacts per day with infected (grp_contatc_freq_median_p_day)
grp_fin_mod_mfeq <- glmer(sars ~ agegrp9 + smokecignow1 + bmicat + ixesarsvarf1 +
                            grp_contatc_freq_median_p_day + 
                            (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                          data = grp_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_mfeq,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote=TRUE, noSpaces = TRUE)

#Average duration of contacts overall with infected (grp_contact_ave_dur_full)
grp_fin_mod_adur <- glmer(sars ~ agegrp9 + smokecignow1 + bmicat + ixesarsvarf1 +
                            grp_contact_ave_dur_full + 
                            (1 | site) + (1 | hhid) + offset(nodes_in_contact_d_grp), 
                          data = grp_contact_summary,
                          family = binomial("logit"))
print(data.frame(tidy(grp_fin_mod_adur,conf.int=TRUE,exponentiate=TRUE,effects="fixed")) %>% 
        mutate_if(is.numeric, round, 4), quote=TRUE, noSpaces = TRUE, scientific = FALSE)

#Also assess association with Wilcoxon rank-sum test
#Individual analysis
#Median duration
ind_mdur_wilcox <- wilcox.test(ind_contatc_dur_median_p_day ~ sars, data = ind_contact_summary,
                               paired = FALSE, 
                               exact = FALSE)
#Median frequency
ind_mfreq_wilcox <- wilcox.test(ind_contatc_freq_median_p_day ~ sars, data = ind_contact_summary,
                                paired = FALSE, 
                                exact = FALSE)
#Average duration
ind_adur_wilcox <- wilcox.test(ind_contact_ave_dur_full ~ sars, data = ind_contact_summary,
                               paired = FALSE, 
                               exact = FALSE)

#Group analysis
#Median duration
grp_mdur_wilcox <-wilcox.test(grp_contatc_dur_median_p_day ~ sars, data = grp_contact_summary,
                              paired = FALSE, 
                              exact = FALSE)
#Median frequency
grp_mfreq_wilcox <- wilcox.test(grp_contatc_freq_median_p_day ~ sars, data = grp_contact_summary,
                                paired = FALSE, 
                                exact = FALSE)
#Average duration
grp_adur_wilcox <- wilcox.test(grp_contact_ave_dur_full ~ sars, data = grp_contact_summary,
                               paired = FALSE, 
                               exact = FALSE)

wilcox_results <- cbind(ind_mdur_wilcox, ind_mfreq_wilcox, ind_adur_wilcox,
                        grp_mdur_wilcox, grp_mfreq_wilcox, grp_adur_wilcox)
wilcox_results
rm(ind_mdur_wilcox, ind_mfreq_wilcox, ind_adur_wilcox,
   grp_mdur_wilcox, grp_mfreq_wilcox, grp_adur_wilcox)
