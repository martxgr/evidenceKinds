# libraries
library(tidyverse)
library(lme4)
library(psych)
library(corrplot)
library(mclust)
library(patchwork)
library(GGally)
library(ade4)
library(lavaan)
#for the notation
options(scipen=999)

d <- read.csv("rb.p2.csv") 

#### Clean Data ####

DVs <- c("jus_1", "mot_1",
         "really_1", "consider_1", "predict_1", "context_1", 
         "surprise_1", "bet_1", "double_1")

dimensions <- c("jus_1", "mot_1",
         "really_1", "consider_1", "predict_1", "context_1", 
         "surprise_1", "bet_1", "double_1")


  mutate(domain = case_when(superstition != "" ~ "superstition",
                            religious != "" ~ "religious",
                            superlative != "" ~ "superlative",
                            political != "" ~ "political")) %>% 
  filter(age != "") %>% 
  mutate(across(all_of(DVs), as.numeric)) %>% 
  rename_all(~ gsub("_1$", "", .)) %>% 
  dplyr::select(prolific, text, 
                jus, mot,
                really, consider, predict, context,
                surprise, bet, double, 
                age, gen, race, domain)

#reverse score items so that high scores are more evidency/really believe
d$double <- -1*(d$double - 100)
d$context <- -1*(d$context - 100)

# prepare for clustering
d.c <- d %>% 
  dplyr::select(really, consider, predict, context,
                surprise, bet, mot, jus, double, domain)