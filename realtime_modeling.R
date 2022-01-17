# vim:foldmethod=marker
# model training and prediction for realtime-analysis

### SETUP {{{
rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(randomForest))


# dates formatted as "YYYY-MM-DD"
args <- commandArgs(trailingOnly = TRUE)
date <- as.character(args[1]) # the date on the sequence header in format YYYY-MM-DD
new_patient <- as.numeric(args[2]) # the patient-ID that's used in the prediction
use_n <- as.numeric(args[3]) # number of past entries to use as training data for the models
lab_interv <- as.numeric(args[4]) # number of days that are allowed between sequence timestamp and sample_day of lab resuls

# for testing only:
#date <- "2006-01-03"
#new_patient <- 19056
#use_n <- 5000
#lab_interv <- 30

data_dir <- "/home/cluster/mlabar/data/realtime"
preds_table_path <- paste(data_dir, "realtime_predictions.csv", sep = "/")
cc_attr_path <- paste(data_dir, "cc_attributes.csv", sep = "/")
cc_bl_path <- paste(data_dir, "cc_bl_attributes.csv", sep = "/")

patient_dir <- paste(paste(data_dir, date, sep = "/"), new_patient, sep = "_")
models_path <- paste(patient_dir, "models.RData", sep = "/")
pred_data_path <- paste(patient_dir, "predict.csv", sep = "/")
train_data_all_path <- paste(patient_dir, "training_all.csv", sep = "/")
train_data_clustered_path <- paste(patient_dir, "training_clustered.csv", sep = "/")
train_data_singletons_path <- paste(patient_dir, "training_singletons.csv", sep = "/")
train_data_all_lab_path <- paste(patient_dir, "training_all_lab.csv", sep = "/")
train_data_clustered_lab_path <- paste(patient_dir, "training_clustered_lab.csv", sep = "/")
train_data_singletons_lab_path <- paste(patient_dir, "training_singletons_lab.csv", sep = "/")
#}}}

### DATA IMPORT {{{
attr <- read.csv(cc_attr_path)
attr_bl <- read.csv(cc_bl_path)
attr_bl <- attr_bl[order(attr_bl$sequence_timestamp, decreasing = T), ]

attr_bl <- attr_bl[1:use_n, ] # choose the use_n most recent entries
attr_bl <- filter(attr_bl, !is.na(nd_gr_bool)) # remove entries with no outcome recorded
attr_bl <- filter(attr_bl, pat_ethnicity != 9)

# variable recoding
# attr
{
    attr$pat_risk[attr$pat_risk > 4] <- 0
    attr$pat_risk[attr$pat_risk == 4] <- 3
    
    attr$pat_sex <- as.factor(attr$pat_sex)
    attr$pat_center1 <- as.factor(attr$pat_center1)
    attr$pat_ethnicity <- as.factor(attr$pat_ethnicity)
    attr$pat_education <- as.factor(attr$pat_education)
    attr$pat_profession <- as.factor(attr$pat_profession)
    attr$pat_risk <- as.factor(attr$pat_risk)
    
    # dates:
    attr <- attr %>% mutate(lab_labdate = as.Date(lab_labdate, "%F", optional = T))
    attr <- attr %>% mutate(fup_fupdate = as.Date(fup_fupdate, "%F", optional = T))
    attr <- attr %>% mutate(adhe_ad_date = as.Date(adhe_ad_date, "%F", optional = T))
    attr <- attr %>% mutate(sequence_timestamp = as.Date(sequence_timestamp, "%F", optional = T))
    
    seq_years <- sapply(sapply(attr$sequence_timestamp, substr, start = 0, stop = 4), as.integer)
    attr$pat_age <- seq_years - attr$pat_born
    
    # age:
    attr <- mutate(attr, pat_age_cat = case_when(pat_age < 15 ~ "0-14", 
                                                       pat_age >= 15 & pat_age < 30 ~ "15-29", 
                                                       pat_age >= 30 & pat_age < 45 ~ "30-44", 
                                                       pat_age >= 45 & pat_age < 60 ~ "45-59", 
                                                       pat_age >= 60 ~ "60+"))
    
    # rna:
    attr <- mutate(attr, lab_rna_cat = case_when(lab_rna == 0 ~ "0", 
                                                       lab_rna >= 1 & lab_rna < 100 ~ "1-99", 
                                                       lab_rna >= 100 & lab_rna < 1000 ~ "100-999", 
                                                       lab_rna >= 1000 ~ "1000+"))
    
    # cd4:
    attr <- mutate(attr, lab_cd4_cat = case_when(lab_cd4 < 200 ~ "0-199", 
                                                       lab_cd4 >= 200 & lab_cd4 < 350 ~ "200-349", 
                                                       lab_cd4 >= 350 & lab_cd4 < 500 ~ "350-499", 
                                                       lab_cd4 >= 500 ~ "500+"))
    
    # alcohol:
    attr <- mutate(attr, fup_alc_comb = case_when(fup_alcohol == 0 | fup_alc_freq == 1 ~ 0,  fup_alcohol == 1 | fup_alc_freq %in% c(2, 3, 4, 5) ~ 1, fup_alcohol == 9 | fup_alc_freq == 9 ~ 9))
    
    # IV-drugs:
    attr <- mutate(attr, fup_iv_comb = case_when(fup_coca_iv == 0 & fup_hero_iv == 0 & fup_other_iv == 0 ~ 0, fup_coca_iv == 1 | fup_hero_iv == 1 | fup_other_iv == 1 ~ 1))
    
    # Non-IV-drugs:
    attr <- mutate(attr, fup_ni_comb = case_when(fup_cana_ni == 0 & fup_coca_ni == 0 & fup_hero_ni == 0 & fup_other_ni == 0 ~ 0, fup_cana_ni == 1 | fup_coca_ni == 1 | fup_hero_ni == 1 | fup_other_ni == 1 ~ 1))
    
    attr <- mutate(attr, clustered = if_else(cl_size > 1, T, F))
    attr <- attr[ , order(names(attr))]
    front_vars <- c("id", "patient_id", "sequence_timestamp", "clustered") # variables to be displayed before all others
    attr <- select(attr, all_of(c(front_vars, names(attr)[!names(attr) %in% front_vars])))
}

# attr_bl
{
    attr_bl$pat_risk[attr_bl$pat_risk > 4] <- 0
    attr_bl$pat_risk[attr_bl$pat_risk == 4] <- 3
    
    attr_bl$pat_sex <- as.factor(attr_bl$pat_sex)
    attr_bl$pat_center1 <- as.factor(attr_bl$pat_center1)
    attr_bl$pat_ethnicity <- as.factor(attr_bl$pat_ethnicity)
    attr_bl$pat_education <- as.factor(attr_bl$pat_education)
    attr_bl$pat_profession <- as.factor(attr_bl$pat_profession)
    attr_bl$pat_risk <- as.factor(attr_bl$pat_risk)
    
    # outcome from 0/1 to False/True: IMPORTANT!!!
    attr_bl <- mutate(attr_bl, nd_gr_bool = nd_gr_bool == 1)
    
    # dates:
    attr_bl <- attr_bl %>% mutate(lab_labdate = as.Date(lab_labdate, "%F", optional = T))
    attr_bl <- attr_bl %>% mutate(fup_fupdate = as.Date(fup_fupdate, "%F", optional = T))
    attr_bl <- attr_bl %>% mutate(adhe_ad_date = as.Date(adhe_ad_date, "%F", optional = T))
    attr_bl <- attr_bl %>% mutate(sequence_timestamp = as.Date(sequence_timestamp, "%F", optional = T))
    
    seq_years <- sapply(sapply(attr_bl$sequence_timestamp, substr, start = 0, stop = 4), as.integer)
    attr_bl$pat_age <- seq_years - attr_bl$pat_born
    
    # age:
    attr_bl <- mutate(attr_bl, pat_age_cat = case_when(pat_age < 15 ~ "0-14", 
                                                       pat_age >= 15 & pat_age < 30 ~ "15-29", 
                                                       pat_age >= 30 & pat_age < 45 ~ "30-44", 
                                                       pat_age >= 45 & pat_age < 60 ~ "45-59", 
                                                       pat_age >= 60 ~ "60+"))
    
    # rna:
    attr_bl <- mutate(attr_bl, lab_rna_cat = case_when(lab_rna == 0 ~ "0", 
                                                       lab_rna >= 1 & lab_rna < 100 ~ "1-99", 
                                                       lab_rna >= 100 & lab_rna < 1000 ~ "100-999", 
                                                       lab_rna >= 1000 ~ "1000+"))
    
    # cd4:
    attr_bl <- mutate(attr_bl, lab_cd4_cat = case_when(lab_cd4 < 200 ~ "0-199", 
                                                       lab_cd4 >= 200 & lab_cd4 < 350 ~ "200-349", 
                                                       lab_cd4 >= 350 & lab_cd4 < 500 ~ "350-499", 
                                                       lab_cd4 >= 500 ~ "500+"))
    
    # alcohol:
    attr_bl <- mutate(attr_bl, fup_alc_comb = case_when(fup_alcohol == 0 | fup_alc_freq == 1 ~ 0,  fup_alcohol == 1 | fup_alc_freq %in% c(2, 3, 4, 5) ~ 1, fup_alcohol == 9 | fup_alc_freq == 9 ~ 9))
    
    # IV-drugs:
    attr_bl <- mutate(attr_bl, fup_iv_comb = case_when(fup_coca_iv == 0 & fup_hero_iv == 0 & fup_other_iv == 0 ~ 0, fup_coca_iv == 1 | fup_hero_iv == 1 | fup_other_iv == 1 ~ 1))
    
    # Non-IV-drugs:
    attr_bl <- mutate(attr_bl, fup_ni_comb = case_when(fup_cana_ni == 0 & fup_coca_ni == 0 & fup_hero_ni == 0 & fup_other_ni == 0 ~ 0, fup_cana_ni == 1 | fup_coca_ni == 1 | fup_hero_ni == 1 | fup_other_ni == 1 ~ 1))
    
    attr_bl <- mutate(attr_bl, clustered = if_else(cl_size > 1, T, F))
    attr_bl <- attr_bl[ , order(names(attr_bl))]
    front_vars <- c("id", "patient_id", "sequence_timestamp", "clustered") # variables to be displayed before all others
    attr_bl <- select(attr_bl, all_of(c(front_vars, names(attr_bl)[!names(attr_bl) %in% front_vars])))
}

attr_bl_clustered <- filter(attr_bl, nd_degree > 0) # all clustered sequences
attr_bl_singletons <- filter(attr_bl, nd_degree == 0) # all singletons

# remove entries with sample_day too far from sequence_timestamp for models that deped on lab data:
attr_bl_all_lab <- filter(attr_bl, as.vector(abs(as.Date(attr_bl$sequence_timestamp, "%F", optional = T) - as.Date(attr_bl$lab_sample_day, "%F", optional = T))) <= lab_interv)
attr_bl_clustered_lab <- filter(attr_bl_clustered, as.vector(abs(as.Date(attr_bl_clustered$sequence_timestamp, "%F", optional = T) - as.Date(attr_bl_clustered$lab_sample_day, "%F", optional = T))) <= lab_interv)
attr_bl_singletons_lab <- filter(attr_bl_singletons, as.vector(abs(as.Date(attr_bl_singletons$sequence_timestamp, "%F", optional = T) - as.Date(attr_bl_singletons$lab_sample_day, "%F", optional = T))) <= lab_interv)

# preparing a directory and saving modeling data
dir.create(patient_dir)
write.csv(attr, pred_data_path, row.names = F)
write.csv(attr_bl, train_data_all_path, row.names = F)
write.csv(attr_bl_clustered, train_data_clustered_path, row.names = F)
write.csv(attr_bl_singletons, train_data_singletons_path, row.names = F)
write.csv(attr_bl_all_lab, train_data_all_lab_path, row.names = F)
write.csv(attr_bl_clustered_lab, train_data_clustered_lab_path, row.names = F)
write.csv(attr_bl_singletons_lab, train_data_singletons_lab_path, row.names = F)

### MODEL FITTING {{{
# model formulae:
frml_mix <-             formula(as.factor(nd_gr_bool) ~ nd_degree + cl_size + cl_past_gr + pat_risk + pat_center1 + pat_age_cat + pat_sex + lab_rna_cat)
frml_mix_cl <-          formula(as.factor(nd_gr_bool) ~ nd_degree + cl_size + cl_past_gr)
frml_mix_pat <-         formula(as.factor(nd_gr_bool) ~ pat_risk + pat_center1 + pat_age_cat + pat_sex + lab_rna_cat)
frml_cl <-              formula(as.factor(nd_gr_bool) ~ nd_degree + cl_size + cl_closeness_median + nd_closeness + cl_density + cl_distance_median + cl_past_gr)
frml_pat <-             formula(as.factor(nd_gr_bool) ~ pat_risk + pat_center1 + pat_age_cat + pat_sex + lab_rna_cat)
frml_vsurf_thres <-     formula(as.factor(nd_gr_bool) ~ cl_closeness_median + nd_closeness + cl_density + nd_degree + cl_distance_median + nd_betweenness + cl_size + cl_past_gr + pat_center1 + pat_age_cat + pat_ethnicity + pat_sex + lab_cd4_cat + lab_rna_cat + cl_betweenness_median + cl_degree_median + pat_x_pref)
frml_vsurf_interp <-    formula(as.factor(nd_gr_bool) ~ cl_closeness_median + nd_closeness + nd_degree + cl_density + cl_distance_median + cl_past_gr + cl_size + nd_betweenness + cl_degree_median)
frml_vsurf_pred <-      formula(as.factor(nd_gr_bool) ~ nd_closeness + cl_closeness_median + nd_degree + cl_size)

# logistic regressions:
mix_all_logr <-               glm(frml_mix, data = attr_bl_all_lab, family = binomial, maxit = 100)
mix_clustered_logr <-         glm(frml_mix, data = attr_bl_clustered_lab, family = binomial, maxit = 100)
mix_cl_logr <-                glm(frml_mix_cl, data = attr_bl_clustered, family = binomial, maxit = 100)
mix_pat_logr <-               glm(frml_mix_pat, data = attr_bl_all_lab, family = binomial, maxit = 100)
cl_logr <-                    glm(frml_cl, data = attr_bl_clustered, family = binomial, maxit = 100)
pat_all_logr <-               glm(frml_pat, data = attr_bl_all_lab, family = binomial, maxit = 100)
pat_singletons_logr <-        glm(frml_pat, data = attr_bl_singletons_lab, family = binomial, maxit = 100)
vsurf_thres_all_logr <-       glm(frml_vsurf_thres, data = attr_bl_all_lab, family = binomial, maxit = 100)
vsurf_thres_clustered_logr <- glm(frml_vsurf_thres, data = attr_bl_clustered_lab, family = binomial, maxit = 100)
vsurf_interp_logr <-          glm(frml_vsurf_interp, data = attr_bl_clustered, family = binomial, maxit = 100)
vsurf_pred_logr <-            glm(frml_vsurf_pred, data = attr_bl_clustered, family = binomial, maxit = 100)

# random forests:
wt_true <-                  sum(attr_bl[["nd_gr_bool"]]) / length(attr_bl[["nd_gr_bool"]])
mix_all_rf <-               randomForest(frml_mix, data = attr_bl_all_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
mix_clustered_rf <-         randomForest(frml_mix, data = attr_bl_clustered_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
mix_cl_rf <-                randomForest(frml_mix_cl, data = attr_bl_clustered, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
mix_pat_rf <-               randomForest(frml_mix_pat, data = attr_bl_all_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
cl_rf <-                    randomForest(frml_cl, data = attr_bl_clustered, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
pat_all_rf <-               randomForest(frml_pat, data = attr_bl_all_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
pat_singletons_rf <-        randomForest(frml_pat, data = attr_bl_singletons_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
vsurf_thres_all_rf <-       randomForest(frml_vsurf_thres, data = attr_bl_all_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
vsurf_thres_clustered_rf <- randomForest(frml_vsurf_thres, data = attr_bl_clustered_lab, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
vsurf_interp_rf <-          randomForest(frml_vsurf_interp, data = attr_bl_clustered, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)
vsurf_pred_rf <-            randomForest(frml_vsurf_pred, data = attr_bl_clustered, classwt = c(1 - wt_true, wt_true), importance = T, na.action = na.omit)

# saving the models
save(mix_all_logr, mix_clustered_logr, mix_cl_logr, mix_pat_logr, cl_logr, pat_all_logr, pat_singletons_logr, 
    vsurf_thres_all_logr, vsurf_thres_clustered_logr, vsurf_interp_logr, vsurf_pred_logr,
    mix_all_rf, mix_clustered_rf, mix_cl_rf, mix_pat_rf, cl_rf, pat_all_rf, pat_singletons_rf, 
    vsurf_thres_all_rf, vsurf_thres_clustered_rf, vsurf_interp_rf, vsurf_pred_rf, 
    file = models_path)
#}}}

### PREDICTION {{{
# a dirty hack to even out the factor levels and prevent errors in the prediction
attr <- rbind(select(attr, -any_of(c("nd_gr", "nd_gr_bool", "cl_gr"))), select(attr_bl, -any_of(c("nd_gr", "nd_gr_bool", "cl_gr"))))
attr <- attr[1, ] # keep only the top row containing the new patient's data
row.names(attr) <- 1

pr_mix_all_logr <-                predict(mix_all_logr, type = "response", newdata = attr)
pr_mix_clustered_logr <-          predict(mix_clustered_logr, type = "response", newdata = attr)
pr_mix_cl_logr <-                 predict(mix_cl_logr, type = "response", newdata = attr)
pr_mix_pat_logr <-                predict(mix_pat_logr, type = "response", newdata = attr)
pr_cl_logr <-                     predict(cl_logr, type = "response", newdata = attr)
pr_pat_all_logr <-                predict(pat_all_logr, type = "response", newdata = attr)
pr_pat_singletons_logr <-         predict(pat_singletons_logr, type = "response", newdata = attr)
pr_vsurf_thres_all_logr <-        predict(vsurf_thres_all_logr, type = "response", newdata = attr)
pr_vsurf_thres_clustered_logr <-  predict(vsurf_thres_clustered_logr, type = "response", newdata = attr)
pr_vsurf_interp_logr <-           predict(vsurf_interp_logr, type = "response", newdata = attr)
pr_vsurf_pred_logr <-             predict(vsurf_pred_logr, type = "response", newdata = attr)

pr_mix_all_rf <-                  as.data.frame(predict(mix_all_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_mix_clustered_rf <-            as.data.frame(predict(mix_clustered_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_mix_cl_rf <-                   as.data.frame(predict(mix_cl_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_mix_pat_rf <-                  as.data.frame(predict(mix_pat_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_cl_rf <-                       as.data.frame(predict(cl_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_pat_all_rf <-                  as.data.frame(predict(pat_all_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_pat_singletons_rf <-           as.data.frame(predict(pat_singletons_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_vsurf_thres_all_rf <-          as.data.frame(predict(vsurf_thres_all_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_vsurf_thres_clustered_rf <-    as.data.frame(predict(vsurf_thres_clustered_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_vsurf_interp_rf <-             as.data.frame(predict(vsurf_interp_rf, type = "prob", newdata = attr))[["TRUE"]]
pr_vsurf_pred_rf <-               as.data.frame(predict(vsurf_pred_rf, type = "prob", newdata = attr))[["TRUE"]]

# make a new row for the prediction
new_pred_df <- data.frame("patient_id" =                   new_patient,
                "sequence_timestamp" =                     date,
                "clustered" =                              attr[["clustered"]],
                "prediction_timestamp" =                   as.character(Sys.time()),
                "prediction_mix_all_logr" =                pr_mix_all_logr,
                "prediction_mix_all_rf" =                  pr_mix_all_rf,
                "prediction_mix_clustered_logr" =          pr_mix_clustered_logr,
                "prediction_mix_clustered_rf" =            pr_mix_clustered_rf,
                "prediction_mix_cl_logr" =                 pr_mix_cl_logr,
                "prediction_mix_cl_rf" =                   pr_mix_cl_rf,
                "prediction_mix_pat_logr" =                pr_mix_pat_logr,
                "prediction_mix_pat_rf" =                  pr_mix_pat_rf,
                "prediction_cl_logr" =                     pr_cl_logr,
                "prediction_cl_rf" =                       pr_cl_rf,
                "prediction_pat_all_logr" =                pr_pat_all_logr,
                "prediction_pat_all_rf" =                  pr_pat_all_rf,
                "prediction_pat_singletons_logr" =         pr_pat_singletons_logr,
                "prediction_pat_singletons_rf" =           pr_pat_singletons_rf,
                "prediction_vsurf_thres_all_logr" =        pr_vsurf_thres_all_logr,
                "prediction_vsurf_thres_all_rf" =          pr_vsurf_thres_all_rf,
                "prediction_vsurf_thres_clustered_logr" =  pr_vsurf_thres_clustered_logr,
                "prediction_vsurf_thres_clustered_rf" =    pr_vsurf_thres_clustered_rf,
                "prediction_vsurf_interp_logr" =           pr_vsurf_interp_logr,
                "prediction_vsurf_interp_rf" =             pr_vsurf_interp_rf,
                "prediction_vsurf_pred_logr" =             pr_vsurf_pred_logr,
                "prediction_vsurf_pred_rf" =               pr_vsurf_pred_rf,
                "outcome_timestamp" =                      NA,
                "outcome" =                                NA)

print(new_pred_df)

# read prediction table and check old predictions for outcomes
if (file.exists(preds_table_path)) {
    pred_df <- read.csv(preds_table_path)
    NA_pids <- filter(pred_df, is.na(outcome))[["patient_id"]]
    
    for (pid in NA_pids) {
        if (pid %in% attr_bl[["patient_id"]]) {
            outc <- filter(attr_bl, patient_id == pid)[["nd_gr_bool"]]
            pred_df[pred_df["patient_id"] == pid, ]["outcome"] <- outc
            pred_df[pred_df["patient_id"] == pid, ]["outcome_timestamp"] <- date
        }
    }
} else {
    pred_df <- data.frame(matrix(nrow = 0, ncol = 28))
    names(pred_df) <- names(new_pred_df)
}

# add new prediction and write to file
pred_df <- rbind(pred_df, new_pred_df)
write.csv(pred_df, preds_table_path, row.names = F)
#}}}
