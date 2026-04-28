
library(dplyr)
library(ggplot2)
library(viridis)
library(effsize) # stats
library(tidyverse)

transform_data_evo_v_evo <- function(df){
  small <- df %>%
    filter(tag_perm != "none") %>%
    mutate(sym_tag = case_when(sym_tag=="" ~ NA, .default = sym_tag),
           tag_distance = as.numeric(tag_distance)) %>%
    filter(!is.na(sym_tag) & (sym_int >= 0.6 | sym_int < -0.6)) %>%
    mutate(cond = case_when(sym_int < -0.6 ~ "para", .default = "mut"),
           tag_perm = as.factor(tag_perm), param = as.factor(param))
  
  small <- small %>% group_by(tag_perm, param) %>% mutate(para_pres = "para"%in%cond)
  small <- small %>% group_by(tag_perm, param) %>% mutate(mut_pres = "mut"%in%cond)
  small <- small %>% filter(para_pres == T & mut_pres == T)
  small$cond <- as.factor(small$cond)
  return(small)
}

calculate_p_vals <- function(df, param_prefix="",csv_name="unnamed_csv.csv"){
  nm <- c()
  md <- c()
  n_test <- 0
  for (tp in unique(df$tag_perm)){
    for(p in unique(df$param)){
      a <- subset(df, tag_perm == tp & param == p)
      if("para" %in% a$cond){
        n_test <- n_test + 1
        nm <- append(nm, paste0("tp",tp,param_prefix,p))
        md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
        
        print("\n*******************")
        print(nm[n_test])
        print(md[n_test])
      }
    }
  }
  print(n_test)
  
  
  models <- data.frame(names = nm, p_vals=md)
  write_csv(models, csv_name)
}

calculate_p_vals_and_eff_sizes <- function(df, param_prefix="",csv_name="unnamed_csv.csv"){
  nm <- c()
  md <- c()
  effect <- c()
  n_test <- 0
  for (tp in unique(df$tag_perm)){
    for(p in unique(df$param)){
      a <- subset(df, tag_perm == tp & param == p)
      if("para" %in% a$cond){
        n_test <- n_test + 1
        nm <- append(nm, paste0("tp",tp,param_prefix,p))
        md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
        effect <- append(effect, cliff.delta(tag_distance ~ cond, data=a)$estimate)
        
        print("\n*******************")
        print(nm[n_test])
        print(md[n_test])
        print(effect[n_test])
      }
    }
  }
  print(n_test)
  
  
  models <- data.frame(names = nm, p_vals=md, deltas=effect)
  write_csv(models, csv_name)
}

## exp 1
print("doing VT... ")
filepath <- "../../data/exp_1_base_vt_evolving/"
dump_df <- read_csv(paste0(filepath,"org_dump.dat"),col_names=T)
dump_df <- dump_df %>% rename(param = vt)
head(dump_df)

print("loaded...")
small <- transform_data_evo_v_evo(dump_df)
calculate_p_vals(small, "_vt","vt_evolve_paravsmut_pvals.csv")
calculate_p_vals_and_eff_sizes(small, "_vt","vt_evolve_paravsmut_deltas.csv")
print("done with VT!")

## exp 2
print("doing TM... ")
filepath <- "../../data/exp_2_base_tagmut_evolving/"
dump_df <- read_csv(paste0(filepath,"org_dump.dat"),col_names=T)
dump_df <- dump_df %>% rename(param = tag_mut) %>% filter(param != 0)
head(dump_df)

print("loaded...")
small <- transform_data_evo_v_evo(dump_df)
calculate_p_vals(small, "_tm","tagmut_evolve_paravsmut_pvals.csv")
calculate_p_vals_and_eff_sizes(small, "_tm","tagmut_evolve_paravsmut_deltas.csv")
print("done with TM!")

