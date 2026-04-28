
library(dplyr)
library(ggplot2)
library(viridis)
library(effsize) # stats
library(tidyverse)

transform_data <- function(a, b){
  print("binding")
  a$sym_tag[a$sym_tag==""]<-NA
  a$cond<-"mut"# <- org_dump_df
  
  b$sym_tag[b$sym_tag==""]<-NA
  b$cond<-"para"
  
  small <- rbind(subset(b, !is.na(sym_tag)), subset(a, !is.na(sym_tag))) %>% filter(tag_perm != "none")
  rm(a)
  rm(b)
  
  print("done binding and deleting; about to transform")
  
  small <- small %>% group_by(tag_perm,param) %>% mutate(para_pres = "para"%in%cond)
  small <- small %>% group_by(tag_perm,param) %>% mutate(mut_pres = "mut"%in%cond)
  small <- small %>% filter(para_pres == T & mut_pres == T) %>% 
    mutate(tag_perm = as.factor(tag_perm), param = as.factor(param), cond = as.factor(cond))
  return(small)
}
calculate_p_vals_and_eff_sizes <- function(df, param_prefix="",csv_name="unnamed_csv.csv"){
  nm <- c()
  md <- c()
  effect <- c()
  n_test <- 0
  for (tp in unique(df$tag_perm)){
    for(v in unique(df$param)){
      a <- subset(df, tag_perm == tp & param == v)
      if("para" %in% a$cond){
        n_test <- n_test + 1
        print(paste0("tp",tp,param_prefix,v))
        nm <- append(nm, paste0("tp",tp,param_prefix,v))
        md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
        effect <- append(effect, cliff.delta(tag_distance ~ cond, data=a)$estimate)
  
        print("\n*******************")
        print(nm)
        print(md)
        print(effect)
	models <- data.frame(names = nm, p_vals=md, deltas=effect)

        write_csv(models, csv_name)
      }
    }
  }
  print(n_test)
  
  
  models <- data.frame(names = nm, p_vals=md, deltas=effect)
  
  write_csv(models, csv_name)
}


## exp 1

print("** VT ** ")
print("about to read mut")

filepath<- "../../data/exp_3_base_vt_mut_fixed/"
mut <- read_csv(paste0(filepath,"org_dump.dat")) %>% rename(param = vt)

print("done with mut; about to read para")
filepath<- "../../data/exp_3_base_vt_para_fixed/"
para <- read_csv(paste0(filepath,"org_dump.dat")) %>% rename(param = vt)

combo <- transform_data(mut, para)
calculate_p_vals_and_eff_sizes(combo, "_vt", "vt_fix_deltas.csv")



## exp 2
print("** Tag mutation ** ")
print("about to read mut")

filepath<- "../../data/exp_4_base_tagmut_mut_fixed/"
mut <- read_csv(paste0(filepath,"org_dump.dat")) %>% rename(param = tag_mut) %>%filter(param != 0)

print("done with mut; about to read para")
filepath<- "../../data/exp_4_base_tagmut_para_fixed/"
para <- read_csv(paste0(filepath,"org_dump.dat")) %>% rename(param = tag_mut) %>%filter(param != 0)

combo <- transform_data(mut, para)
calculate_p_vals_and_eff_sizes(combo, "_tm", "tm_fix_deltas.csv")
