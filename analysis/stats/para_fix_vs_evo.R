
library(dplyr)
library(ggplot2)
library(viridis)
library(effsize) # stats
library(tidyverse)


transform_data <- function(a, b){
  print("binding")
  
  a <- a %>% mutate(sym_tag = case_when(sym_tag == "" ~ NA, .default = sym_tag)) %>%
    filter(sym_int < -0.6 & !is.na(sym_tag) & tag_perm != "none") %>% 
    mutate(cond = "evolved")
  
  b <- b %>% mutate(sym_tag = case_when(sym_tag == "" ~ NA, .default = sym_tag)) %>%
    filter(!is.na(sym_tag)) %>% mutate(cond = "fixed")
  
  small <- rbind(b, a)
  rm(a)
  rm(b)
  
  print("done binding and deleting; about to transform")
  
  small <- small %>% group_by(tag_perm,param) %>% mutate(fixed_pres = "fixed"%in%cond)
  small <- small %>% group_by(tag_perm,param) %>% mutate(evolved_pres = "evolved"%in%cond)
  small <- small %>% filter(evolved_pres == T & fixed_pres == T) %>% 
    mutate(tag_perm = as.factor(tag_perm), param = as.factor(param), cond = as.factor(cond))
  return(small)
}
calculate_p_vals_and_eff_sizes  <- function(df, param_prefix="",csv_name="unnamed_csv.csv"){
  nm <- c()
  md <- c()
  effect <- c()
  n_test <- 0
  for (td in unique(df$tag_perm)){
    for(v in unique(df$param)){
      a <- subset(df, tag_perm == td & param == v)
      if("fixed" %in% a$cond & "evolved"%in% a$cond){
        n_test <- n_test + 1
        nm <- append(nm, paste0("td",td,param_prefix,v))
        md <- append(md, wilcox.test(tag_distance ~ cond, data=a)$p.value)
        effect <- append(effect, cliff.delta(tag_distance ~ cond, data=a)$estimate)
        
        print("\n*******************")
        print(nm[n_test])
        print(md[n_test])
        print(effect[n_test])

	models <- data.frame(names = nm, p_vals=md, deltas=effect)

        write_csv(models, csv_name)
      }
    }
  }
  print(n_test)
  
  print(str(df$cond))
  models <- data.frame(names = nm, p_vals=md, deltas=effect)
  
  write_csv(models, csv_name)
  
  
}


## exp 1
print("** VT **")
evo <- read_csv("../../data/exp_1_base_vt_evolving/org_dump.dat") %>% rename(param = vt)
fix <- read_csv("../../data/exp_3_base_vt_para_fixed/org_dump.dat") %>% rename(param = vt)
combo <- transform_data(evo, fix)
rm(evo, fix)
calculate_p_vals_and_eff_sizes(combo, param_prefix = "_vt", "vt_para_evo_v_fix_deltas.csv")
rm(combo)

## exp 2
print("** TM **")
evo <- read_csv("../../data/exp_2_base_tagmut_evolving/org_dump.dat") %>% rename(param = tag_mut) %>% filter(param !=0)
fix <- read_csv("../../data/exp_4_base_tagmut_para_fixed/org_dump.dat") %>% rename(param = tag_mut) %>% filter(param !=0)
combo <- transform_data(evo, fix)
rm(evo, fix)
calculate_p_vals_and_eff_sizes(combo, param_prefix = "_tm", "tm_para_evo_v_fix_deltas.csv")
rm(combo)
