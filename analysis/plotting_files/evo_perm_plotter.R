library(ggplot2) # data visualization 
library(viridis) # data visualization 
library(ggh4x) # data visualization 
library(dplyr) # data manipulation
library(tidyverse)
library(cowplot)
library(gridExtra)
library(ggplotify)
library(patchwork)
library(ggpattern)

## themes and colors
options(scipen=999)
facet_nested_theme <- theme(
  strip.background = element_rect(fill = "white", colour = "grey", linetype="dotted", linewidth=0.2), 
  panel.background = element_rect(fill='white', color='grey'),
  panel.border = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.key = element_rect(color = "transparent"),
  panel.spacing.x = unit(0,"line"),
  panel.spacing.y = unit(0,"line"),
  axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.2))

make_combo_plot <- function(vt_sweep_plot, tag_mut_plot){
  legend <- cowplot::get_legend(vt_sweep_plot+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
  plot <- wrap_plots(((vt_sweep_plot+ theme(legend.position = "none")) + (tag_mut_plot+ theme(legend.position = "none"))), 
             legend, 
             heights = c(1, 0.1)) + plot_annotation(tag_levels  = list(c("a)", "b)", "")))
  return(plot)
}

########### binned behavior ###########
classify_behavior <- function(df){
  small <- df %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                            mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
  small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                              ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                              ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                              .default = "Coexistence"))
  return(small)
}
binned_behavior_plotter <- function(df, facet_prefix = ""){
  small <- df %>% group_by(evo_mut, param, class) %>% summarize(count = n())
  
  plot <- small %>%
    ggplot(aes(x=as.factor(evo_mut), fill = class, y=count)) +
    geom_col() +
    facet_wrap( ~ paste0(facet_prefix, param)) + 
    facet_nested_theme +
    scale_fill_manual(name="Behavior",values=c("Exclusive parasitism"="#891901",
                                               "Exclusive mutualism"="#c5f6fc",
                                               "Coexistence"="#66ff99",
                                               "Symbiont extinction"="#cc9900")) + 
    labs(x = "Permissiveness mutation size", y = "Replicate count") + 
    theme(panel.spacing.x = unit(0.5,"line"),
          panel.spacing.y = unit(0.5,"line"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  return(plot)
} 

## parastab ##
folder <- "../../data/exp_6_evoperm_tagmut/"
para_syms_df <- read_csv(paste0(folder,"sym_counts.dat"))
para_syms_df <- subset(para_syms_df, update==max(update))
para_syms_df <- para_syms_df %>% rename(param = tag_mut)
para_syms_df <- classify_behavior(para_syms_df)
tag_mut_plot <- binned_behavior_plotter(para_syms_df, "TM rate ")

## vt ## 
folder <- "../../data/exp_5_evoperm_vt/"
vt_syms_df <- read_csv(paste0(folder,"sym_counts.dat"))
vt_syms_df <- subset(vt_syms_df, update==max(update))
vt_syms_df <- vt_syms_df %>% rename(param = vt)
vt_syms_df <- classify_behavior(vt_syms_df)
vt_sweep_plot <- binned_behavior_plotter(vt_syms_df, "VT rate ")

## combo! ##
# 12 x 4in
make_combo_plot(vt_sweep_plot, tag_mut_plot)

## looking at counts!
para_syms_df$param_name <- "tag_mut"
vt_syms_df$param_name <- "vt"
combo <- rbind(vt_syms_df, para_syms_df)

# Symbionts do not evolve exclusive parasitism when permissiveness mutation size is greater than 0.0001
unique((combo %>% filter(class == "Exclusive parasitism"))$evo_mut)

# count of other classes 
combo<- combo %>% group_by(param, param_name, evo_mut) %>%  mutate(count = 1) %>%
  pivot_wider(names_from = class, values_from = count, values_fill = 0) 

combo <- combo %>% group_by(param, param_name) %>% 
  filter(sum(`Exclusive parasitism`) > 0) %>%
  filter(evo_mut != 0 & evo_mut != 0.0001) 

ep <- sum(combo$`Exclusive parasitism`)
em <- sum(combo$`Exclusive mutualism`)
c <- sum(combo$Coexistence)
se <- sum(combo$`Symbiont extinction`)
total <- ep + em + c + se

round(ep / total*100,2) # para prop 
round(se / total*100,2) # extinction prop
round(em / total*100,2) # mutualism prop
round(c / total*100,2)  # coexistence prop

rm(ep, em, c, se, total)



########### mean permissiveness ###########
classify_permissivenss <- function(df){
  small <- df %>% group_by(seed,evo_mut,param) %>%
    filter(symbiont_tag_richness > 0) %>%
    group_by(seed,evo_mut,param) %>%
    filter(update==max(update))
  small <- small %>% group_by(seed) %>% mutate(class = case_when(mean_host_perm <= 0 ~ "x <= 0",
                                                              mean_host_perm <= 0.25 & mean_host_perm > 0 ~ "0 < x <= 0.25",
                                                              mean_host_perm <= 0.5 & mean_host_perm > 0.25 ~ "0.25 < x <= 0.5",
                                                              mean_host_perm <= 0.75 & mean_host_perm > 0.5 ~ "0.5 < x <= 0.75",
                                                              mean_host_perm <= 1 & mean_host_perm > 0.75 ~ "0.75 < x <= 1",
                                                              mean_host_perm > 1 ~ "x > 1", .default = "ERROR"))
  small$class <- factor(small$class, levels=c("x <= 0", 
                                              "0 < x <= 0.25", 
                                              "0.25 < x <= 0.5", 
                                              "0.5 < x <= 0.75", 
                                              "0.75 < x <= 1",             
                                              "x > 1" ))
  return(small)
}
lsu_tag_dist <- function(df, facet_prefix = ""){
  small <- df %>% group_by(evo_mut, param, class) %>% summarize(count = n())
  
  plot <- small %>%
    ggplot(aes(x=as.factor(evo_mut), fill = class, y=count)) +
    geom_col() +
    facet_wrap( ~ paste0(facet_prefix, param)) + 
    facet_nested_theme +
    scale_fill_viridis_d(name="Mean final\npermissiveness") + 
    labs(x = "Permissiveness mutation size", y = "Replicate count") + 
    theme(panel.spacing.x = unit(0.5,"line"),
          panel.spacing.y = unit(0.5,"line"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  return(plot)
}

## parastab ## 
folder <- "../../data/exp_6_evoperm_tagmut/"
para_tags_df <- read_csv(paste0(folder,"tag_dists.dat"))
para_tags_df <- para_tags_df %>% rename(param = tag_mut)
para_tags_df <- classify_permissivenss(para_tags_df)
tag_mut_plot <- lsu_tag_dist(para_tags_df, "TM rate ")

## vt ## 
folder <- "../../data/exp_5_evoperm_vt/"
vt_tags_df <- read_csv(paste0(folder,"tag_dists.dat"))
vt_tags_df <- vt_tags_df %>% rename(param = vt)
vt_tags_df <- classify_permissivenss(vt_tags_df)
vt_sweep_plot <- lsu_tag_dist(vt_tags_df, "VT rate ")

## combo! ##
# 12 x 4in
make_combo_plot(vt_sweep_plot, tag_mut_plot)


########### combining behavior and mean permissiveness ###########
para_combo_df <- para_syms_df %>% rename(behavior_class = class) %>%
  select(-c(update, `Hist_-1`:Hist_0.9, mean_intval:min_intval)) %>%
  left_join(para_tags_df %>% select(seed, evo_mut, param, mean_host_perm, class)) %>%
  group_by(param) %>% filter(sum(behavior_class == "Exclusive parasitism") > 0)

vt_combo_df <- vt_syms_df %>% rename(behavior_class = class) %>%
  select(-c(update, `Hist_-1`:Hist_0.9, mean_intval:min_intval)) %>%
  left_join(vt_tags_df %>% select(seed, evo_mut, param, mean_host_perm, class)) %>%
  group_by(param) %>% filter(sum(behavior_class == "Exclusive parasitism") > 0)

unique((para_combo_df %>% filter(behavior_class == "Exclusive parasitism"))$class)
unique((vt_combo_df %>% filter(behavior_class == "Exclusive parasitism"))$class)

small <- para_combo_df %>% filter(behavior_class == "Symbiont extinction" & evo_mut != 0 & evo_mut != 0.0001) 
table(small$class)
l0 <- nrow(small %>% filter(class == "x <= 0"))
l0.25 <- nrow(small %>% filter(class == "0 < x <= 0.25"))
l0.5 <- nrow(small %>% filter(class == "0.25 < x <= 0.5"))

small <- vt_combo_df %>% filter(behavior_class == "Symbiont extinction" & evo_mut != 0 & evo_mut != 0.0001) 
table(small$class)
l0 <- l0 +nrow(small %>% filter(class == "x <= 0"))
l0.25 <- l0.25 + nrow(small %>% filter(class == "0 < x <= 0.25"))
l0.5 <- l0.5 + nrow(small %>% filter(class == "0.25 < x <= 0.5"))

round(l0 / (l0 + l0.25 + l0.5)*100,2)
round(l0.25 / (l0 + l0.25 + l0.5)*100,2)
round(l0.5 / (l0 + l0.25 + l0.5)*100,2)
l0 + l0.25 + l0.5 # should be 1197



############ High permissiveness mutation can drive mutualisms extinct ############
# count of classes
small <- para_syms_df %>% filter(param >= 0.0001 & param <= 0.001) %>% group_by(param, evo_mut, class) %>% 
  summarise(count = n())# %>% pivot_wider(names_from = class, values_from = count,values_fill = 0)

small%>% filter(evo_mut <= 0.001) %>% group_by(class) %>% summarize(total = sum(count))
small%>% filter(evo_mut > 0.001) %>% group_by(class) %>% summarize(total = sum(count)) %>% ungroup() %>%
  mutate(percent = as.character(round(total / sum(total)*100, 2)))

# grab LSU from the behavior csv
lsu_para_syms_df <- read_csv(paste0("../../data/exp_6_evoperm_tagmut/","sym_counts.dat")) %>% rename(param = tag_mut) %>%
  group_by(seed,evo_mut,param) %>%
  filter(count > 0) %>%
  filter(update==max(update)) %>% 
  left_join(para_tags_df %>% rename(perm_class = class)) %>%
  filter(param >= 0.0001 & param <= 0.001)
lsu_para_syms_df <- classify_behavior(lsu_para_syms_df)
binned_behavior_plotter(lsu_para_syms_df)

tags_small <- lsu_para_syms_df %>% rename(pre_ex_class = class) %>% 
  select(-c(update:Hist_0.9,mean_tag_distance:mean_host_perm, para_count, mut_count)) %>%
  group_by(param, evo_mut, seed) %>%
  left_join(para_syms_df %>% select(param, evo_mut, seed, class) %>% rename(end_exp_class = class))

nrow(tags_small %>% filter(end_exp_class == "Symbiont extinction")) # should be 267

tags_small %>% filter(end_exp_class == "Symbiont extinction") %>% group_by(perm_class) %>% 
  summarize(class_count = n()) %>% ungroup() %>% mutate(percent = as.character(round(class_count/sum(class_count) * 100, 2))) 
tags_small %>% filter(evo_mut <= 0.001) %>% group_by(perm_class) %>% 
  summarize(class_count = n()) %>% ungroup() %>% mutate(percent = as.character(round(class_count/sum(class_count) * 100, 2))) 
rm(small, tags_small)

############ Evolving permissiveness increases mutualist--parasite coexistence ############
combo_df <- rbind(para_syms_df %>% filter(param == 0.0005) %>% mutate(param_name = "tag_mut"),
                  vt_syms_df %>% filter(param == 0.5) %>% mutate(param_name = "vt")) %>%
                  select(-c(update:Hist_0.9)) %>% rename(behavior_class = class)

### counts of behaviors 
# mutualism
table((combo_df %>% filter(param_name == "vt" & param == 0.5 & evo_mut <= 0.001))$behavior_class)
table((combo_df %>% filter(param_name == "tag_mut" & param == 0.0005 & evo_mut < 0.001))$behavior_class)

# coexistence
table((combo_df %>% filter(param_name == "vt" & param == 0.5 & evo_mut > 0.001))$behavior_class)
table((combo_df %>% filter(param_name == "tag_mut" & param == 0.0005 & evo_mut >= 0.001 & evo_mut < 0.01))$behavior_class)


### counts of tag distances 
combo_tag_df <- rbind(para_tags_df %>% filter(param == 0.0005) %>% mutate(param_name = "tag_mut"),
                       vt_tags_df  %>% filter(param == 0.5) %>% mutate(param_name = "vt")) %>%
   select(-c(mean_tag_distance:mean_host_perm)) %>% rename(perm_class = class)

combo_df <- combo_df %>% group_by(param_name, param, evo_mut, seed) %>% 
  left_join(combo_tag_df)

table((combo_df %>% filter((param_name == "vt" & param == 0.5 & evo_mut <= 0.001 & behavior_class == "Exclusive mutualism") | 
                           (param_name == "tag_mut" & param == 0.0005 & evo_mut < 0.001 & behavior_class == "Exclusive mutualism")))$perm_class)

table((combo_df %>% filter((param_name == "vt" & param == 0.5 & evo_mut > 0.001 & behavior_class == "Coexistence") | 
                             (param_name == "tag_mut" & param == 0.0005 & evo_mut >= 0.001 & evo_mut < 0.01 & behavior_class == "Coexistence")))$perm_class)


sum((combo_df %>% filter((param_name == "vt" & param == 0.5 & evo_mut <= 0.001 & behavior_class == "Exclusive mutualism") | 
                             (param_name == "tag_mut" & param == 0.0005 & evo_mut < 0.001 & behavior_class == "Exclusive mutualism")))$mut_count)

small <- combo_df
small <- small %>% mutate(mut_behavior_class = case_when(
  ((param_name == "vt" & param == 0.5 & evo_mut <= 0.001 & behavior_class == "Exclusive mutualism") | 
  (param_name == "tag_mut" & param == 0.0005 & evo_mut < 0.001 & behavior_class == "Exclusive mutualism")) ~ "Low-mutation\nmutualism",
  ((param_name == "vt" & param == 0.5 & evo_mut > 0.001 & behavior_class == "Coexistence") | 
     (param_name == "tag_mut" & param == 0.0005 & evo_mut >= 0.001 & evo_mut < 0.01 & behavior_class == "Coexistence")) ~ "High-mutation\ncoexistence",
  .default = "other"))
mut_plot <- small %>% filter(mut_behavior_class != "other") %>%
  ggplot(aes(x = as.factor(evo_mut), y = mut_count, color=mut_behavior_class)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.2) +
  labs(x = "Permissiveness mutation size", y = "Mutualist count", color = "Condition") + facet_nested_theme + 
  ylim(0, 6000) 

para_plot <- small %>% filter(mut_behavior_class != "other") %>%
  ggplot(aes(x = as.factor(evo_mut), y = para_count, color=mut_behavior_class)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.2) +
  labs(x = "Permissiveness mutation size", y = "Parasite count") + facet_nested_theme + 
  ylim(0, 6000) 

make_combo_plot(mut_plot, para_plot)
