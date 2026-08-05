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

## themes and colors ##
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

perm_color_map <- c("x <= 0" = "#440154FF",
                    "0 < x <= 0.25" = "#414487FF", 
                    "0.25 < x <= 0.5" = "#2A788EFF", 
                    "0.5 < x <= 0.75" = "#22A884FF", 
                    "0.75 < x <= 1" = "#7AD151FF",             
                    "x > 1" = "#FDE725FF") 

## classification methods ##
classify_behavior <- function(df){
  small <- df %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                            mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
  small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                              ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                              ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                              .default = "Coexistence"))
  return(small)
}
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

####### interaction values ####### 
## parastab ##
folder <- "../../data/exp_6_evoperm_tagmut/"
para_syms_df <- read_csv(paste0(folder,"sym_counts.dat")) %>% 
  filter(update == max(update)) %>% rename(param = tag_mut) %>% mutate(param_name = "tag_mut")
para_syms_df <- classify_behavior(para_syms_df)

## vt ## 
folder <- "../../data/exp_5_evoperm_vt/"
vt_syms_df <- read_csv(paste0(folder,"sym_counts.dat")) %>% 
  filter(update == max(update)) %>% rename(param = vt) %>% mutate(param_name = "vt")
vt_syms_df <- classify_behavior(vt_syms_df)

## combo ## 
combo_syms_df <- rbind(vt_syms_df, para_syms_df) %>% select(param_name, param, evo_mut, seed, class) %>% rename(behavior_class = class)

####### tags ####### 
## parastab ## 
folder <- "../../data/exp_6_evoperm_tagmut/"
para_tags_df <- read_csv(paste0(folder,"tag_dists.dat")) %>% 
  rename(param = tag_mut) %>% mutate(param_name = "tag_mut")
para_tags_df <- classify_permissivenss(para_tags_df)

## vt ## 
folder <- "../../data/exp_5_evoperm_vt/"
vt_tags_df <- read_csv(paste0(folder,"tag_dists.dat")) %>% 
  rename(param = vt) %>% mutate(param_name = "vt")
vt_tags_df <- classify_permissivenss(vt_tags_df)

## combo ## 
combo_tags_df <- rbind(para_tags_df, vt_tags_df) %>% select(param_name, param, evo_mut, seed, class)


##### combined behavior and tags ##### 
combo_df <- combo_syms_df %>% group_by(param_name, param, evo_mut, seed) %>% 
  left_join(combo_tags_df) %>%
  group_by(param_name, param, evo_mut, class, behavior_class) %>%
  summarize(count = n()) %>% 
  ungroup() %>%
  group_by(param_name, param, evo_mut, behavior_class) %>%
  mutate(behavior_class_count = sum(count), evo_mut = as.factor(evo_mut))


##### plot replicate class and LSU permissivness #####
rep_behavior_perm_plotter <- function(param_name_val = "tag_mut", 
                                      facet_prefix = "TM rate ",
                                      behavior_class_val = "Coexistence"){
  return(
      combo_df %>% filter(param_name == param_name_val, behavior_class == behavior_class_val) %>%
      ggplot(aes(x = as.factor(evo_mut), y = count/behavior_class_count, fill = class)) +
      geom_col() + 
      scale_fill_manual(name="Mean final\npermissiveness", values = perm_color_map) + 
      labs(x = "Permissiveness mutation size", 
           y = paste0("Proportion of ",str_to_lower(behavior_class_val)," replicates")) +
      facet_wrap(~paste0(facet_prefix, param)) +
      facet_nested_theme +
      theme(panel.spacing.x = unit(0.5,"line"),
            panel.spacing.y = unit(0.5,"line"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  )
}

### Tag mutation plot ###
coex_plot <- rep_behavior_perm_plotter("tag_mut", "TM rate ", "Coexistence")
mut_plot <- rep_behavior_perm_plotter("tag_mut", "TM rate ", "Exclusive mutualism")
extinct_plot <- rep_behavior_perm_plotter("tag_mut", "TM rate ", "Symbiont extinction")
para_plot <- rep_behavior_perm_plotter("tag_mut", "TM rate ", "Exclusive parasitism")
# order: coex, mut, para, extinction

legend <- cowplot::get_legend(coex_plot+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
plot <- wrap_plots((((coex_plot+ theme(legend.position = "none")) + 
                    (mut_plot+ theme(legend.position = "none"))) /
                    ((para_plot+ theme(legend.position = "none")) +
                    (extinct_plot+ theme(legend.position = "none")))), 
                   legend, 
                   heights = c(1, 0.1)) + plot_annotation(tag_levels  = list(c("a)", "b)", "c)", "d)", "")))
plot

### VT plot ###
coex_plot <- rep_behavior_perm_plotter("vt", "VT rate ", "Coexistence")
mut_plot <- rep_behavior_perm_plotter("vt", "VT rate ", "Exclusive mutualism")
extinct_plot <- rep_behavior_perm_plotter("vt", "VT rate ", "Symbiont extinction")
para_plot <- rep_behavior_perm_plotter("vt", "VT rate ", "Exclusive parasitism")
# order: coex, mut, para, extinction

legend <- cowplot::get_legend(coex_plot+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
plot <- wrap_plots((((coex_plot+ theme(legend.position = "none")) + 
                       (mut_plot+ theme(legend.position = "none"))) /
                      ((para_plot+ theme(legend.position = "none")) +
                         (extinct_plot+ theme(legend.position = "none")))), 
                   legend, 
                   heights = c(1, 0.1)) + plot_annotation(tag_levels  = list(c("a)", "b)", "c)", "d)", "")))
plot