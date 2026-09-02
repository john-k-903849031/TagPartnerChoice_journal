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
binned_behavior_plotter <- function(df){
  plot <- df %>%
    ggplot(aes(x=as.factor(tag_perm), fill = class, y=count)) +
    geom_col() +
    facet_wrap( ~ param, ncol = 3) + 
    facet_nested_theme +
    scale_fill_manual(name="Behavior",
                      values=c("Exclusive parasitism"="#891901",
                               "Exclusive mutualism"="#c5f6fc",
                               "Coexistence"="#66ff99",
                               "Symbiont extinction"="#cc9900"),
                      labels = function(x) str_wrap(x, width = 5)) + 
    labs(x = "Tag permissiveness", y = "Replicate count") + 
    theme(panel.spacing.x = unit(0.5,"line"),
          panel.spacing.y = unit(0.5,"line"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
          legend.direction = "horizontal",legend.justification="center" ,
          legend.box.just = "bottom", legend.position = "bottom")
  return(plot)
} 

## Experiment 2 ## 
para_syms_df <- read_csv("../../data/exp_2_base_tagmut_evolving/sym_counts.dat") %>%
  filter(update == max(update)) %>% rename(param = tag_mut)
para_syms_df <- classify_behavior(para_syms_df)

small <- para_syms_df %>% group_by(tag_perm, param, class) %>% summarize(count = n())
small <- rbind(small %>% filter(param == "none") %>% 
                 slice(rep(row_number(), length(unique(small$param)) - 1)) %>%
                 mutate(param = unique((small %>% filter(param != "none"))$param)), 
               small %>% filter(param != "none")) %>% 
  mutate(tag_perm = case_when(tag_perm == "none" ~ "no tags", .default = tag_perm),
         param = factor(paste0("TM rate ", as.numeric(param) * 100, "%")))
tag_mut_plot <- binned_behavior_plotter(small)
# 4 x 6
tag_mut_plot


## Experiment 1 ## 
vt_syms_df <- read_csv("../../data/exp_1_base_vt_evolving/sym_counts.dat") %>%
  filter(update == max(update)) %>% rename(param = vt)
vt_syms_df <- classify_behavior(vt_syms_df)

small <- vt_syms_df %>% group_by(tag_perm, param, class) %>% summarize(count = n()) %>%
  mutate(tag_perm = case_when(tag_perm == "none" ~ "no tags", .default = tag_perm),
         param = factor(paste0("VT rate ", param * 100, "%")))
vt_plot <- binned_behavior_plotter(small)
#5 x 6
vt_plot
