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
  small <- df %>% group_by(tag_perm, param, class) %>% summarize(count = n())
  
  plot <- small %>%
    ggplot(aes(x=as.factor(tag_perm), fill = class, y=count)) +
    geom_col() +
    facet_wrap( ~ paste0(facet_prefix, param)) + 
    facet_nested_theme +
    scale_fill_manual(name="Behavior",values=c("Exclusive parasitism"="#891901",
                                               "Exclusive mutualism"="#c5f6fc",
                                               "Coexistence"="#66ff99",
                                               "Symbiont extinction"="#cc9900")) + 
    labs(x = "Tag permissiveness", y = "Replicate count") + 
    theme(panel.spacing.x = unit(0.5,"line"),
          panel.spacing.y = unit(0.5,"line"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  return(plot)
} 

## parastab ##
folder <- "../../data/exp_2_base_tagmut_evolving/"
para_syms_df <- read_csv(paste0(folder,"sym_counts.dat"))
para_syms_df <- subset(para_syms_df, update==max(update))
para_syms_df <- para_syms_df %>% rename(param = tag_mut)
para_syms_df <- classify_behavior(para_syms_df)
tag_mut_plot <- binned_behavior_plotter(para_syms_df %>% filter(tag_perm != "none"), "TM rate ")

folder <- "../../data/exp_1_base_vt_evolving/"
vt_syms_df <- read_csv(paste0(folder,"sym_counts.dat"))
vt_syms_df <- subset(vt_syms_df, update==max(update))
vt_syms_df <- vt_syms_df %>% rename(param = vt)
vt_syms_df <- classify_behavior(vt_syms_df)
vt_plot <- binned_behavior_plotter(vt_syms_df %>% filter(tag_perm != "none"), "VT rate ")

tag_mut_plot
vt_plot
## combo! ##
# 12 x 4in
make_combo_plot(vt_plot, tag_mut_plot)