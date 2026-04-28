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
library(ggthemes)
library(paletteer)
library(ggsignif)


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

######## classify behavior function #######
classify_behavior <- function(df){
  small <- df %>% filter(update == max(update))
  
  small <- small %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                            mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
  
  small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                              ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                              ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                              .default = "Coexistence"),
                            metric = case_when(metric == 0 ~ "Hamming", metric == 1 ~ "Streak", metric == 2 ~ "Hash",
                                               metric == 3 ~ "Raw Hamming", metric == 4 ~ "Raw Streak", metric == 5 ~ "Raw Hash", 
                                               .default = "ERROR"),
                            metric = factor(metric, levels = c("Raw Streak", "Raw Hamming", "Raw Hash", "Streak", "Hamming", "Hash")))
  return(small)
}

classify_tags <- function(df){
  small <- df %>% filter(symbiont_tag_richness > 0) %>%
    group_by(seed, tag_perm, param) %>% 
    filter(update == max(update))# grab the last surviving sym update
  
  small <- small %>% mutate(tag_dist_class = case_when(mean_tag_distance < 0.1 ~ "x < 0.1",
                                                       mean_tag_distance >= 0.1 & mean_tag_distance < 0.2 ~ "0.1 <= x < 0.2",
                                                       mean_tag_distance >= 0.2 & mean_tag_distance < 0.3 ~ "0.2 <= x < 0.3",
                                                       .default = "x >= 0.3"),
                            metric = case_when(metric == 0 ~ "Hamming", metric == 1 ~ "Streak", metric == 2 ~ "Hash",
                                               .default = "ERROR"),
                            metric = factor(metric, levels = c("Streak", "Hamming", "Hash")))
  return(small)
}


######## load data #######
## behavior - tag mutation
para_syms_df <- rbind(read_csv("../../data/exp_8_metrics_tagmut_corrected/sym_counts.dat"),
                      read_csv("../../data/exp_8_metrics_tagmut_base/sym_counts.dat") %>%
                        mutate(metric = metric + 3)) %>% 
                          rename(param = tag_mut)

para_syms_df <- classify_behavior(para_syms_df)

unique(para_syms_df$metric)
## behavior - vertical transmission
vt_syms_df <- rbind(read_csv("../../data/exp_7_metrics_vt_corrected/sym_counts.dat"),
                    read_csv("../../data/exp_7_metrics_vt_base/sym_counts.dat") %>%
                      mutate(metric = metric + 3)) %>% 
                        rename(param = vt)

vt_syms_df <- classify_behavior(vt_syms_df)

## tag distance - tag mutation
para_tags_df <- read_csv("../../data/exp_8_metrics_tagmut_corrected/tag_dists.dat") %>%
                          rename(param = tag_mut)

para_tags_df <- classify_tags(para_tags_df)

## tag distance - vertical transmission
vt_tags_df <- read_csv("../../data/exp_7_metrics_vt_corrected/tag_dists.dat") %>%
                    rename(param = vt)

vt_tags_df <- classify_tags(vt_tags_df)


######### barplot: compare raw and corrected hamming ######### 
compare_raw_vs_cor <- function(m = "Hamming"){
  plot <- rbind(vt_syms_df %>% mutate(param_name = "vt") %>% filter(metric == m | metric == paste0("Raw ",m)),
                para_syms_df %>% mutate(param_name = "tag_mut") %>% filter(metric == m | metric == paste0("Raw ",m))) %>%
    group_by(param, param_name, metric, tag_perm, class) %>%
    summarise(count = n()) %>% 
    mutate(param = as.factor(param), tag_perm = as.factor(tag_perm)) %>%
    group_by(metric, class) %>% summarize(total = sum(count)) %>%
    ggplot(aes(x = class, y = total, fill = metric)) + geom_col(position = "dodge") + 
    facet_nested_theme + scale_fill_manual(values = c(plasma(4)[3], plasma(4)[1])) + 
    labs(x = "Behavior", y = "Count of replicates", fill = "Metric") +
    theme(axis.text.x = element_text(angle = 22.5, hjust = 1, vjust=1),
          legend.position="bottom") 
  return(plot)
}
compare_raw_vs_cor("Hamming") +
compare_raw_vs_cor("Hash") +
compare_raw_vs_cor("Streak")

# 4 x 4 in


####### percent behavior class by metric ####### 
syms_combo <- rbind(vt_syms_df %>% group_by(tag_perm, metric, param, class) %>% summarize(count = n()) %>% filter(metric != "Raw Hamming"), 
                    para_syms_df %>% group_by(tag_perm, metric, param, class) %>% summarize(count = n()) %>% filter(metric != "Raw Hamming")) %>% 
  group_by(metric, class) %>% summarize(total = sum(count)) 

syms_combo %>% filter(metric == "Hamming") %>% mutate(percent = format(round(total / sum(total) * 100, 2), nsmall = 2))

syms_combo %>% filter(metric == "Streak") %>% mutate(percent = format(round(total / sum(total) * 100, 2), nsmall = 2))

syms_combo %>% filter(metric == "Hash") %>% mutate(percent = format(round(total / sum(total) * 100, 2), nsmall = 2))




######### barplot: compare behaviors ######### 
syms_plot <- syms_combo %>% group_by(metric, class) %>%
  ggplot(aes(x = class, y = total, fill = metric)) + geom_col(position = "dodge") + 
  facet_nested_theme + scale_fill_manual(values = c(turbo(4))) + 
  labs(x = "Behavior", y = "Count of replicates", fill = "Metric") +
  theme(axis.text.x = element_text(angle = 22.5, hjust = 1, vjust=1),
        legend.position="bottom")


add_p_row <- function(class_chr, source_df, p_df){
  matrix_mut <- source_df %>% filter((metric == "Hamming" | metric == "Hash")) %>%
    group_by(metric) %>% summarize(ex_mu = sum(case_when(class == class_chr ~ total, .default=0)), 
                                   not_ex_mut = sum(case_when(class != class_chr ~ total, .default=0)))
  
  mat <- as.matrix(matrix_mut[,2:3], header=TRUE,
                   row.names=1)
  rownames(mat) <- matrix_mut$metric
  return(p_df %>% add_row(class=class_chr, ham_hash_p=fisher.test(mat, alternative="two.sided")$p.value))
}
p_vals_df <- tibble(class = character(), ham_hash_p = numeric())
p_vals_df <- add_p_row("Exclusive mutualism", syms_combo, p_vals_df)
p_vals_df <- add_p_row("Exclusive parasitism", syms_combo, p_vals_df)
p_vals_df <- add_p_row("Coexistence", syms_combo, p_vals_df)
p_vals_df <- add_p_row("Symbiont extinction", syms_combo, p_vals_df)

p_vals_df <- p_vals_df %>% mutate(annotation = case_when(ham_hash_p < 0.001/8 ~ "***",
                                                         ham_hash_p < 0.01/8 ~ "**",
                                                         ham_hash_p < 0.05/8 ~ "*",
                                                         .default = "NS")) %>% arrange(class) 
syms_plot + ylim(0,3000) +
  geom_signif(
    y_position = c(1317, 2841, 1525, 493), xmin = c(1, 2, 3, 4), xmax = c(1.3, 2.3, 3.3, 4.3),
    annotation = p_vals_df$annotation, tip_length = 0.01
  )




####### percent tag class by metric ####### 
tags_combo <- rbind(vt_tags_df %>% group_by(tag_perm, metric, param, tag_dist_class) %>% summarize(count = n()) %>% filter(metric != "Raw Hamming"), 
                    para_tags_df %>% group_by(tag_perm, metric, param, tag_dist_class) %>% summarize(count = n()) %>% filter(metric != "Raw Hamming")) %>% 
  group_by(metric, tag_dist_class) %>% summarize(total = sum(count)) %>% 
  mutate(tag_dist_class = factor(tag_dist_class, levels= c("x < 0.1", "0.1 <= x < 0.2", "0.2 <= x < 0.3", "x >= 0.3")))

tags_combo %>% filter(metric == "Hamming") %>% mutate(percent = format(round(total / sum(total) * 100, 2), nsmall = 2))

tags_combo %>% filter(metric == "Streak") %>% mutate(percent = format(round(total / sum(total) * 100, 2), nsmall = 2))

tags_combo %>% filter(metric == "Hash") %>% mutate(percent = format(round(total / sum(total) * 100, 2), nsmall = 2))

######### barplot: compare tags ######### 
tags_plot <- tags_combo %>% group_by(metric, tag_dist_class) %>%
  ggplot(aes(x = tag_dist_class, y = total, fill = metric)) + geom_col(position = "dodge", width=0.9) + 
  facet_nested_theme + scale_fill_manual(values = c(turbo(4))) + 
  labs(x = "Mean final tag distance", y = "Count of replicates", fill = "Metric") +
  theme(axis.text.x = element_text(angle = 22.5, hjust = 1, vjust=1),
        legend.position="bottom")

add_p_row <- function(class_chr, source_df, p_df){
  matrix_mut <- source_df %>% filter((metric == "Hamming" | metric == "Hash")) %>%
    group_by(metric) %>% summarize(ex_mu = sum(case_when(tag_dist_class == class_chr ~ total, .default=0)), 
                                   not_ex_mut = sum(case_when(tag_dist_class != class_chr ~ total, .default=0)))
  
  mat <- as.matrix(matrix_mut[,2:3], header=TRUE,
                   row.names=1)
  rownames(mat) <- matrix_mut$metric
  return(p_df %>% add_row(tag_dist_class=class_chr, ham_hash_p=fisher.test(mat, alternative="two.sided")$p.value))
}
p_vals_df <- tibble(tag_dist_class = character(), ham_hash_p = numeric())
p_vals_df <- add_p_row("x < 0.1", tags_combo, p_vals_df)
p_vals_df <- add_p_row("0.1 <= x < 0.2", tags_combo, p_vals_df)
p_vals_df <- add_p_row("0.2 <= x < 0.3", tags_combo, p_vals_df)
p_vals_df <- add_p_row("x >= 0.3", tags_combo, p_vals_df)

p_vals_df <- p_vals_df %>% mutate(annotation = case_when(ham_hash_p < 0.001/8 ~ "***",
                                                         ham_hash_p < 0.01/8 ~ "**",
                                                         ham_hash_p < 0.05/8 ~ "*",
                                                         .default = "NS")) %>% arrange(tag_dist_class) 


tags_plot + ylim(0,3000) +
  geom_signif(
    y_position = c(2240, 2252, 1556, 288), xmin = c(1, 2, 3, 4), xmax = c(1.3, 2.3, 3.3, 4.3),
    annotation = p_vals_df$annotation, tip_length = 0.01
  )


######### compare tag distances within behavior classes ######### 
combo_syms_tags_df <- rbind(para_syms_df %>% mutate(param_name = "tm"), 
                            vt_syms_df %>% mutate(param_name = "vt")) %>% 
  group_by(tag_perm, metric, param, param_name, seed) %>% select(tag_perm, metric, param, param_name, seed, class) %>%
  left_join( rbind(vt_tags_df%>% mutate(param_name = "vt"), 
                   para_tags_df%>% mutate(param_name = "tm")) %>% select(tag_perm, metric, param, param_name, seed, tag_dist_class, mean_tag_distance) ) %>%
  mutate(tag_dist_class = factor(tag_dist_class, levels= c("x < 0.1", "0.1 <= x < 0.2", "0.2 <= x < 0.3", "x >= 0.3"))) 

work_df <- combo_syms_tags_df %>%
  filter((metric == "Hash" | metric == "Hamming") & (class == "Exclusive mutualism" | class == "Exclusive parasitism")) %>%
  group_by(metric, class, tag_dist_class) %>% summarize(count = n()) %>%
  group_by(metric, class) %>%
  mutate(percent = format(round(count / sum(count) * 100, 2), nsmall = 2))

work_df %>% filter(tag_dist_class == "x < 0.1")

