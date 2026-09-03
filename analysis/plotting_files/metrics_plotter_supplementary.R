library(ggplot2) # data visualization 
library(viridis) # data visualization 
library(ggh4x) # data visualization 
library(dplyr) # data manipulation
library(tidyverse)
library(patchwork)

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

######## load data #######
## behavior - tag mutation
## behavior - tag mutation
para_syms_df <- rbind(read_csv("../../data/exp_8_metrics_tagmut_corrected/sym_counts.dat"),
                      read_csv("../../data/exp_8_metrics_tagmut_base/sym_counts.dat") %>%
                        mutate(metric = metric + 3)) %>% 
  rename(param = tag_mut)

para_syms_df <- classify_behavior(para_syms_df)

## behavior - vertical transmission
vt_syms_df <- rbind(read_csv("../../data/exp_7_metrics_vt_corrected/sym_counts.dat"),
                    read_csv("../../data/exp_7_metrics_vt_base/sym_counts.dat") %>%
                      mutate(metric = metric + 3)) %>% 
  rename(param = vt)

vt_syms_df <- classify_behavior(vt_syms_df)


####### summary stats - behavior (Hamming - Streak) #######
syms_full <- rbind(vt_syms_df %>% group_by(tag_perm, metric, param, class) %>% summarize(count = n()) %>% 
                     mutate(param_name = "vt"), 
                   para_syms_df %>% group_by(tag_perm, metric, param, class) %>% summarize(count = n()) %>% 
                     mutate(param_name = "tag_mut")) %>%
  ungroup %>% 
  filter(metric == "Hamming" | metric == "Streak")%>% 
  mutate(metric = as.character(metric))


stats <- syms_full %>% ungroup() %>% group_by(param_name) %>%
  complete(metric, param, class, tag_perm, fill = list(count = 0)) %>%
  pivot_wider(id_cols = c(param, param_name, tag_perm,class), names_from = metric, values_from = count) %>%
  mutate_at(c("tag_perm","param"), as.factor) %>%
  filter(!(Hamming == 0 & Streak == 0)) %>%
  mutate(ham_advantage = case_when(Hamming > 0 | Streak > 0 ~ Hamming - Streak, .default = NaN)) 

stats %>% group_by(param_name, class, param, tag_perm) %>% 
  mutate(x_label = case_when(param_name == "tag_mut" ~ "Experiment 8\n(starting parasitic)", .default = "Experiment 7\n(starting mixed)")) %>%
  filter(ham_advantage != 0) %>%
  ggplot(aes(x = class, y = ham_advantage)) +
  geom_boxplot() + 
  geom_point(alpha = 0.2, position = position_jitterdodge()) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  facet_nested_theme + 
  labs(x = "Behavior class", y = str_wrap("Per-treatment difference in count of replicates (Hamming - Streak)",width=50)) +
  theme(axis.text.x = element_text(angle = 10, hjust = 0.5, vjust=0.75),
        legend.position="bottom") +
  facet_wrap(~x_label, nrow= 2) 

stats %>% 
  filter(ham_advantage != 0) %>%
  mutate(facet_param_name = str_wrap(case_when(param_name == "vt" ~ "Vertical transmission rate", 
                                               .default = "Tag mutation rate"), width = 20),
         param_percent = factor(paste0(as.numeric(as.character(param))*100, "%"),
                                levels = c('0%', '0.01%', '0.05%', '0.1%', '0.5%', '1%', '5%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%',  '100%'))) %>%
  ggplot(aes(x = tag_perm, y = param_percent, fill = ham_advantage)) + 
  geom_tile() + 
  facet_nested_wrap(class~facet_param_name, nrow=2, scales="free") + 
  scale_fill_gradient2(high = "#1AE4B6FF", mid = "black", low = "#FABA39FF") + 
  facet_nested_theme + 
  labs(x = "Tag permissiveness", y = "Parameter rate", fill = "Replicate\ncount\ndifference:\nHamming\n- Streak")

####### summary stats - behavior (Hamming - Hash) #######
syms_full <- rbind(vt_syms_df %>% group_by(tag_perm, metric, param, class) %>% summarize(count = n()) %>% 
                     mutate(param_name = "vt"), 
                   para_syms_df %>% group_by(tag_perm, metric, param, class) %>% summarize(count = n()) %>% 
                     mutate(param_name = "tag_mut")) %>%
  ungroup %>% 
  filter(metric == "Hamming" | metric == "Hash")%>% 
  mutate(metric = as.character(metric))


stats <- syms_full %>% ungroup() %>% group_by(param_name) %>%
  complete(metric, param, class, tag_perm, fill = list(count = 0)) %>%
  pivot_wider(id_cols = c(param, param_name, tag_perm,class), names_from = metric, values_from = count) %>%
  mutate_at(c("tag_perm","param"), as.factor) %>%
  filter(!(Hamming == 0 & Hash == 0)) %>%
  mutate(ham_advantage = case_when(Hamming > 0 | Hash > 0 ~ Hamming - Hash, .default = NaN)) 

stats %>% 
  filter(ham_advantage != 0) %>%
  mutate(facet_param_name = str_wrap(case_when(param_name == "vt" ~ "Vertical transmission rate", 
                                               .default = "Tag mutation rate"), width = 20),
         param_percent = factor(paste0(as.numeric(as.character(param))*100, "%"),
                                levels = c('0%', '0.01%', '0.05%', '0.1%', '0.5%', '1%', '5%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%',  '100%'))) %>%
  ggplot(aes(x = tag_perm, y = param_percent, fill = ham_advantage)) + 
  geom_tile() + 
  facet_nested_wrap(class~facet_param_name, nrow=2, scales="free") + 
  scale_fill_gradient2(high = "#1AE4B6FF", mid = "black", low = "#FABA39FF") + 
  facet_nested_theme + 
  labs(x = "Tag permissiveness", y = "Parameter rate", fill = "Replicate\ncount\ndifference:\nHamming\n- Hash")


######### barplot: compare raw and corrected ######### 
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


######### summary stats - tags ######### 
tags_full <- rbind(read_csv("../../data/exp_7_metrics_vt_corrected/tag_dists.dat") %>%
                     rename(param = vt) %>%
                     mutate(param_name = "vt"), 
                   read_csv("../../data/exp_8_metrics_tagmut_corrected/tag_dists.dat") %>%
                     rename(param = tag_mut) %>% 
                     mutate(param_name = "tag_mut")) %>%
  mutate(metric = case_when(metric == 0 ~ "Hamming", metric == 1 ~ "Streak", metric == 2 ~ "Hash",
                            .default = "ERROR")) %>%
  filter(metric == "Hamming" | metric == "Hash") 

work <- tags_full %>% 
  filter(update == max(update) & symbiont_tag_richness > 0) %>% 
  group_by(param, param_name, metric, tag_perm) %>%
  summarize(mean_td = mean(mean_tag_distance), count = n())%>%
  pivot_wider(id_cols = c(param, param_name, tag_perm), names_from = metric, values_from = c(count, mean_td)) %>%
  mutate(ham_sub_hash = mean_td_Hamming - mean_td_Hash) %>%
  filter(!(param_name == "tag_mut" & param == 0), 
         !is.na(count_Hamming) & !is.na(count_Hash),
         count_Hamming >= 15 & count_Hash >= 15)

stats <- work %>% group_by(param_name) %>% 
  mutate(test = wilcox.test(ham_sub_hash, mu = 0, conf.int = TRUE, conf.level = 0.95)$p.value,
         estimate = wilcox.test(ham_sub_hash, mu = 0, conf.int = TRUE, conf.level = 0.95)$estimate,
         low = wilcox.test(ham_sub_hash, mu = 0, conf.int = TRUE, conf.level = 0.95)$conf.int[1],
         high = wilcox.test(ham_sub_hash, mu = 0, conf.int = TRUE, conf.level = 0.95)$conf.int[2]) %>%
  mutate(annotation = case_when(test / 5 < 0.0001 ~ "***", .default= "other"))


stats %>% mutate(x_label = case_when(param_name == "tag_mut" ~ "Experiment 8\n(starting parasitic)", .default = "Experiment 7\n(starting mixed)")) %>%
  ggplot(aes(x = x_label, y = ham_sub_hash, group = x_label)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.2) +
  geom_hline(aes(yintercept = 0), linetype= "dashed") +
  facet_nested_theme + 
  labs(x = "", y = str_wrap("Per-treatment difference in Hamming and Hash mean tag distance (Hamming - Hash)", width = 40))  +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5)) +
  geom_text(mapping = aes(y = 0.01, label = annotation), size=5)  




############### binned metrics behavior ############### 
classify_behavior <- function(df){
  small <- df %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                            mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
  small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                              ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                              ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                              .default = "Coexistence"))
  return(small)
}
binned_behavior_plotter <- function(df, colcount = 3){
  plot <- df %>%
    ggplot(aes(x=as.factor(tag_perm), fill = class, y=count)) +
    geom_col() +
    facet_wrap( ~ param, ncol = colcount) + 
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

## Experiment 8 ## 
para_syms_df <- read_csv("../../data/exp_8_metrics_tagmut_corrected/sym_counts.dat") %>%
  filter(update == max(update)) %>% rename(param = tag_mut)
para_syms_df <- classify_behavior(para_syms_df)

small <- para_syms_df %>% group_by(tag_perm, param, class, metric) %>% summarize(count = n()) %>%
  mutate(metric = case_when(metric == 0 ~ "Hamming", metric == 1 ~ "Streak", metric == 2 ~ "Hash"),
         param = factor(paste0("TM rate ", as.numeric(param) * 100, "%")))

streak_tag_mut_plot <- binned_behavior_plotter(small %>% filter(metric == "Streak")) + ggtitle ("Streak metric")
ham_tag_mut_plot <- binned_behavior_plotter(small %>% filter(metric == "Hamming")) + ggtitle ("Hamming metric")
hash_tag_mut_plot <- binned_behavior_plotter(small %>% filter(metric == "Hash")) + ggtitle ("Hash metric")

plot <- wrap_plots(streak_tag_mut_plot+ theme(legend.position = "none"),
                   ham_tag_mut_plot+ theme(legend.position = "none"),
                   hash_tag_mut_plot+ theme(legend.position = "none"),
                   cowplot::get_legend(streak_tag_mut_plot), 
                   heights = c(1, 1,1, 0.1))
plot

## Experiment 7 ## 
vt_syms_df <- read_csv("../../data/exp_7_metrics_vt_corrected/sym_counts.dat") %>%
  filter(update == max(update)) %>% rename(param = vt)
vt_syms_df <- classify_behavior(vt_syms_df)

small <- vt_syms_df %>% group_by(tag_perm, param, class, metric) %>% summarize(count = n()) %>%
  mutate(metric = case_when(metric == 0 ~ "Hamming", metric == 1 ~ "Streak", metric == 2 ~ "Hash"),
         param = factor(paste0("VT rate ", as.numeric(param) * 100, "%")))

streak_tag_mut_plot <- binned_behavior_plotter(small %>% filter(metric == "Streak"), 4) + ggtitle ("Streak metric")
ham_tag_mut_plot <- binned_behavior_plotter(small %>% filter(metric == "Hamming"), 4) + ggtitle ("Hamming metric")
hash_tag_mut_plot <- binned_behavior_plotter(small %>% filter(metric == "Hash"), 4) + ggtitle ("Hash metric")

plot <- wrap_plots(streak_tag_mut_plot+ theme(legend.position = "none"),
                   ham_tag_mut_plot+ theme(legend.position = "none"),
                   hash_tag_mut_plot+ theme(legend.position = "none"),
                   cowplot::get_legend(streak_tag_mut_plot), 
                   heights = c(1, 1,1, 0.1))
plot


############### mean tag distance across treatments and behavior ############### 
syms_full <- rbind(read_csv("../../data/exp_7_metrics_vt_corrected/sym_counts.dat") %>%
                     rename(param = vt) %>%
                     mutate(param_name = "vt"), 
                   read_csv("../../data/exp_8_metrics_tagmut_corrected/sym_counts.dat") %>%
                     rename(param = tag_mut) %>% 
                     mutate(param_name = "tag_mut")) 

tags_full <- rbind(read_csv("../../data/exp_7_metrics_vt_corrected/tag_dists.dat") %>%
                     rename(param = vt) %>%
                     mutate(param_name = "vt"), 
                   read_csv("../../data/exp_8_metrics_tagmut_corrected/tag_dists.dat") %>%
                     rename(param = tag_mut) %>% 
                     mutate(param_name = "tag_mut")) %>%
  mutate(metric = case_when(metric == 0 ~ "Hamming", metric == 1 ~ "Streak", metric == 2 ~ "Hash",
                            .default = "ERROR")) %>%
  filter(metric == "Hamming" | metric == "Hash") 

combo_df <- classify_behavior(syms_full) 
combo_df <- combo_df %>% left_join(tags_full %>% filter(update == max(update)))

work <- combo_df %>% select(c(seed, param, tag_perm, param_name, metric, count, mean_tag_distance, class)) %>%
  filter(metric == "Hamming" | metric == "Hash") %>%
  group_by(param, tag_perm, param_name, metric, class) %>%
  group_by(param, tag_perm, param_name, class) %>%
  mutate(ham_count = sum(metric == "Hamming"), hash_count = sum(metric == "Hash")) %>%
  filter(ham_count >= 5, hash_count >= 5)

exp78_td_full_theme <- facet_nested_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust =1, vjust=1),
        legend.position="bottom") 

nested_tag_distance_plotter <- function(df, fancy_strip){
  df %>% ggplot(aes(x = as.factor(tag_perm), y = mean_tag_distance,  color = metric)) +
    geom_boxplot() + geom_point(position=position_jitterdodge(), alpha = 0.2) +
    exp78_td_full_theme +
    facet_nested(~class +  param, scales="free_x", space = "free", 
                 strip  = fancy_strip) +
    labs(y = str_wrap("Difference in tag distance (Hamming - Hash)",30),x="Tag permissiveness", color = "Metric")+
    scale_color_manual(values=c("Hash"="#0D0887FF", "Hamming"="#F48849FF" )) 
}


top <- nested_tag_distance_plotter(
  work %>% filter(param_name == "vt") %>% 
    filter(!is.na(mean_tag_distance), class != "Exclusive mutualism")%>%
    mutate(param = factor(paste0("VT rate ", as.numeric(param) * 100, "%"))),
  strip_nested(text_x = elem_list_text(colour = c("black", "white","black","black","black","black", "black")),
               background_x = elem_list_rect(fill = c("#66ff99", "#891901", "white", "white", "white", "white","white", "white")))
) + ggtitle("Experiment 7 (starting mixed)")


bottom <-nested_tag_distance_plotter(
  work %>% filter(param_name == "vt") %>% 
    filter(!is.na(mean_tag_distance), class == "Exclusive mutualism")%>%
    mutate(param = factor(paste0("VT rate ", as.numeric(param) * 100, "%"))),
  strip_nested(text_x = elem_list_text(colour = c("black","black","black","black","black", "black")),
               background_x = elem_list_rect(fill = c("#c5f6fc", "white", "white", "white", "white","white", "white")))
)


exp7_plot <- wrap_plots(cowplot::get_title(top),
                        top+ theme(legend.position = "none") + ggtitle("") + xlab(""),
                        bottom+ theme(legend.position = "none"),
                        cowplot::get_legend(top), 
                        heights = c(0.1, 1, 1, 0.1))
exp7_plot


top <- nested_tag_distance_plotter(
  work %>% filter(param_name == "tag_mut") %>% 
    filter(param != 0, !is.na(mean_tag_distance), class !="Exclusive parasitism")%>%
    mutate(param = factor(paste0("TM rate ", as.numeric(param) * 100, "%"))),
  strip_nested(text_x = elem_list_text(colour = c("black", "black","black","black","black", "black")),
    background_x = elem_list_rect(fill = c("#66ff99", "#c5f6fc", "white", "white", "white","white", "white")))
) + ggtitle("Experiment 8 (starting parasitic)")


bottom <- nested_tag_distance_plotter(
  work %>% filter(param_name == "tag_mut") %>% 
    filter(param != 0, !is.na(mean_tag_distance), class =="Exclusive parasitism")%>%
    mutate(param = factor(paste0("TM rate ", as.numeric(param) * 100, "%"))),
  strip_nested(text_x = elem_list_text(colour = c("white", "black","black","black","black", "black")),
    background_x = elem_list_rect(fill = c("#891901", "white", "white", "white","white", "white")))
)

exp8_plot <- wrap_plots(cowplot::get_title(top),
                        top+ theme(legend.position = "none") + ggtitle("") + xlab(""),
                        bottom+ theme(legend.position = "none"),
                        cowplot::get_legend(top), 
                        heights = c(0.1, 1, 1, 0.1))
exp8_plot










