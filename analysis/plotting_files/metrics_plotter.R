library(ggplot2) # data visualization 
library(viridis) # data visualization 
library(ggh4x) # data visualization 
library(dplyr) # data manipulation
library(tidyverse)

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
                                              .default = "ERROR"),
                            metric = factor(metric, levels = c("Streak", "Hamming", "Hash")))
  return(small)
}

######## load data #######
## behavior - tag mutation
para_syms_df <- read_csv("../../data/exp_8_metrics_tagmut_corrected/sym_counts.dat") %>% rename(param = tag_mut)
para_syms_df <- classify_behavior(para_syms_df)

## behavior - vertical transmission
vt_syms_df <- read_csv("../../data/exp_7_metrics_vt_corrected/sym_counts.dat") %>% rename(param = vt)
vt_syms_df <- classify_behavior(vt_syms_df)


####### summary stats - behavior #######
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
  mutate(ham_advantage = case_when(Hamming > 0 | Hash > 0 ~ Hamming - Hash, .default = NaN))  %>%
  group_by(class, param_name) %>%
  mutate(test =  wilcox.test(ham_advantage, mu = 0)$p.value,
         total = 5,
         annotation = case_when(test < 0.001/total ~"***",test < 0.01/total ~"**",test < 0.05/total ~"*",.default = "NS"))

stats %>% group_by(param_name, class, param, tag_perm) %>% 
  mutate(x_label = case_when(param_name == "tag_mut" ~ "Experiment 8\n(starting parasitic)", .default = "Experiment 7\n(starting mixed)")) %>%
  filter(ham_advantage != 0) %>%
  ggplot(aes(x = class, y = ham_advantage)) +
  geom_boxplot() + 
  geom_point(alpha = 0.2, position = position_jitterdodge()) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  facet_nested_theme + 
  labs(x = "Behavior class", y = str_wrap("Per-treatment difference in count of replicates (Hamming - Hash)",width=80)) +
  theme(axis.text.x = element_text(angle = 10, hjust = 0.5, vjust=0.75),
        legend.position="bottom") +
  facet_wrap(~x_label, nrow= 2)  +
  geom_text(mapping = aes(y = 33, label = annotation), size=5)  
  
