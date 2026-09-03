library(ggplot2) # data visualization 
library(viridis) # data visualization 
library(ggh4x) # data visualization 
library(dplyr) # data manipulation
library(tidyverse)
library(effsize)

options(scipen=999)
facet_nested_theme <- theme(
  strip.background = element_rect(fill = "white", colour = "grey", linetype="dotted", linewidth=0.2), 
  panel.background = element_rect(fill='white', color='grey'),
  panel.border = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.key = element_rect(color = "transparent"),
  panel.spacing.x = unit(0,"line"),
  panel.spacing.y = unit(0,"line"))

colors <- c("Parasites"= "#891901", "Mutualists"="#c5f6fc")


########### Experiment 1 ###########
#small_combo <- rbind(read_csv(paste0("../../data/exp_1_base_vt_evolving/","org_dump.dat"),col_names=T) %>%
#                       filter(sym_tag != ""))

#small_combo <- small_combo %>% filter(sym_int < -0.2 | sym_int >= 0.2) %>% 
#  mutate(sym_bin = case_when(sym_int < -0.2 ~ "para", sym_int >= 0.2  ~ "mut", .default = "error"))

#smallsmall <- small_combo %>% group_by(sym_bin, tag_perm, vt, seed) %>%
#  summarise(mean_td = mean(tag_distance), count = n()) 

exp1_tag_df <- read_csv("../../data/exp_1_base_vt_evolving/summarized_org_dump.dat")
exp1_tag_df <- exp1_tag_df %>%
  pivot_wider(names_from = sym_bin, values_from = c(mean_td, count)) %>%
  filter(count_mut >= 50, count_para >= 50) %>%
  group_by(tag_perm, vt) %>%
  mutate(count_coexistent_reps = n(), mut_sub_para = mean_td_mut - mean_td_para) %>%
  filter(count_coexistent_reps >= 15) %>%
  mutate(vt = factor(paste0(vt*100, "%"), levels = c("20%", "30%", "40%", "50%")))

exp1_stats <- exp1_tag_df %>% 
  select(-c(count_coexistent_reps, count_para, count_mut)) %>%
  group_by(tag_perm, vt) %>%
  summarise(test = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$p.value,
            estimate = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$estimate,
            low = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$conf.int[1],
            high = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$conf.int[2]) %>%
  ungroup() %>% mutate(total = n()) %>% mutate(annotation = case_when(test < 0.0001/total ~ "p < 0.0001",
                                                                      test < 0.01/total ~ "0.0001 <= p < 0.01",
                                                                      test < 0.05/total ~ "0.01 <= p < 0.05",
                                                                      .default = "NS"))
exp1_tag_df %>% 
  ggplot(aes(x = as.factor(tag_perm), y = mut_sub_para))+ 
  geom_boxplot() + geom_point(alpha = 0.2, position = position_jitterdodge()) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  facet_nested_theme +
  facet_nested( ~ "Vertical transmission chance" + vt, scales="free_x", space = 'free') +
  theme(axis.text.x = element_text(angle = 22.5, hjust = 0.5, vjust=1),
        legend.position="bottom") + 
  labs(y = str_wrap("Difference in host tag similarity (mutualist - parasite)", width = 40),
       x="Tag permissiveness") +
  guides(color = "none") +
  geom_text(data = exp1_tag_df %>% ungroup() %>% slice_head(n=1) %>%
              mutate(vt = "20%", tag_perm = 0.5, mut_sub_para = -0.13), 
            label = str_wrap("Mutualists more similar", width = 12), 
            size = 3, hjust = 0.25) +
  geom_segment(data = exp1_tag_df %>% ungroup() %>% slice_head(n=1) %>%
                 mutate(vt = "20%", tag_perm = 0.4375, mut_sub_para = -0.075),
               aes(yend = -0.137),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "red", size = 1
  ) 


########### Experiment 2 ###########
#small_combo <- rbind(read_csv(paste0("../../data/exp_2_base_tagmut_evolving/","org_dump.dat"),col_names=T) %>% filter(sym_tag != ""))

#small_combo <- small_combo %>% filter(sym_int < -0.2 | sym_int >= 0.2) %>% 
#  mutate(tag_distance = as.numeric(tag_distance),
#         sym_bin = case_when(sym_int < -0.2 ~ "para",sym_int >= 0.2  ~ "mut",.default = "error"))

#smallsmall <- small_combo %>% group_by(sym_bin, tag_perm, tag_mut, seed) %>%
#  summarise(mean_td = mean(tag_distance), count = n())

exp2_tag_df<- read_csv( "../../data/exp_2_base_tagmut_evolving/summarized_org_dump.dat")%>%
  pivot_wider(names_from = sym_bin, values_from = c(mean_td, count)) %>%
  filter(count_mut >= 50, count_para >= 50, tag_mut != 0) %>%
  group_by(tag_perm, tag_mut) %>%
  mutate(count_coexistent_reps = n(), mut_sub_para = mean_td_mut - mean_td_para) %>%
  filter(count_coexistent_reps >= 15) %>%
  mutate(tag_mut_label = factor(paste0(tag_mut*100, "%"), levels = c("0.01%", "0.05%", "0.1%", "0.5%", "1%", "5%", "10%")))

exp2_stats <- exp2_tag_df %>% 
  select(-c(count_coexistent_reps, count_para, count_mut)) %>%
  group_by(tag_perm, tag_mut) %>%
  summarise(test = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$p.value,
            estimate = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$estimate,
            low = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$conf.int[1],
            high = wilcox.test(mut_sub_para, mu = 0, conf.int = TRUE, conf.level = 0.95)$conf.int[2]) %>%
  ungroup() %>% mutate(total = n()) %>% mutate(annotation = case_when(test < 0.0001/total ~ "p < 0.0001",
                                                                      test < 0.01/total ~ "0.0001 <= p < 0.01",
                                                                      test < 0.05/total ~ "0.01 <= p < 0.05",
                                                                      .default = "NS"))

exp2_tag_df %>% 
  ggplot(aes(x = as.factor(tag_perm), y = mut_sub_para)) + 
  geom_boxplot() + geom_point(alpha = 0.2, position = position_jitterdodge()) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  facet_nested_theme +
  facet_nested( ~ "Tag mutation chance" + tag_mut_label, scales="free_x", space = 'free') +
  theme(axis.text.x = element_text(angle = 22.5, hjust = 0.5, vjust=1)) + 
  labs(y = str_wrap("Difference in host tag similarity (mutualist - parasite)", width = 40),
       x="Tag permissiveness") +
  scale_y_continuous(breaks = c( 0, -0.07, -0.14)) +
  geom_text(data = exp2_tag_df %>%  ungroup() %>% slice_head(n=1) %>%
              mutate(tag_mut = "0.01%", tag_perm = 0.1875, mut_sub_para = -0.14), 
            label = str_wrap("Mutualists more similar", width = 12), 
            size = 3, hjust = 0.25) +
  geom_segment(data = exp2_tag_df %>% ungroup() %>% slice_head(n=1) %>%
                 mutate(tag_mut = "0.01%", tag_perm = 0.125, mut_sub_para = -0.075),
               aes(yend = -0.147),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "red", size = 1
  ) 


########### Experiment 3 ###########
#small_combo <- rbind(read_csv(paste0("../../data/exp_3_base_vt_mut_fixed/","org_dump.dat"),col_names=T) %>%
#                       filter(sym_tag != "") %>% mutate(cond = "mut"),
#                     read_csv(paste0("../../data/exp_3_base_vt_para_fixed/","org_dump.dat"),col_names=T) %>%
#                       filter(sym_tag != "") %>% mutate(cond = "para"))

# smallsmall <- small_combo %>% group_by(cond, tag_perm, vt, seed) %>%
#  summarise(mean_td = mean(tag_distance)) %>%
#  group_by(tag_perm, vt) %>%
#  add_count(cond, name = "cond_count") %>%
#  mutate(para_count = case_when(cond == "para" ~ cond_count, .default = 0),
#         mut_count = case_when(cond == "mut" ~ cond_count, .default = 0))

exp3_tag_df <- read_csv("../../data/exp_3_base_vt_mut_fixed/summarized_org_dump.dat") %>% 
  filter(vt == 0 | vt == 0.1 | vt == 0.5 | vt == 0.9 | vt == 1) %>%
  filter(tag_perm != 0.4375, tag_perm != 0.5625) %>%
  group_by(tag_perm, vt) %>% filter(max(para_count) >= 15 & max(mut_count) >= 15) %>%
  mutate(vt = factor(paste0(vt*100,"%"), levels = c( "0%", "10%", "50%", "90%", "100%"))) %>%
  mutate(cond_labels= case_when(cond == "mut" ~"Mutualists", .default ="Parasites"))

exp3_stats <- exp3_tag_df %>% 
  select(-c(para_count, mut_count, cond_count)) %>%
  group_by(tag_perm, vt) %>%
  summarise(test = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$p.value,
            estimate = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$estimate,
            low = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$conf.int[1],
            high = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$conf.int[2],
            delta = cliff.delta(mean_td ~ cond)$estimate)%>%
  ungroup() %>% mutate(total = n()) %>% mutate(annotation = case_when(test < 0.0001/total ~ "***",
                                                                      test < 0.01/total ~ "**",
                                                                      test < 0.05/total ~ "*",
                                                                      .default = "NS"))

exp3_tag_df %>%
  mutate(cond_labels= case_when(cond == "mut" ~"Mutualists", .default ="Parasites"))%>%
  ggplot(aes(x = as.factor(tag_perm), y = mean_td, color = cond_labels)) +
  geom_boxplot() + 
  geom_point(alpha = 0.2, position = position_jitterdodge()) + 
  scale_color_manual(values = colors) + 
  facet_nested_theme +
  facet_nested( ~ "Vertical transmission chance" + vt) +
  theme(axis.text.x = element_text(angle = 22.5, hjust = 1, vjust=1),
        legend.position="bottom") + 
  labs(x = "Tag permissiveness", y = "Mean host-symbiont tag distance", color = "Condition") 



########### Experiment 4 ###########
# small_combo <- rbind(read_csv(paste0("../../data/exp_4_base_tagmut_mut_fixed/","org_dump.dat"),col_names=T) %>%
#                        filter(sym_tag != "") %>% mutate(cond = "mut"),
#                      read_csv(paste0("../../data/exp_4_base_tagmut_para_fixed/","org_dump.dat"),col_names=T) %>%
#                        filter(sym_tag != "") %>% mutate(cond = "para"))
# 
# small_df <- small_combo %>% group_by(cond, tag_perm, tag_mut, seed) %>%
#  summarise(mean_td = mean(tag_distance), mean_symrepro = mean(sym_repro_count), mean_hostrepro = mean(host_repro_count )) %>%
#  group_by(tag_perm, tag_mut) %>%
#  add_count(cond, name = "cond_count") %>%
#  mutate(para_count = case_when(cond == "para" ~ cond_count, .default = 0),
#         mut_count = case_when(cond == "mut" ~ cond_count, .default = 0))
# write_csv(small_df, "../../data/exp_4_base_tagmut_mut_fixed/summarized_org_dump.dat")

exp4_tag_df <- read_csv("../../data/exp_4_base_tagmut_mut_fixed/summarized_org_dump.dat") %>%
  filter(tag_mut != 0, tag_mut != 0.001, tag_mut != 0.01) %>%
  filter(tag_perm != 0.4375, tag_perm != 0.5625) %>%
  group_by(tag_perm, tag_mut) %>% filter(max(para_count) >= 15 & max(mut_count) >= 15) %>%
  mutate(tag_mut = as.factor(paste0(as.numeric(tag_mut)*100, "%")))

exp4_stats <- exp4_tag_df %>% 
  select(-c(para_count, mut_count, cond_count)) %>%
  group_by(tag_perm, tag_mut) %>%
  summarise(test = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$p.value,
            estimate = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$estimate,
            low = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$conf.int[1],
            high = wilcox.test(mean_td ~ cond, conf.int = TRUE, conf.level = 0.95)$conf.int[2],
            delta = cliff.delta(mean_td ~ cond)$estimate) %>%
  ungroup() %>% mutate(total = n()) %>% mutate(annotation = case_when(test < 0.0001/total ~ "***",
                                                                      .default = "NS"))

exp4_tag_df %>%
  mutate(cond_labels= case_when(cond == "mut" ~"Mutualists", .default ="Parasites"))%>%
  ggplot(aes(x = as.factor(tag_perm), y = mean_td, color = cond_labels)) +
  geom_boxplot() + geom_point(alpha = 0.2, position = position_jitterdodge()) + 
  scale_color_manual(values = colors) + 
  facet_nested_theme +
  facet_nested( ~ "Tag mutation chance" + tag_mut) +
  theme(axis.text.x = element_text(angle = 22.5, hjust = 1, vjust=1),
        legend.position="bottom") + 
  labs(x = "Tag permissiveness", y = "Mean host-symbiont tag distance", color = "Condition") 

