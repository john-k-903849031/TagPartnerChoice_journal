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

tenhelix <- c("#891901", "#B50142", "#D506AD", "#AB08FF", "#5731FD", "#4755FF", "#42bcf5", "#42e0f5", "#9af1fc", "#c5f6fc")
cat_names <- c("-1 to -0.8 (Parasitic)", "-0.8 to -0.6 (Parasitic)", "-0.6 to -0.4 (Detrimental)", 
                       "-0.4 to -0.2 (Detrimental)", "-0.2 to 0 (Nearly Neutral)", "0 to 0.2 (Nearly Neutral)", 
                       "0.2 to 0.4 (Positive)", "0.4 to 0.6 (Positive)", 
                       "0.6 to 0.8 (Mutualistic)", "0.8 to 1.0 (Mutualistic)")

#################################################
### mean tag distance over time, experiment 1 ###
#################################################

folder <- "../../data/exp_1_base_vt_evolving/"
tags_df <- read_csv(paste0(folder, "tag_dists.dat")) %>% mutate_at(c("seed", "vt", "tag_perm"), as.numeric)
syms_df <- read_csv(paste0(folder, "sym_counts.dat"))%>% filter(tag_perm != "none")%>% mutate_at(c("seed", "vt", "tag_perm"), as.numeric)

combo_df <- left_join(tags_df %>% mutate_at(c("seed", "vt", "tag_perm"), as.factor), 
                      syms_df %>% mutate_at(c("seed", "vt", "tag_perm"), as.factor))


combo_df <- combo_df %>% filter(count > 0) %>%
  mutate(bin = cut(mean_intval, breaks=seq(-1, 1, 0.2), right = FALSE),
         bin = case_when(mean_intval == 1 | bin == "[0.8,1)" ~ "[0.8,1]", .default = bin),
         bin = factor(bin, levels = c("[-1,-0.8)", "[-0.8,-0.6)", "[-0.6,-0.4)", 
                                      "[-0.4,-0.2)", "[-0.2,0)", "[0,0.2)", "[0.2,0.4)", 
                                      "[0.4,0.6)", "[0.6,0.8)", "[0.8,1]")))

combo_df %>% 
  filter(count > 0) %>%
  ggplot(aes(x = update, y = mean_tag_distance, group = seed, color = bin)) + 
  geom_line(na.rm=T, linewidth=0.2) + 
  labs(x = "Timesteps", y = "Mean tag distance") + 
  scale_y_continuous(breaks=c(0, 0.2, 0.4)) + 
  scale_x_continuous(breaks=c(0, 25000, 50000, 75000)) + 
  facet_nested("Tag permissiveness" + tag_perm ~ "Vertical transmission rate" + vt) + 
  facet_nested_theme +
  scale_color_manual(name="Mean symbiont\ninteraction values",values = tenhelix) +
  theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom") +
  theme(legend.position = 'bottom')


#################################################
### mean tag distance over time, experiment 2 ###
#################################################

folder <- "../../data/exp_2_base_tagmut_evolving/"
tags_df <- read_csv(paste0(folder, "tag_dists.dat")) %>% mutate_at(c("seed", "tag_mut", "tag_perm"), as.numeric)
syms_df <- read_csv(paste0(folder, "sym_counts.dat"))%>% filter(tag_perm != "none")%>% mutate_at(c("seed", "tag_mut", "tag_perm"), as.numeric)

combo_df <- left_join(tags_df %>% mutate_at(c("seed", "tag_mut", "tag_perm"), as.factor), 
                      syms_df %>% mutate_at(c("seed", "tag_mut", "tag_perm"), as.factor))

combo_df <- combo_df %>% filter(count > 0) %>%
  mutate(bin = cut(mean_intval, breaks=seq(-1, 1, 0.2), right = FALSE),
         bin = case_when(mean_intval == 1 | bin == "[0.8,1)" ~ "[0.8,1]", .default = bin),
         bin = factor(bin, levels = c("[-1,-0.8)", "[-0.8,-0.6)", "[-0.6,-0.4)", 
                                         "[-0.4,-0.2)", "[-0.2,0)", "[0,0.2)", "[0.2,0.4)", 
                                         "[0.4,0.6)", "[0.6,0.8)", "[0.8,1]")))
combo_df %>% 
  filter(count > 0) %>%
  ggplot(aes(x = update, y = mean_tag_distance, group = seed, color = bin)) + 
  geom_line(na.rm=T, linewidth=0.2) + 
  labs(x = "Timesteps", y = "Mean tag distance") + 
  scale_y_continuous(breaks=c(0, 0.2, 0.4)) + 
  scale_x_continuous(breaks=c(0, 25000, 50000, 75000)) + 
  facet_nested("Tag permissiveness" + tag_perm ~ "Tag mutation rate" + tag_mut) + 
  facet_nested_theme +
  scale_color_manual(name="Mean symbiont\ninteraction values",values = tenhelix) +
  theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom") +
  theme(legend.position = 'bottom')


#################################################
### mean tag distance over time, experiment 3 ###
#################################################

mut_df <- read_csv("../../data/exp_3_base_vt_mut_fixed/tag_dists.dat") %>% 
  mutate_at(c("seed", "vt", "tag_perm"), as.numeric) %>%
  mutate(bin = "Fixed mutualist [0.8,1]")
para_df <- read_csv("../../data/exp_3_base_vt_para_fixed/tag_dists.dat") %>% 
  mutate_at(c("seed", "vt", "tag_perm"), as.numeric) %>%
  mutate(bin = "Fixed parasite [-1,-0.8)")

combo_df <- rbind(mut_df, para_df) %>% mutate_at(c("seed", "vt", "tag_perm"), as.factor)

combo_df %>% 
  filter(symbiont_tag_richness > 0) %>% 
  ggplot(aes(x = update, y = mean_tag_distance, group = seed, color = bin)) + 
  geom_line(na.rm=T, linewidth=0.2) + 
  labs(x = "Timesteps", y = "Mean tag distance") + 
  scale_y_continuous(breaks=c(0, 0.2, 0.4)) + 
  scale_x_continuous(breaks=c(0, 25000, 50000, 75000)) + 
  facet_nested("Tag permissiveness" + tag_perm ~ "Vertical transmission rate" + vt) + 
  facet_nested_theme +
  scale_color_manual(name="Symbiont\ninteraction values",values = c("Fixed mutualist [0.8,1]"="#c5f6fc",
                                                                         "Fixed parasite [-1,-0.8)"="#891901")) +
  theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom") +
  theme(legend.position = 'bottom')

#################################################
### mean tag distance over time, experiment 4 ###
#################################################

mut_df <- read_csv("../../data/exp_4_base_tagmut_mut_fixed/tag_dists.dat") %>% 
  mutate_at(c("seed", "tag_mut", "tag_perm"), as.numeric) %>%
  mutate(bin = "Fixed mutualist [0.8,1]")
para_df <- read_csv("../../data/exp_4_base_tagmut_para_fixed/tag_dists.dat") %>% 
  mutate_at(c("seed", "tag_mut", "tag_perm"), as.numeric) %>%
  mutate(bin = "Fixed parasite [-1,-0.8)")

combo_df <- rbind(mut_df, para_df) %>% mutate_at(c("seed", "tag_mut", "tag_perm"), as.factor)

combo_df %>% 
  filter(symbiont_tag_richness > 0) %>%
  ggplot(aes(x = update, y = mean_tag_distance, group = seed, color = bin)) + 
  geom_line(na.rm=T, linewidth=0.2) + 
  labs(x = "Timesteps", y = "Mean tag distance") + 
  scale_y_continuous(breaks=c(0, 0.2, 0.4)) + 
  scale_x_continuous(breaks=c(0, 25000, 50000, 75000)) + 
  facet_nested("Tag permissiveness" + tag_perm ~ "Tag mutation rate" + tag_mut) + 
  facet_nested_theme +
  scale_color_manual(name="Symbiont\ninteraction values",values = c("Fixed mutualist [0.8,1]"="#c5f6fc",
                                                                    "Fixed parasite [-1,-0.8)"="#891901")) +
  theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom") +
  theme(legend.position = 'bottom')
