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

cutoffs_df <- read_csv("methods_fig_data/means.csv") %>% 
  rename(perm_mean = mean) %>%
  group_by(perm_mean, cutoff) %>% 
  summarize(cutoff_count = n()) %>% 
  ungroup() %>% 
  complete(perm_mean, cutoff, fill = list(cutoff_count = 0)) %>%
  group_by(perm_mean) %>%
  arrange(desc(cutoff)) %>%
  mutate(cumulative_count = cumsum(cutoff_count)) %>%
  filter(cutoff < 33)

small <- cutoffs_df %>%  
  mutate(cumulative_rate = cumulative_count/1000000,
         perm_mean = as.factor(perm_mean))

small %>% 
  ggplot(aes(x = cutoff/32, y = perm_mean, fill = cumulative_rate, color=cumulative_rate)) +
  geom_tile() + scale_fill_viridis() + 
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1)) + 
  labs(x = "Tag distance", y = "Permissiveness", fill= "Acceptance\nrate") +
  theme_classic() + guides(color="none") + 
  scale_color_viridis()
