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
  labs(x = "Host-symbiont tag distance", y = "Host permissiveness", fill= "Symbiont\nacceptance\nrate") +
  theme_classic() + guides(color="none") + 
  scale_color_viridis() +
  scale_y_discrete(labels=c('0.125'=expression(bold("0.125")), '0.1875'=expression(bold("0.1875")), 
                            '0.25'=expression(bold("0.25")), '0.3125'=expression(bold("0.3125")), 
                            '0.375'=expression(bold("0.375")), '0.4375'=expression(bold("0.4375")), 
                            '0.5'=expression(bold("0.5")), '0.5625'=expression(bold("0.5625")), 
                            '0.625'=expression(bold("0.625")), parse=TRUE)) + 
  geom_hline(aes(yintercept = 1.5), linetype= "dashed", color = "white") + 
  geom_hline(aes(yintercept = 19.5), linetype= "dashed", color = "white") 
