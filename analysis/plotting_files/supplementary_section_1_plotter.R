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

tenhelix <- rev(c("#891901", "#B50142", "#D506AD", "#AB08FF", "#5731FD", "#4755FF", "#42bcf5", "#42e0f5", "#9af1fc", "#c5f6fc"))


base_parastab_behavior_plot <- function(){
  ## SYMBIONT
  # build dataframe
  hist_df <- syms_df %>% 
    group_by(tag_perm, tag_mut) %>% 
    mutate(com_rep = match(seed, unique(seed))) %>% 
    select(c(com_rep, tag_perm, tag_mut), `Hist_-1`:Hist_0.9) %>%
    pivot_longer(-c(com_rep, tag_perm, tag_mut)) %>%
    mutate(variable = case_when(name == 'Hist_-1'   | name == 'Hist_-0.9' ~ "-1 to -0.8 (Parasitic)", 
                                name == 'Hist_-0.8' | name == 'Hist_-0.7' ~ "-0.8 to -0.6 (Parasitic)", 
                                name == 'Hist_-0.6' | name == 'Hist_-0.5' ~ "-0.6 to -0.4 (Detrimental)", 
                                name == 'Hist_-0.4' | name == 'Hist_-0.3' ~ "-0.4 to -0.2 (Detrimental)", 
                                name == 'Hist_-0.2' | name == 'Hist_-0.1' ~ "-0.2 to 0 (Nearly Neutral)", 
                                name == 'Hist_0.0'  | name == 'Hist_0.1'  ~ "0 to 0.2 (Nearly Neutral)", 
                                name == 'Hist_0.2'  | name == 'Hist_0.3'  ~ "0.2 to 0.4 (Positive)", 
                                name == 'Hist_0.4'  | name == 'Hist_0.5'  ~ "0.4 to 0.6 (Positive)", 
                                name == 'Hist_0.6'  | name == 'Hist_0.7'  ~ "0.6 to 0.8 (Mutualistic)", 
                                name == 'Hist_0.8'  | name == 'Hist_0.9'  ~ "0.8 to 1.0 (Mutualistic)",
                                .default = "ERROR"))
  
  # force set the order of the variable column
  hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))
  
  # compress bin counts
  hist_df <- hist_df %>%
    group_by(com_rep, tag_perm, tag_mut, variable) %>% 
    summarise(value = sum(value))
  
  hist_df$tag_perm <- factor(hist_df$tag_perm, levels=rev(unique(hist_df$tag_perm)))
  
  # order by proportion of ubermutualists
  hist_df$max_mut_count <- rep(hist_df$value[hist_df$variable=="0.8 to 1.0 (Mutualistic)"],1,each=10)
  hist_df <- hist_df %>% 
    group_by(tag_perm, tag_mut,com_rep) %>% 
    mutate(max_mut_prop = max_mut_count/sum(value)) %>% 
    group_by(tag_perm, tag_mut) %>% 
    arrange(desc(max_mut_prop)) %>% mutate(group_id = rep(1:30,1,each=10))
  
  a <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value)) + 
    geom_col(position="fill", na.rm = TRUE, width = 1) +
    scale_fill_manual(name="Behavior",values=tenhelix,guide="none")+
    xlab("Replicate") + ylab("Proportion of symbiont population") +
    facet_nested("Tag permissiveness" + tag_perm ~ "Tag mutation rate" + tag_mut) + 
    facet_nested_theme
  
  syms_hist <- hist_df
  
  ## HOSTS
  hist_df <- hosts_df %>% 
    group_by(tag_perm, tag_mut) %>% 
    mutate(com_rep = match(seed, unique(seed))) %>% 
    select(c(com_rep, tag_perm, tag_mut), `Hist_-1`:Hist_0.9) %>%
    pivot_longer(-c(com_rep, tag_perm, tag_mut)) %>%
    mutate(variable = case_when(name == 'Hist_-1'   | name == 'Hist_-0.9' ~ "-1 to -0.8 (Parasitic)", 
                                name == 'Hist_-0.8' | name == 'Hist_-0.7' ~ "-0.8 to -0.6 (Parasitic)", 
                                name == 'Hist_-0.6' | name == 'Hist_-0.5' ~ "-0.6 to -0.4 (Detrimental)", 
                                name == 'Hist_-0.4' | name == 'Hist_-0.3' ~ "-0.4 to -0.2 (Detrimental)", 
                                name == 'Hist_-0.2' | name == 'Hist_-0.1' ~ "-0.2 to 0 (Nearly Neutral)", 
                                name == 'Hist_0.0'  | name == 'Hist_0.1'  ~ "0 to 0.2 (Nearly Neutral)", 
                                name == 'Hist_0.2'  | name == 'Hist_0.3'  ~ "0.2 to 0.4 (Positive)", 
                                name == 'Hist_0.4'  | name == 'Hist_0.5'  ~ "0.4 to 0.6 (Positive)", 
                                name == 'Hist_0.6'  | name == 'Hist_0.7'  ~ "0.6 to 0.8 (Mutualistic)", 
                                name == 'Hist_0.8'  | name == 'Hist_0.9'  ~ "0.8 to 1.0 (Mutualistic)",
                                .default = "ERROR"))
  
  # force set the order of the variable column
  hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))
  
  # compress bin counts
  hist_df <- hist_df %>%
    group_by(com_rep, tag_perm, tag_mut, variable) %>% 
    summarise(value = sum(value))
  
  hist_df$tag_perm <- factor(hist_df$tag_perm, levels=rev(unique(hist_df$tag_perm)))
  
  hist_df <- hist_df %>% left_join(select(syms_hist, c(tag_perm, tag_mut, group_id, variable, com_rep)), 
                                   by= join_by(tag_perm,variable,tag_mut,com_rep))
  
  b <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value)) + 
    geom_col(position="fill", na.rm=T, width = 1) +
    scale_fill_manual(name="Behavior",values=tenhelix)+
    xlab("Replicate") + ylab("Proportion of host population") +
    facet_nested("Tag permissiveness" + tag_perm ~ "Tag mutation rate" + tag_mut) + 
    facet_nested_theme
  
  legend <- cowplot::get_legend(b+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
  plots <- grid.arrange(a+ theme(axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank()),
                        b+theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank()),
                        ncol=2,widths=c(1,1))
  grid.arrange(plots,legend, ncol=1, nrow=2,heights=c(1,0.2))
  
  rm(hist_df, a, b, syms_hist)
} 

base_vt_behavior_plot <- function(){
  ## SYMBIONT
  # build dataframe
  hist_df <- syms_df %>% 
    group_by(tag_perm, vt) %>% 
    mutate(com_rep = match(seed, unique(seed))) %>% 
    select(c(com_rep, tag_perm, vt), `Hist_-1`:Hist_0.9) %>%
    pivot_longer(-c(com_rep, tag_perm, vt)) %>%
    mutate(variable = case_when(name == 'Hist_-1'   | name == 'Hist_-0.9' ~ "-1 to -0.8 (Parasitic)", 
                                name == 'Hist_-0.8' | name == 'Hist_-0.7' ~ "-0.8 to -0.6 (Parasitic)", 
                                name == 'Hist_-0.6' | name == 'Hist_-0.5' ~ "-0.6 to -0.4 (Detrimental)", 
                                name == 'Hist_-0.4' | name == 'Hist_-0.3' ~ "-0.4 to -0.2 (Detrimental)", 
                                name == 'Hist_-0.2' | name == 'Hist_-0.1' ~ "-0.2 to 0 (Nearly Neutral)", 
                                name == 'Hist_0.0'  | name == 'Hist_0.1'  ~ "0 to 0.2 (Nearly Neutral)", 
                                name == 'Hist_0.2'  | name == 'Hist_0.3'  ~ "0.2 to 0.4 (Positive)", 
                                name == 'Hist_0.4'  | name == 'Hist_0.5'  ~ "0.4 to 0.6 (Positive)", 
                                name == 'Hist_0.6'  | name == 'Hist_0.7'  ~ "0.6 to 0.8 (Mutualistic)", 
                                name == 'Hist_0.8'  | name == 'Hist_0.9'  ~ "0.8 to 1.0 (Mutualistic)",
                                .default = "ERROR"))
  
  # force set the order of the variable column
  hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))
  
  # compress bin counts
  hist_df <- hist_df %>%
    group_by(com_rep, tag_perm, vt, variable) %>% 
    summarise(value = sum(value))
  
  hist_df$tag_perm <- factor(hist_df$tag_perm, levels=rev(unique(hist_df$tag_perm)))
  
  # order by proportion of ubermutualists
  hist_df$max_mut_count <- rep(hist_df$value[hist_df$variable=="0.8 to 1.0 (Mutualistic)"],1,each=10)
  hist_df <- hist_df %>% 
    group_by(tag_perm, vt,com_rep) %>% 
    mutate(max_mut_prop = max_mut_count/sum(value)) %>% 
    group_by(tag_perm, vt) %>% 
    arrange(desc(max_mut_prop)) %>% mutate(group_id = rep(1:30,1,each=10))
  
  a <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value)) + 
    geom_col(position="fill", na.rm = TRUE, width = 1) +
    scale_fill_manual(name="Behavior",values=tenhelix,guide="none")+
    xlab("Replicate") + ylab("Proportion of symbiont population") +
    facet_nested("Tag permissiveness" + tag_perm ~ "Vertical transmission rate" + vt) + 
    facet_nested_theme
  
  syms_hist <- hist_df
  
  ## HOSTS
  hist_df <- hosts_df %>% 
    group_by(tag_perm, vt) %>% 
    mutate(com_rep = match(seed, unique(seed))) %>% 
    select(c(com_rep, tag_perm, vt), `Hist_-1`:Hist_0.9) %>%
    pivot_longer(-c(com_rep, tag_perm, vt)) %>%
    mutate(variable = case_when(name == 'Hist_-1'   | name == 'Hist_-0.9' ~ "-1 to -0.8 (Parasitic)", 
                                name == 'Hist_-0.8' | name == 'Hist_-0.7' ~ "-0.8 to -0.6 (Parasitic)", 
                                name == 'Hist_-0.6' | name == 'Hist_-0.5' ~ "-0.6 to -0.4 (Detrimental)", 
                                name == 'Hist_-0.4' | name == 'Hist_-0.3' ~ "-0.4 to -0.2 (Detrimental)", 
                                name == 'Hist_-0.2' | name == 'Hist_-0.1' ~ "-0.2 to 0 (Nearly Neutral)", 
                                name == 'Hist_0.0'  | name == 'Hist_0.1'  ~ "0 to 0.2 (Nearly Neutral)", 
                                name == 'Hist_0.2'  | name == 'Hist_0.3'  ~ "0.2 to 0.4 (Positive)", 
                                name == 'Hist_0.4'  | name == 'Hist_0.5'  ~ "0.4 to 0.6 (Positive)", 
                                name == 'Hist_0.6'  | name == 'Hist_0.7'  ~ "0.6 to 0.8 (Mutualistic)", 
                                name == 'Hist_0.8'  | name == 'Hist_0.9'  ~ "0.8 to 1.0 (Mutualistic)",
                                .default = "ERROR"))
  
  # force set the order of the variable column
  hist_df$variable <- factor(hist_df$variable, levels=rev(unique(hist_df$variable)))
  
  # compress bin counts
  hist_df <- hist_df %>%
    group_by(com_rep, tag_perm, vt, variable) %>% 
    summarise(value = sum(value))
  
  hist_df$tag_perm <- factor(hist_df$tag_perm, levels=rev(unique(hist_df$tag_perm)))
  
  hist_df <- hist_df %>% left_join(select(syms_hist, c(tag_perm, vt, group_id, variable, com_rep)), 
                                   by= join_by(tag_perm,variable,vt,com_rep))
  
  b <- ggplot(hist_df, aes(fill=variable, x=group_id,y=value)) + 
    geom_col(position="fill", na.rm=T, width = 1) +
    scale_fill_manual(name="Behavior",values=tenhelix)+
    xlab("Replicate") + ylab("Proportion of host population") +
    facet_nested("Tag permissiveness" + tag_perm ~ "Vertical transmission rate" + vt) + 
    facet_nested_theme
  
  legend <- cowplot::get_legend(b+ theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
  plots <- grid.arrange(a+ theme(axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank()),
                        b+theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank()),
                        ncol=2,widths=c(1,1))
  grid.arrange(plots,legend, ncol=1, nrow=2,heights=c(1,0.2))
  
  rm(hist_df, a, b, syms_hist)
} 



# Base, parastab
folder <- "../../data/exp_2_base_tagmut_evolving/"
syms_df <- read_csv(paste0(folder,"sym_counts.dat"))
hosts_df <- read_csv(paste0(folder,"host_counts.dat"))
syms_df <- subset(syms_df, update==max(update))
hosts_df <- subset(hosts_df, update==max(update))
hosts_df <- hosts_df %>% mutate(tag_perm = case_when(tag_perm == "none" ~ "no tag", .default=tag_perm),
                                tag_mut = case_when(tag_mut == "none" ~ "no tag", .default=tag_mut))
syms_df <- syms_df %>% mutate(tag_perm = case_when(tag_perm == "none" ~ "no tag", .default=tag_perm),
                              tag_mut = case_when(tag_mut == "none" ~ "no tag", .default=tag_mut))
hosts_df$tag_perm <- factor(hosts_df$tag_perm, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))
syms_df$tag_perm <- factor(syms_df$tag_perm, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))
hosts_df$tag_mut <- factor(hosts_df$tag_mut, levels=c("no tag", "0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1"))
syms_df$tag_mut <- factor(syms_df$tag_mut, levels=c("no tag", "0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1"))

#6x12in
base_parastab_behavior_plot() 
small <-  syms_df %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                                   mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                            ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                            ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                           .default = "Coexistence"))
small %>% filter(tag_perm == "no tag") %>% group_by(class) %>% summarize(val = n())

# Base, vtsweep
folder <- "../../data/exp_1_base_vt_evolving/"
syms_df <- read_csv(paste0(folder,"sym_counts.dat"))
hosts_df <- read_csv(paste0(folder,"host_counts.dat"))
syms_df <- subset(syms_df, update==max(update))
hosts_df <- subset(hosts_df, update==max(update))
hosts_df <- hosts_df %>% mutate(tag_perm = case_when(tag_perm == "none" ~ "no tag", .default=tag_perm))
syms_df <- syms_df %>% mutate(tag_perm = case_when(tag_perm == "none" ~ "no tag", .default=tag_perm))
hosts_df$tag_perm <- factor(hosts_df$tag_perm, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))
syms_df$tag_perm <- factor(syms_df$tag_perm, levels=c("0.125",  "0.1875", "0.25",   "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625" ,"no tag"))

#6x12in
base_vt_behavior_plot() 
small <-  syms_df %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                                mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                            ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                            ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                            .default = "Coexistence"))
small %>% filter(tag_perm == "no tag") %>% group_by(class, vt) %>% summarize(val = n())
