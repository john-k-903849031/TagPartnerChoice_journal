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

td_colors <- rev(c("#3E9BFEFF", "#46F884FF", "#E1DD37FF", "#F05B12FF", "#7A0403FF"))

########### TAG DISTANCE PLOTS ###########

plot_fixed_data <- function(mut_folder, para_folder, param_facet_title="Vertical transmission rate"){
  a <- read_csv(paste0("../../data/",mut_folder,"org_dump.dat"),col_names=T)
  a$sym_tag[a$sym_tag==""]<-NA
  a$cond<-"mut"
  print("Done reading mut!")
  
  b <- read_csv(paste0("../../data/",para_folder,"org_dump.dat"),col_names=T)
  b$sym_tag[b$sym_tag==""]<-NA
  b$cond<-"para"
  print("Done reading para!")
  
  small <- rbind(b %>% filter(!is.na(sym_tag) & tag_perm != 0.125 & tag_perm != 0.1875 & tag_perm != 0.25), 
                 a %>% filter(!is.na(sym_tag) & tag_perm != 0.125 & tag_perm != 0.1875 & tag_perm != 0.25))
  rm(b, a)
  
  small$tag_perm <- factor(small$tag_perm, levels=c("0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625"))
  
  if(param_facet_title == "Tag mutation rate") {
    small <- small %>% filter(tag_mut > 0) %>% rename(param = tag_mut)
  } else {
    small <- small %>% rename(param = vt)
  }
  small <- small %>% 
    mutate(sym_bin = case_when(sym_int < -0.6 ~ "-1 to -0.6 (Parasitic)", 
                               sym_int <= 1 & sym_int >= 0.6 ~ "0.6 to 1.0 (Mutualistic)",
                               .default = "error"),
           param = as.factor(param)) %>%
    group_by(sym_bin,tag_perm,param) %>% mutate(td_mean =  mean(tag_distance)) 
  print("Done data manipulating!")
  
  plot <- ggplot(small, aes(x=tag_distance)) + 
    geom_histogram(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill="Fixed parasites",y=..count../10000),alpha=0.7,  bins=16)+
    geom_histogram(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill="Fixed mutualists",y=..count../10000), alpha=0.7, bins=16) +
    geom_vline(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(xintercept = td_mean), color=td_colors[5],linetype="dotted") +
    geom_vline(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(xintercept = td_mean), color=td_colors[1],linetype="dotted") +
    xlab("Tag proportional mismatch") + 
    scale_fill_manual(name="Symbiont behavior",
                      values=c("Fixed parasites"=td_colors[1], "Fixed mutualists"=td_colors[5]),
                      breaks=c("Fixed parasites", "Fixed mutualists")) + 
    ylab("Pooled count across replicates / 10,000")+
    #facet_nested("a" + tag_perm ~ param_facet_title + param, drop=TRUE) + 
    facet_nested_theme + theme(legend.position = 'bottom') + 
    theme(strip.text.y.right = element_text(angle = 0)) +
    scale_y_continuous(breaks = c(0,4,8))+
    scale_x_continuous(breaks = c(0,0.4,0.8))
  rm(small)
  
  print("Returning plot!")
  return(plot)
}

#### tag mut control #### 
tm_fix_plot <- plot_fixed_data(mut_folder = "exp_4_base_tagmut_mut_fixed/", 
                               para_folder = "exp_4_base_tagmut_para_fixed/", 
                               param_facet_title = "Tag mutation rate")

effect_sizes <- read_csv("../stats/tm_fix_deltas.csv") %>%
  mutate(annotation = case_when(abs(deltas) < 0.11 ~ "x", 
                                abs(deltas) < 0.28 ~ "+", 
                                abs(deltas) < 0.43 ~ "++", 
                                .default = "+++")) %>%
  separate_wider_delim(col = names, "_tm", names =  c("tag_perm", "param")) %>% 
  separate_wider_delim(col = tag_perm, "tp", names =  c(NA, "tag_perm")) %>%
  mutate(tag_perm = as.factor(tag_perm), param = as.factor(param))  %>%
  filter(param != 0)

combo_tm_fix_plot <- tm_fix_plot + 
  facet_nested("Tag\npermiss-\niveness" + tag_perm ~ "Tag mutation rate" + param, drop=TRUE) + 
  theme(strip.text.y.right = element_text(angle = 0)) +
  geom_text(data = effect_sizes, na.rm = T,
                                     mapping = aes(x = Inf,
                                                   y = Inf,
                                                   vjust=1,
                                                   hjust=1,
                                                   label = annotation), size=3)  

ggplot2::ggsave("parastab_base_control_distances.pdf", combo_tm_fix_plot, width = 6, height = 4, units = "in")


########### vt control ###########

vt_fix_plot <- plot_fixed_data("exp_3_base_vt_mut_fixed/","exp_3_base_vt_para_fixed/", "Vertical transmission rate")  ########## SHOULD BE AT THIS POINT

effect_sizes <- read_csv("../stats/vt_fix_deltas.csv") %>%
  mutate(annotation = case_when(abs(deltas) < 0.11 ~ "x", 
                                abs(deltas) < 0.28 ~ "+", 
                                abs(deltas) < 0.43 ~ "++", 
                                .default = "+++")) %>%
  separate_wider_delim(col = names, "_vt", names =  c("tag_perm", "param")) %>% 
  separate_wider_delim(col = tag_perm, "tp", names =  c(NA, "tag_perm")) %>%
  mutate(tag_perm = as.factor(tag_perm), param = as.factor(param))

combo_vt_fix_plot <- vt_fix_plot +
  facet_nested("Tag\npermiss-\niveness" + tag_perm ~ "Vertical transmission rate" + param, drop=TRUE) + 
  theme(strip.text.y.right = element_text(angle = 0)) +
  geom_text(data = effect_sizes, na.rm = T,
                        mapping = aes(x = Inf,
                                      y = Inf,
                                      vjust=1,
                                      hjust=1,
                                      label = annotation), size=3)  
combo_vt_fix_plot

ggplot2::ggsave("vt_base_control_distances.pdf", combo_vt_fix_plot, width = 6, height = 4, units = "in")

plot_evo_data <- function(folder, cond){
  small <- read_csv(paste0("../../data/",folder,"org_dump.dat"),col_names=T)
  small$sym_tag[small$sym_tag==""]<-NA
  if(cond == "tm") {
    small <- small %>% filter(tag_mut > 0) %>% rename(param = tag_mut)
  } else {
    small <- small %>% filter(vt <= 0.5) %>% rename(param = vt)
  }
  
  small <- small %>% filter(!is.na(sym_tag) & tag_perm != "-1" & tag_perm != "none") %>% 
                              mutate(tag_distance = as.numeric(tag_distance),
                                     tag_perm = factor(tag_perm, levels=c("0.125",  "0.1875", "0.25",  
                                                                        "0.3125", "0.375",  "0.4375", "0.5",    "0.5625", "0.625")),
                              sym_bin = case_when(sym_int < -0.6 ~ "-1 to -0.6 (Parasitic)",
                                                  sym_int < -0.2 & sym_int >= -0.6 ~ "-0.6 to -0.2 (Detrimental)",
                                                  sym_int < 0.2  & sym_int >= -0.2 ~ "-0.2 to 0.2 (Nearly neutral)",
                                                  sym_int < 0.6  & sym_int >= 0.2  ~ "0.2 to 0.6 (Positive)",
                                                  sym_int <= 1   & sym_int >= 0.6  ~ "0.6 to 1.0 (Mutualistic)",
                                                  .default = "error"
                                                  ))

  
  plot <- ggplot(small, aes(x=tag_distance)) + 
    geom_histogram(data=subset(small, sym_bin=="-1 to -0.6 (Parasitic)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16)+
    geom_histogram(data=subset(small, sym_bin=="-0.6 to -0.2 (Detrimental)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
    geom_histogram(data=subset(small, sym_bin=="-0.2 to 0.2 (Nearly neutral)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
    geom_histogram(data=subset(small, sym_bin=="0.2 to 0.6 (Positive)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
    geom_histogram(data=subset(small, sym_bin=="0.6 to 1.0 (Mutualistic)"), aes(fill=sym_bin,y=..count../10000), alpha=0.7, bins=16) +
    xlab("Tag proportional mismatch") + 
    scale_fill_manual(name="Symbiont\nbehavior",
                      values=c("-1 to -0.6 (Parasitic)"=td_colors[1], '-0.6 to -0.2 (Detrimental)'=td_colors[2], 
                               '-0.2 to 0.2 (Nearly neutral)'=td_colors[3],
                               "0.2 to 0.6 (Positive)"=td_colors[4], "0.6 to 1.0 (Mutualistic)"=td_colors[5]),
                      breaks=c("-1 to -0.6 (Parasitic)",'-0.6 to -0.2 (Detrimental)',
                               '-0.2 to 0.2 (Nearly neutral)',
                               "0.2 to 0.6 (Positive)", "0.6 to 1.0 (Mutualistic)")) + 
    ylab("Pooled count across replicates / 10,000")+
    #facet_nested("Tag permissiveness" + tag_perm ~ "Tag mutation rate" + tag_mut, drop=TRUE) + 
    scale_y_continuous(breaks = c(0,5)) + 
    facet_nested_theme+ theme(legend.position = 'bottom')+
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(legend.key.size = unit(0.4, 'cm'), 
          legend.title = element_text(size=8), #change legend title font size
          legend.text = element_text(size=8)) #change legend text font size
  plot
}


########### parastab evo ########### 

tagmut_evo_plot <- plot_evo_data("exp_2_base_tagmut_evolving/", "tm") 

effect_sizes <- read_csv("../stats/tagmut_evolve_paravsmut_deltas.csv") %>%
  mutate(annotation = case_when(abs(deltas) < 0.11 ~ "x", 
                                abs(deltas) < 0.28 ~ "+", 
                                abs(deltas) < 0.43 ~ "++", 
                                .default = "+++")) %>%
  separate_wider_delim(col = names, "_tm", names =  c("tag_perm", "param")) %>% 
  separate_wider_delim(col = tag_perm, "tp", names =  c(NA, "tag_perm")) %>%
  mutate(tag_perm = as.factor(tag_perm), param = as.factor(param))

combo_tm_evo_plot <- tagmut_evo_plot + 
  facet_nested("Tag\npermiss-\niveness" + tag_perm ~ "Tag mutation rate" + param, drop=TRUE) + 
  theme(strip.text.y.right = element_text(angle = 0)) +
  geom_text(data = effect_sizes, na.rm = T,
            mapping = aes(x = Inf,
                          y = Inf,
                          vjust=1,
                          hjust=1,
                          label = annotation), size=3)  

ggplot2::ggsave("parastab_base_evo_distances.pdf", combo_tm_evo_plot, width = 6, height = 4, units = "in")



########### vt evo ########### 
vt_evo_plot <- plot_evo_data("exp_1_base_vt_evolving/", "vt") 


effect_sizes <- read_csv("../stats/vt_evolve_paravsmut_deltas.csv") %>%
  select(tag_perm:deltas) %>%
  mutate(tag_perm = as.factor(tag_perm), param = as.double(param)) %>%
  mutate(annotation = case_when(abs(deltas) < 0.11 ~ "x", 
                                abs(deltas) < 0.28 ~ "+", 
                                abs(deltas) < 0.43 ~ "++", 
                                .default = "+++")) 

combo_vt_evo_plot <- vt_evo_plot + 
  facet_nested("Tag\npermiss-\niveness" + tag_perm ~ "Vertical transmission rate" + param, drop=TRUE) + 
  theme(strip.text.y.right = element_text(angle = 0)) +
  geom_text(data = effect_sizes, na.rm = T,
                                      mapping = aes(x = Inf,
                                                    y = Inf,
                                                    vjust=1,
                                                    hjust=1,
                                                    label = annotation), size=3)  
  
ggplot2::ggsave("vt_base_evo_distances.pdf", combo_vt_evo_plot, width = 6, height = 4, units = "in")
