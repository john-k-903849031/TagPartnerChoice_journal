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

## themes and colors ##
options(scipen=999)
facet_nested_theme <- theme(
  strip.background = element_rect(fill = "white", colour = "grey", linetype="dotted", linewidth=0.2), 
  panel.background = element_rect(fill='white', color='grey'),
  panel.border = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.key = element_rect(color = "transparent"),
  panel.spacing.x = unit(0.5,"line"),
  panel.spacing.y = unit(0.5,"line"),
  legend.direction = "horizontal",legend.justification="left" ,
  axis.text.x = element_blank(),# element_text(angle = 12.5, hjust = 0.5, vjust=0.5),
  legend.box.just = "bottom", legend.position = "bottom")
  #legend.margin=margin(0,0,0,0),
  #legend.box.margin=margin(-10,-10,-10,-10)) 

perm_color_map <- c("x <= 0" = "#440154FF",
                    "0 < x <= 0.25" = "#414487FF", 
                    "0.25 < x <= 0.5" = "#2A788EFF", 
                    "0.5 < x <= 0.75" = "#22A884FF", 
                    "0.75 < x <= 1" = "#7AD151FF",             
                    "x > 1" = "#FDE725FF") 
behavior_color_map <- c("Exclusive\nparasitism"="#891901",
                        "Exclusive\nmutualism"="#c5f6fc",
                        "Coexistence"="#66ff99",
                        "Symbiont\nextinction"="#cc9900")
## classification methods ##
classify_behavior <- function(df){
  small <- df %>% group_by(seed) %>% mutate(para_count = rowSums(across(`Hist_-1`:`Hist_-0.3`)), 
                                            mut_count = rowSums(across(Hist_0.2:Hist_0.9)))
  small <- small %>% mutate(class = case_when(count == 0 ~ "Symbiont extinction",
                                              ((para_count > 0) & (count - para_count) == 0) ~ "Exclusive parasitism",
                                              ((mut_count > 0) & (count - mut_count) == 0) ~ "Exclusive mutualism",
                                              .default = "Coexistence"))
  return(small)
}
classify_permissivenss <- function(df){
  small <- df %>% group_by(seed,evo_mut,param) %>%
    filter(symbiont_tag_richness > 0) %>%
    group_by(seed,evo_mut,param) %>%
    filter(update==max(update))
  small <- small %>% group_by(seed) %>% mutate(class = case_when(mean_host_perm <= 0 ~ "x <= 0",
                                                                 mean_host_perm <= 0.25 & mean_host_perm > 0 ~ "0 < x <= 0.25",
                                                                 mean_host_perm <= 0.5 & mean_host_perm > 0.25 ~ "0.25 < x <= 0.5",
                                                                 mean_host_perm <= 0.75 & mean_host_perm > 0.5 ~ "0.5 < x <= 0.75",
                                                                 mean_host_perm <= 1 & mean_host_perm > 0.75 ~ "0.75 < x <= 1",
                                                                 mean_host_perm > 1 ~ "x > 1", .default = "ERROR"))
  small$class <- factor(small$class, levels=c("x <= 0", 
                                              "0 < x <= 0.25", 
                                              "0.25 < x <= 0.5", 
                                              "0.5 < x <= 0.75", 
                                              "0.75 < x <= 1",             
                                              "x > 1" ))
  return(small)
}

####### interaction values ####### 
## parastab ##
folder <- "../../data/exp_6_evoperm_tagmut/"
para_syms_df <- read_csv(paste0(folder,"sym_counts.dat")) %>% 
  filter(update == max(update)) %>% rename(param = tag_mut) %>% mutate(param_name = "tag_mut")
para_syms_df <- classify_behavior(para_syms_df)

## combo ## 
combo_syms_df <- rbind(para_syms_df) %>% select(param_name, param, evo_mut, seed, class) %>% rename(behavior_class = class)

####### tags ####### 
## parastab ## 
folder <- "../../data/exp_6_evoperm_tagmut/"
para_tags_df <- read_csv(paste0(folder,"tag_dists.dat")) %>% 
  rename(param = tag_mut) %>% mutate(param_name = "tag_mut")
para_tags_df <- classify_permissivenss(para_tags_df)

## combo ## 
combo_tags_df <- rbind(para_tags_df) %>% select(param_name, param, evo_mut, seed, class)

####### interaction values pre-extinction ####### 
lsu_para_syms_df <- classify_behavior(
  read_csv(paste0("../../data/exp_6_evoperm_tagmut/","sym_counts.dat")) %>% rename(param = tag_mut) %>%
  group_by(seed,evo_mut,param) %>% filter(count > 0) %>% filter(update==max(update)) %>% 
  left_join(para_tags_df %>% rename(perm_class = class))) %>%
  select(param_name, param, evo_mut, perm_class, class, update) %>%
  rename(pre_extinction_update = update, 
         pre_extinction_behavior_class = class,
         class = perm_class) 

##### combined behavior and tags ##### 
combo_df <- combo_syms_df %>% group_by(param_name, param, evo_mut, seed) %>% 
  left_join(combo_tags_df) %>%
  mutate(evo_mut = as.factor(evo_mut)) %>% left_join(lsu_para_syms_df%>% mutate(evo_mut = as.factor(evo_mut))) %>%
  group_by(param_name, param, evo_mut, class, behavior_class, pre_extinction_behavior_class) %>%
  summarize(row_count = n())%>%
  mutate(behavior_class = str_wrap(behavior_class, width = 10),
         pre_extinction_behavior_class = str_wrap(pre_extinction_behavior_class, width = 10))


subset_combo_df <- combo_df %>%
  filter(param_name == "tag_mut" & param >= 0.0001 & param <= 0.001) %>%
  mutate(evo_mut = as.numeric(as.character(evo_mut)),
         category = case_when(evo_mut <= 0.001 ~ "lower", 
                              evo_mut > 0.001 ~ "higher",
                              .default = "other"),
         param_label = case_when(param_name == "tag_mut" ~ "0.01% <= TM rate <= 0.1%",
                                 .default="other")) %>%
  filter(category != "other") 

View(
  subset_combo_df %>% 
    group_by(category) %>% mutate(category_count = sum(row_count)) %>%
    group_by(category, behavior_class,category_count) %>% 
    summarize(behavior_count = sum(row_count)) %>%
    mutate(percent = round(behavior_count/category_count, digits = 3)*100)
)

View(
  subset_combo_df %>% 
    filter((behavior_class == "Symbiont\nextinction") |
             (category=="lower"))%>%
    group_by(category) %>% mutate(category_count = sum(row_count)) %>%
    group_by(category, class, behavior_class, category_count) %>% 
    summarize(perm_count = sum(row_count)) %>%
    mutate(percent = round(perm_count/category_count, digits = 3)*100)
)

View(
  subset_combo_df %>% 
    filter((behavior_class == "Symbiont\nextinction"))%>%
    group_by(behavior_class) %>% mutate(category_count = sum(row_count)) %>%
    group_by(pre_extinction_behavior_class, behavior_class, category_count) %>% 
    summarize(perm_count = sum(row_count)) %>%
    mutate(percent = round(perm_count/category_count, digits = 3)*100)
)


##### plots #####
head(subset_combo_df)
a <- subset_combo_df %>% 
  ungroup()%>%
  add_row(param_name = "other", 
          param = 0.0001,
          evo_mut = 0,
          class = "0.25 < x <= 0.5",            
          behavior_class = "Exclusive\nparasitism",       
          pre_extinction_behavior_class = "Exclusive\nparasitism",                          
          row_count = 0,
          category = "lower",
          param_label = "0.01% <= TM rate <= 0.1%")%>%
  mutate(facet_label = case_when(category == "lower"~"Lower\npermissiveness\nmutation",
                                 category == "higher"~"Higher\npermissiveness\nmutation",
                                 .default = "other"),
         facet_label = factor(facet_label, levels = c("Lower\npermissiveness\nmutation", "Higher\npermissiveness\nmutation", "other"))) %>%
  group_by(facet_label, param_label, behavior_class) %>% summarize(total = sum(row_count))%>%
  ggplot(aes(x=param_label, fill = behavior_class, y=total)) +
  geom_col() +
  facet_wrap( ~ as.factor(facet_label)) + 
  facet_nested_theme +
  scale_fill_manual(name="Behavior",values=behavior_color_map) + 
  labs(x = "", y = "Replicate count") + 
  guides(fill = guide_legend(nrow = 2))
a

b <- subset_combo_df %>% 
  mutate(facet_label = case_when((category == "lower" & (behavior_class == "Coexistence" | behavior_class == "Exclusive\nmutualism"))~"Lower\npermissiveness\nmutation\nand coexistence or\nexclusive mutualism",
                                 (category == "higher" & behavior_class == "Symbiont\nextinction")~"Higher\npermissiveness\nmutation and\nsymbiont extinction",
                                 .default = "other"),
         facet_label = factor(facet_label, levels = c("Lower\npermissiveness\nmutation\nand coexistence or\nexclusive mutualism", "Higher\npermissiveness\nmutation and\nsymbiont extinction", "other"))) %>%
  filter(facet_label!= "other") %>%
  group_by(facet_label, param_label, class) %>% summarize(total = sum(row_count))%>%
  ggplot(aes(x = param_label, y = total, fill = class)) +
  geom_col() + 
  scale_fill_manual(name="Mean final\npermissiveness", values = perm_color_map) + 
  labs(x = "", 
       y = "Replicate count") +
  facet_wrap(~as.factor(facet_label)) +
  facet_nested_theme + 
  guides(fill = guide_legend(nrow = 2))

a / b + plot_annotation(tag_levels  = list(c("a)", "b)")))


##### pre-extinction behavior #####


c <- subset_combo_df %>% 
  filter(behavior_class == "Symbiont\nextinction") %>%
  mutate(facet_label = category) %>%
  group_by(facet_label, param_label, pre_extinction_behavior_class) %>% summarize(total = sum(row_count))%>%
  ggplot(aes(x = param_label, y = total, fill = as.factor(pre_extinction_behavior_class))) +
  geom_col() + 
  scale_fill_manual(name="Mean final\npermissiveness", values = behavior_color_map) + 
  labs(x = "", 
       y = "Replicate count") +
  facet_wrap(~"Pre-extinction behavior") +
  facet_nested_theme + 
  guides(fill = guide_legend(nrow = 2))
c


#legend.margin=margin(0,0,0,0),
#legend.box.margin=margin(-10,-10,-10,-10)
?margin
#t = 0, r = 0, b = 0, l = 0
## combined plot 
legend_a <- cowplot::get_legend(a + theme(legend.title.position = "top",
                                          legend.margin=margin(-10, l = 10),
                                          legend.box.margin=margin(-50)
                                          ) + 
                                  guides(fill = guide_legend(nrow = 2,
                                                             label.theme = element_text(size = 8)))) 
legend_b <- cowplot::get_legend(b + theme(legend.title.position = "top",
                                          legend.margin=margin(-10, l = 10),
                                          legend.box.margin=margin(-50)
                                          )+ 
                                  guides(fill = guide_legend(nrow = 3,
                                                             label.theme = element_text(size = 8))))
plot <- wrap_plots(a+ theme(legend.position = "none"),
                   c + theme(legend.position = "none"),
                   b+ theme(legend.position = "none"), 
                   wrap_plots(A=legend_a,  B=legend_b, design="A
                                                               B"), 
                   ncol = 2,
                   nrow = 2,
                   widths = c(1.5,1),
                   heights = c(1, 1)) +
  plot_annotation(tag_levels  = list(c("a)", "c)", "b)")))
plot
