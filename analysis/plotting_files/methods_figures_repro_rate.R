library(patchwork)
#### themes ####
facet_nested_theme <- theme(
  strip.background = element_rect(fill = "white", colour = "grey", linetype="dotted", size = 0.2, linewidth=30), 
  panel.background = element_rect(fill='white', color='grey'),
  panel.border = element_blank(),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.key = element_rect(color = "transparent"),
  panel.spacing.x = unit(0,"line"),
  panel.spacing.y = unit(0,"line"))

#### functions #### 
host_res_f <- function(host_int, sym_int){
  host_res <- 0
  host_donate <- 0
  if(host_int < 0){
    host_res <- res_d - (host_int * res_d * -1)
  } 
  else{
    host_donate <- host_int * res_d
    host_res <- res_d - (host_int * res_d)
  }
  
  if(sym_int < 0){
    if(host_int > 0){
      host_int <- 0
    }
    if(sym_int < host_int){
      host_res <- host_res - ((host_int - sym_int) * host_res)
    }
  }
  else{
    host_res <- host_donate * sym_int * synergy + host_res
  }
  return(host_res)
}


sym_res_f <- function(host_int, sym_int){
  host_res <- 0
  host_donate <- 0
  sym_res <- 0 
  if(host_int < 0){
    host_res <- res_d - (host_int * res_d * -1)
  } 
  else{
    host_donate <- host_int * res_d
    host_res <- res_d - (host_int * res_d)
  }
  
  if(sym_int < 0){
    if(host_int > 0){
      host_int <- 0
    }
    if(sym_int < host_int){
      sym_res <- (host_int - sym_int)*host_res
    }
    sym_res <- sym_res + host_donate
  }
  else{ # mutualistic sym
    sym_res <- host_donate - (host_donate * sym_int)
  }
  return(sym_res)
}


#### parameters #### 
res_d <- 100 # RES_DISTRIBUTE 
synergy <- 3
shtr <- 500 # SYM_HORIZ_TRANS_RES 
hrr <- 1000 # HOST_REPRO_RES 
updates <- 100000 # UPDATES -1
x <- rep(((0:20)/10)-1,21)
y <- sort(rep(((0:20)/10)-1,21))
vt_vals <- seq(0,1,0.5)

#### plotter #### 
intval_m <- data.frame(host_int = x, sym_int = y)

intval_m <- intval_m %>%
  mutate(sym_res = mapply(sym_res_f, host_int, sym_int),
         host_res = mapply(host_res_f, host_int, sym_int),
         host_repro_rate = host_res / hrr,
         sym_horiz_rate = sym_res / shtr)
a <- uncount(intval_m, length(vt_vals)) %>% group_by(host_int,sym_int) %>% 
  mutate(vt = vt_vals[row_number()],
         sym_vt_rate = host_repro_rate * vt,
         sym_repro_rate = sym_horiz_rate + sym_vt_rate)

host_plot <- intval_m %>%
  ggplot(aes(x = host_int, y = sym_int, fill = host_repro_rate)) + 
  geom_tile() + scale_fill_viridis_c(option="plasma", limits=c(0, 0.301)) + 
  facet_nested_theme + 
  labs(x = "Host interaction value", y = "Symbiont interaction value", 
       fill="Host\nreproductive\nrate")

host_plot

max(intval_m$sym_horiz_rate)
sym_plot <- ggplot(a, aes(x=host_int,y=sym_int,fill=sym_repro_rate)) + 
  geom_tile() +
  scale_fill_viridis_c(name="Symbiont\nreproductive\nrate",option="plasma", limits=c(0, 0.301)) +
  xlab("Host interaction value") + ylab("Symbiont interaction value") + 
  facet_wrap(~paste0("VT rate  ",vt)) + facet_nested_theme
sym_plot
