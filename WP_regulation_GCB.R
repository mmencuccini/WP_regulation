#___________________________________________________________________________


# CrUCIAL LINES


#___________________________________________________________________________
# lines 1234, 1247, and block 1273-1337 to be modified to go from individual tree to species average
# lines 1353 to be modified if multiple regression estimates of G7K are needed

# Notes Los Alamos 
# from Sanna Sevanto

{
# The WP data-set contains periodic tree level WP data from 2006 thru 2012. 
# Measurements were made at either predawn or midday. Data Qa/Qc has been performed on these files.
# PJ day refers to days since start of project (i.e., 1/1/2006).  
# The treatment classes provided in the file 
# are as follows; ambient (1), drought (2), cover-control (3), and irrigation (4).
# The experiment used plot aspect as the blocking factor.
# There are 3 different replicate blocks and block classifications designated
# in the files; flat aspect (1), north aspect (2), and south aspect (3).?  This will be obvious when viewing 
# the files.

# Tree numbers are always grouped by species as follows (regardless of plot); 
# Trees 1-5 are original Pinus edulis, 
# Trees 6-10 are original Juniper monosperma.
# When one of these original trees died, an additional tree in the plot was added to retain an adequate 
# sample size over time (i.e., multiple years+).
# These additional trees are grouped as follows; 
# Trees 11-15 are replacement Pinus edulis, Trees 16-20 are replacement Juniper monosperma.

# Replacement is used here in a more restricted sense, as these 
# additional trees have their separate and unique tree designation number.
# So, in differing plots you will have 
# differing numbers of trees depending on; 1) the number of trees for which data was collected, and 2) how many
# additional replacement trees had to be designated due to mortality (or partial mortality) of original trees.
# Many plots have n=10 trees, based on the original T1-T5 & T6-10 designation, as these particular plots did 
# not experience mortality.
# However, a plot like P10 has a total of n=16 trees. In P10, the original T1-5 T6-T10 trees are listed, 
# a replacement Pinon (T11) is listed, and five additional/replacement junipers (T16-T20).
# In some cases you will see data present at the same time for both original and replacement junipers 
# (plots 6 & 10).
# This is fine, as juniper experiences a slow/partial canopy dieback, so we monitored the original and replacement 
# trees at the same time in these two plots. Finally, we only provide data on trees for which data was collected 
# (so, e.g., in some instances you may only have n=4 cols of data for a particular species in a particular plot). 


# List of trees that died from the Sevilleta experiment. 
# These are based on that the team stopped measuring leaf water potential because there was no living branches in the tree.

# Plot 4: Trees 1 and 2 died some time between July and Oct 2009
# Plot 5: Tree 5 died in 2008
# Plot 6: Trees 1, 4 &5 died in 2008
# Plot 6: Trees 2, 3, 11 & 12 died some time between July and Oct 2009
# Plot 6: Tree 9 died in June 2011
# Plot 6: Tree 10 was reported dead in March 2011
# Plot 9: Tree 5 died in 2008
# Plot 10: Trees 1, 2, 4 & 5 died in 2008
# Plot 10: Trees 3 and 10 died some time between July and Oct 2009
# 
# 
# For SUMO experiment
# Tree 109 died in 2012
# 
# After that nothing died until 2016 when tree number 15 died.
  
}

#### Library####
library(conflicted)
library(devtools)
library(RGraphics)
library(gridExtra)
library(optimx)
library(tidyverse)
library(data.table) # loads data.table and set functions from library
library(broom)
library(broom.mixed)
library(janitor)
library(stringr)
library(stringi)
library(RColorBrewer)
library(splines)
library(lattice)
library(latticeExtra)
library(car)
library(scales) # controls formats on figures axes
library(rlang)
library(plotrix)
library(sjPlot)
library(visreg)
library(ggplot2)
library(ggrepel)
library(ggformula)
library(ggforce)
library(ggrepel)
library(plotly)
library(ggmap)
library(maps)
library(cowplot)
library(grid)
library(lmerTest) # this gives lsmeans; otherwise use emmeans package with more options
library(lme4)
library(nlme)
library(MuMIn)
library(multcomp)
library(emmeans) # use CLD for compact letter display (uses 'multcomp::cld')
library(Hmisc)
library(finalfit)
library(pbkrtest)
library(magrittr)
library(quantreg)
library(smatr)
library(splancs) # calculates area of convex hull

meanNA   <- function(x) mean(x, na.rm=TRUE)
sdNA     <- function(x) sqrt(var(x,na.rm=TRUE))
sterr    <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
medianNA <- function(x) median(x, na.rm=TRUE)
maxNA    <- function(x) max(x, na.rm=TRUE)
minNA    <- function(x) min(x, na.rm=TRUE)


# is.NaN method for data frames
is.nan.data.frame <- function(x)
{
  do.call(cbind, lapply(x, is.nan))
} 

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
lsos()

uniqueNA <- function(x) if(length((x))>1) {
  as.character(unique(x[!is.na(x)]))} else {unique(x)}

My_Theme = theme(
  axis.title.x = element_text(size=rel(1.5)),
  axis.text.x = element_text(size=rel(2)),
  axis.title.y = element_text(size=rel(1.5)),
  axis.text.y = element_text(size=rel(2)),
  title = element_text(size=rel(1.5)))


# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df, y, x){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

plt_prob  <- function(df, y, x, z){
  m <- lm(y ~ x * z, df);
  paste("italic(P)==", format(min(summary(m)$coefficients[5,4], 
                                  summary(m)$coefficients[6,4]),
                              digits = 2, trim = TRUE, nsmall=1) )
}

plt_r2 <- function(df, y, x, z){
  m <- lm(y ~ x * z, df);
  paste( "italic(r)^2==", 
         format(summary(m)$r.squared, digits = 2), sep = "")
}

floor_5yrs    = function(value){ return(value - value %% 5) }
ceiling_5yrs  = function(value){ return(floor_decade(value)+5) }
round_to_5yrs = function(value){ return(round(value / 5) * 5) }


# library conflicts  _________####
conflict_prefer('first', 'dplyr')
conflict_prefer('last', 'dplyr')
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("unique", "BiocGenerics")
conflict_prefer("lmer", "lmerTest")
conflict_prefer("cbind", "BiocGenerics")
conflict_prefer("rbind", "BiocGenerics")
conflict_prefer("as.data.frame", "BiocGenerics")
conflict_prefer("scores", "vegan")
conflict_prefer("pvalue", "coin")
conflict_prefer("filter", "dplyr")
conflict_prefer("which", "BiocGenerics")
conflict_prefer("qqnorm", "mgcViz")
conflict_prefer("paste", "base")
conflict_prefer("summarise", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("setdiff", "BiocGenerics")
conflict_prefer("year", "lubridate")
conflict_prefer("rename", "dplyr") 
conflict_prefer("mean", 'base')
conflict_prefer("sd", 'stats')
conflict_prefer("%+%", "ggplot2")
conflict_prefer("arrange", "dplyr")
conflict_prefer("geom_errorbarh", "ggplot2")
conflict_prefer("extract", "tidyr")
conflict_prefer("map", "purrr")
conflict_prefer("mask", "raster")
conflict_prefer('recode', 'dplyr')
conflict_prefer('transpose', 'purrr')
conflict_prefer("date", 'lubridate')
conflict_prefer("intersect", 'lubridate')
conflict_prefer("setdiff", 'lubridate')
conflict_prefer("union", 'lubridate')
conflict_prefer("loadings", "stats")





My_Theme = theme(
  axis.title.x = element_text(size=rel(1.5)),
  axis.text.x = element_text(size=rel(2)),
  axis.title.y = element_text(size=rel(1.5)),
  axis.text.y = element_text(size=rel(2)),
  title = element_text(size=rel(1.5)))


#### Los Alamos SUMO dataset ####

#  data from:

# https://data.ess-dive.lbl.gov/view/doi:10.15485/1454272
# Sevanto S ; Dickman T L ; Collins A ; Grossiord C ; Adams H ; Borrego I ; McDowell N (2018): 
# SUMO Micrometeorology. Vegetation Survival-Mortality (SUMO), ESS-DIVE repository. Dataset. 
# doi:10.15485/1454272 accessed via https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1454272 


# met data
{
  setwd("D:/Dropbox/My Documents/science/papers/my papers/Homoiohydry paper/data/Los Alamos SUMO")
  file <- 'Cerro Montoso Met.csv'
  met_LosAl <- read.csv(file, sep = ",", header = TRUE)
  
  
  # calculate saturated VP using formulas in Yuan et al 2019
  Pmsl <- 1013.25 # in hPa (pressure mean sea level)
  Z <- 1960 # altitude (m)
  
  # Saturation vapour pressure
  met_LosAl$Pmst <- with(met_LosAl,
                         Pmsl * ( (Temp_C + 273.16) / ( (Temp_C + 273.16) + 0.0065 * Z) ) ^ 5.625 )
  met_LosAl$fw <- with(met_LosAl,
                       (1 + 7 * 10^(-4) + 3.46 * 10 ^(-6) * Pmst ) )
  
  met_LosAl$SatVap <- with(met_LosAl,
                           6.112 * fw * exp( (17.67 * Temp_C) / (Temp_C + 243.5)) )
  
  met_LosAl$VPD <- with(met_LosAl, SatVap  - Vapor_Pressure ) / 10 # divide by 10 to go from hPa to in kPa
  
  met_LosAl <- met_LosAl %>%
    select(c(Date, Hour, VPD)) %>%
    filter(Hour == '4' | Hour == '13') %>%
    # mutate(chamber = rep.int(1:12, times = length(met_Drake$Hour)/12)) %>%
    # mutate(chamber = as.character(chamber)) %>%
    # mutate(chamber = paste('C0', chamber, sep ='')) %>%
    # group_by(date, chamber) %>%
    group_by(Date) %>%
    pivot_wider(., names_from = 'Hour', values_from = 'VPD') %>%
    mutate(VPDpd = `4`,
           VPDmd = `13`) %>%
    pivot_longer(., cols = c(`4`, `13`), names_to = 'Hour', values_to = 'VPD') %>%
    mutate(Time = case_when(Hour == '4' ~ 'Ypd',
                            Hour == '13' ~ 'Ymd')) 
  
  
  
  # first split DateTime into date and hours
  met_LosAl <-   met_LosAl %>% 
    mutate(
      Day = strsplit(Date, "/")  %>% 
        as.data.frame() %>% 
        t %>% 
        data.frame(stringsAsFactors = F) %>% 
        pull(1),
      Month = strsplit(Date, "/")  %>% 
        as.data.frame() %>% 
        t %>% 
        data.frame(stringsAsFactors = F) %>% 
        pull(2)
    ) %>%
    mutate(date = parse_date_time(Date, orders = 'dmy')) %>%
    mutate(Month = as.numeric(Month))

  
  
}


# https://data.ess-dive.lbl.gov/view/doi:10.15485/1439886
# Sevanto S ; Dickman T L ; Collins A ; Grossiord C ; Adams H ; Borrego I ; McDowell N (2018): 
# SUMO Leaf Water Potential. Vegetation Survival-Mortality (SUMO), ESS-DIVE repository. Dataset. 
# doi:10.15485/1439886

# WP data
{
  
  file <- 'Leaf_Water_Potential__LWP__2011-17.csv'
  WP_LosAl <- read.csv(file, sep = ",", header = TRUE)
  file <- 'Target tree information.csv'
  Info_LosAl <- read.csv(file, sep = ",", header = TRUE)
  
  
  WP_LosAl <- left_join( WP_LosAl, Info_LosAl, by = 'Tree_ID')
  
  
  # assumes pre-dawn at 4:00 and midday at 13:00
  WP_LosAl_long <- WP_LosAl %>%
    mutate(Ypdreal = Ypd) %>%
    mutate(Ymdreal = Ymd) %>%
    gather(., Time, WP, c("Ypd", "Ymd")) %>%
    mutate(Hour = case_when(Time == 'Ypd' ~ '4',
                            Time == 'Ymd' ~ '13')) %>%
    select(!c(Year, DOY))
  
  # WP_LosAl %>%
  #   +   group_by(Tree_ID) %>%
  #   +   mutate(DeltaY = (Ypd - Ymd)) %>%
  #   +   ggplot(aes(x = DeltaY)) + geom_vline(xintercept=0) +
  #   +   geom_histogram() + facet_wrap(vars(Tree_ID))
  # histograms show some days have Ymd>Ypd, perhaps indicating storms. I kept them.
}

# merge all data
{
  
  prob_slct <- 0.99
  

  Los_Al <- left_join(WP_LosAl_long, met_LosAl, by= c('Date', 'Hour', 'Time'))
  
  predawn <-   Los_Al %>%
    filter(Time == 'Ypd')
    
  # seasonal correction for differences in VPD across plots (see Micrometeorology_2012-2016.xlsx file)
  midday <-   Los_Al %>%
    filter(Time == 'Ymd') %>%
    mutate(VPD = ifelse(Treatment =='Heat' | Treatment =='Heat + Drought', 
                           ifelse(Month < 4 & Month > 9,
                                  VPD + 0.26, VPD + 0.49), VPD ) ) 
  
 Los_Al <- full_join(predawn, midday) 
 
 
 predawn <- Los_Al %>%
   filter(Time == 'Ypd') %>%
   mutate(VPDpd = VPD)
 
 midday <- Los_Al %>%
   filter(Time == 'Ymd') %>%
   mutate(VPDmd = VPD)
 
 Los_Al <- full_join(predawn, midday)
 

  

    # calculates ref GKpd using quantiles over dates by chamber
  pippo <- Los_Al %>%
    filter(Time == 'Ymd') %>%
    #___________________________________________________________________________
    #_CRUCIAL LINE
    # Group(Tree_ID) gets one set of values for all IDs
    # group_by(Tree_ID) %>%
    # an alternative is to calculate the 95% percentile by species 
    group_by(Species) %>%
    #___________________________________________________________________________
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD, + GKmaxpd * VPDpd
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    summarise(GKmaxmd = quantile(GKvary, probs = prob_slct, na.rm = TRUE)) %>%
    ungroup()  
  
  
  #___________________________________________________________________________
  #_CRUCIAL LINE
  # Los_Al <- left_join(Los_Al, pippo, by = 'Tree_ID')
  Los_Al <- left_join(Los_Al, pippo, by = 'Species')
  #___________________________________________________________________________

  
  pippo <- Los_Al %>%
    group_by(Tree_ID) %>%
    select(c(Tree_ID, Species, Treatment, Ymdreal, Ypdreal, VPDmd, VPDpd, GKmaxmd)) %>%
    summarise(VPDmdref = quantile(VPDmd, probs = 1-prob_slct, na.rm = TRUE),
              Ypdref = quantile(Ypdreal, probs = prob_slct, na.rm = TRUE)) %>%
    ungroup()
  
  Los_Al <- left_join(Los_Al, pippo, by = c('Tree_ID')) 
  

  # estimate GKmax at reference Ysoil and VPD
  #_CRUCIAL LINES ALSO INSIDE HERE, depending on group_by(Species) or (Tree_ID)
  {
    
    Ypdref = -0.1
    VPDref = +0.1
    
    pippo1 <-  Los_Al %>%
      filter(Time == 'Ymd',
             Species == 'Juniperus monosperma') %>%
      #___________________________________________________________________________
      group_by(Species) %>%
      # group_by(Tree_ID) %>%
      #___________________________________________________________________________
      mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD,  + GKmaxpd * VPDpd
      mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
      drop_na(GKvary) %>%
      #___________________________________________________________________________
      split(list(.$Species)) %>%
      # split(list(.$Tree_ID)) %>%
      #___________________________________________________________________________
      lapply(., function(x) { 
        predict(lm(GKvary ~ Ypdreal * VPDmd,
                   x[,]),
                newdata = data.frame(Ypdreal = Ypdref, VPDmd = VPDref)) } ) %>%
      do.call(rbind, .) %>%
      data.frame(.) %>%
      #___________________________________________________________________________
      mutate(Species = rownames(.)) %>%
      # mutate(Tree_ID = rownames(.)) %>%
      #___________________________________________________________________________
      ungroup() %>%
      rename(pred_GKmaxmd = X1) %>%
      left_join(Los_Al %>%
                  filter(Time == 'Ymd'), 
                ., 
                #___________________________________________________________________________
                by = c('Species')) %>%
      # by = c('Tree_ID')) %>%
    #___________________________________________________________________________
      ungroup()
    
    
    
    pippo2 <-  Los_Al %>%
      filter(Time == 'Ymd',
             Species == 'Pinus edulis') %>%
      #___________________________________________________________________________
      group_by(Species) %>%
      # group_by(Tree_ID) %>%
      #___________________________________________________________________________
      mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD,  + GKmaxpd * VPDpd
      mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
      drop_na(GKvary) %>%
      #___________________________________________________________________________
      split(list(.$Species)) %>%
      # split(list(.$Tree_ID)) %>%
      #___________________________________________________________________________
      lapply(., function(x) { predict(lm(GKvary ~ Ypdreal + VPDmd,
                                         x[,]), 
                                      newdata = data.frame(Ypdreal = Ypdref, VPDmd = VPDref)) } ) %>%
      do.call(rbind, .) %>%
      data.frame(.) %>%
      #___________________________________________________________________________
      mutate(Species = rownames(.)) %>%
      # mutate(Tree_ID = rownames(.)) %>%
      #___________________________________________________________________________
      ungroup() %>%
      rename(pred_GKmaxmd = X1) %>%
      left_join(Los_Al %>%
                  filter(Time == 'Ymd'), ., 
      #___________________________________________________________________________
      by = c('Species')) %>%
      # by = c('Tree_ID')) %>%
    #___________________________________________________________________________
      ungroup()
    
    pippo <- full_join(pippo1, pippo2) %>% 
      drop_na(pred_GKmaxmd)
    
    Los_Al <- full_join(Los_Al %>% 
                          filter(Time == 'Ypd'),
                        pippo )
    
  }
  
  #___________________________________________________________________________
  #_CRUCIAL LINE
  # decide whether to use quantile GKmax or predicted from multiple regression
  # Los_Al <- Los_Al %>%
  #   mutate(GKmaxmd = pred_GKmaxmd)
  #___________________________________________________________________________
  
  predawn <- Los_Al %>%
    filter(Time == 'Ypd') %>%
    group_by(Tree_ID) %>%
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD, + GKmaxpd * VPDpd
    mutate(Ytf = Ysoil) %>% # i.e., in practice, Ytf = Tpd at pd,  - GKmaxpd * VPDpd
    mutate(Reg_tot = WP - Ytf)%>% 
    mutate(sequence = 1) %>%
    ungroup() 
  
  midday <- Los_Al %>%
    filter(Time == 'Ymd') %>%
    group_by(Tree_ID) %>%
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD, + GKmaxpd * VPDpd
    mutate(Ytf = Ysoil - GKmaxmd * VPD) %>% 
    mutate(Reg_tot = WP - Ytf) %>%
    mutate(sequence = 2) %>%
    ungroup()  
  
  # here forcing does not account for VPD-Ysoil interaction
  Los_Al <- full_join(predawn, midday) %>%
    arrange(sequence, .by_group=TRUE) %>%
    mutate(Species = Tree_ID) %>%
    mutate(Species = substring(Species, first = 1, last = 1)) %>%
    mutate(Species = as.factor(Species)) %>%
    mutate(Spp = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Spp = fct_recode(Spp, 'Juniperus monosperma' = 'J')) %>%
    mutate(Forc_Ysoil = Ysoil) %>%
    mutate(Forc_VPD = ((-1) * GKmaxmd * VPD)) %>%
    mutate(Reg_Ysoil = Ymdreal - Ysoil) %>%
    mutate(Reg_VPD = Ymdreal - ((-1) * GKmaxmd * VPD)) %>%
    ungroup()
    
  
  pippo <- Los_Al %>%
    filter(Time=='Ymd') %>%
    split(list(.$Species, .$Treatment)) %>%
    lapply(., function(x) { lm(VPDmd ~ Ypdreal,
                                x[,])$coefficients } ) %>%
    do.call(rbind, .) %>%
    data.frame(.) %>%
    mutate(SppTreat = rownames(.)) %>%
    mutate(Species = substring(SppTreat, first = 1, last = 1),
           Treatment = substring(SppTreat, first = 3, last = 100)) %>%
    select(!c(X.Intercept., SppTreat)) %>%
    relocate(c(Species, Treatment)) %>%
    rename(dVPDdYpd = Ypdreal) %>%
    ungroup()
  
  # reset slopes to zero for all pines since P always >> 0.05
  pippo <- full_join(pippo %>% 
                       filter(Species!='P'),
                     pippo %>% 
                       filter(Species=='P') %>%
                       mutate(dVPDdYpd = 0))
    
    
    
    Los_Al <- left_join(Los_Al, pippo, by = c('Species', 'Treatment')) %>%
      # (Eqn.5a)
      mutate(Forc_VPD_net = Forc_VPD + GKmaxmd * (Ysoil-Ypdref) * dVPDdYpd) %>%
      # (Eqn.7a)
      mutate(Reg_VPD_net = Ymdreal - Forc_VPD_net)
    
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      split(list(.$Species, .$Treatment)) %>%
      lapply(., function(x) { lm(Ypdreal ~ VPDmd,
                                 x[,])$coefficients } ) %>%
      do.call(rbind, .) %>%
      data.frame(.) %>%
      mutate(SppTreat = rownames(.)) %>%
      mutate(Species = substring(SppTreat, first = 1, last = 1),
             Treatment = substring(SppTreat, first = 3, last = 100)) %>%
      select(!c(X.Intercept., SppTreat)) %>%
      relocate(c(Species, Treatment)) %>%
      rename(dYpddVPD = VPDmd) %>%
      ungroup()

    # reset slopes to zero for all pines since P always >> 0.05
    pippo <- full_join(pippo %>% 
                         filter(Species!='P'),
                       pippo %>% 
                         filter(Species=='P') %>%
                         mutate(dYpddVPD = 0))
    
    
    
    Los_Al <- left_join(Los_Al, pippo, by = c('Species', 'Treatment')) %>%
      # (Eqn.5b)
      mutate(Forc_Ysoil_net = Ysoil - (VPDmd-VPDmdref) * dYpddVPD) %>%
      # (Eqn.5c)
      mutate(Forc_cov_net = (VPDmd-VPDmdref) * dYpddVPD - GKmaxmd * (Ysoil-Ypdref) * dVPDdYpd ) %>%
      # (Eqn.7b)
      mutate(Reg_Ysoil_net = Ymdreal - Forc_Ysoil_net)  %>%
      # (Eqn.7c)
      mutate(Reg_cov_Ysoil_VPD = Ymdreal - Forc_cov_net)
    
    
  }

# total forcing
{

  #___________________________________________________________________
  # FIGURE 3 TOP
  # plot of observations against total forcing

  
  
  ann_text1 <- data.frame(Ytf = c(-5), 
                          WP  = c(-7),
                          lab = "Text",
                          Spp = factor('Pinus edulis',
                                       levels = c("Pinus edulis","Juniperus monosperma")),
                          Treatment = factor('Heat + Drought',
                                             levels = c("Ambient","Chamber Control","Drought",
                                                        'Heat', 'Heat + Drought')) )
  
  ann_text2 <- data.frame(lab = "Text",
                          Spp = factor('Pinus edulis',
                                       levels = c("Pinus edulis","Juniperus monosperma")),
                          Treatment = factor('Heat + Drought',
                                             levels = c("Ambient","Chamber Control","Drought",
                                                        'Heat', 'Heat + Drought')) )
  
  
    
  Los_Al  %>%
    mutate(Spp = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Spp = fct_recode(Spp, 'Juniperus monosperma' = 'J')) %>%
    filter(Time == 'Ymd') %>%
    ggplot(aes(x = Ytf,
               y = WP, group=Treatment)) +
    geom_point(size=2, aes(color = Treatment, alpha = 0.5)) +
    geom_point(size=2, color = 'black',  
               data = function(x) {x[x$Time %in% c("Ypd"), ]} ) +
    geom_quantile(formula =  y ~ x,
                  method = 'rq',
                  quantiles = c(0.95), col = 'black', size = 1.5) +
    geom_abline(intercept=0, slope = 1, linetype = 'dotted', size=1.5) +
    # geom_line(aes(x=Ytf,
    #               y=(-1)*exp(InvRegrWP)),
    #           linetype = 'dashed', size = 1.5, color = 'black') +
    geom_smooth(method='lm',
                formula= y~x,
                se = FALSE,
                aes(x = Ytf,
                    y = WP,
                    group = Treatment,
                    color = Treatment),
                size = 2) +
    theme_classic() +
    My_Theme +
    theme(legend.position="none") +
    labs(title = "Water status regulation Juniper monosperma / Pinus edulis") +
    labs(x = expression(Psi[hf] ~ ', MPa'),
         y = expression(Psi[md] ~ ', MPa')) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    xlim(-10, 0) +
    ylim(-7.5,0) +
    scale_x_continuous(breaks = c(0,-2,-4,-6,-8, -10), limits = c(-10,0)) +
    facet_grid(Spp ~ Treatment) +
    geom_text(data = ann_text1,
              label = expression(A[s]),
              size = 6,
              fontface = "italic",
              color = 'black') +
    geom_segment(data = ann_text2,
                 x    = c(-4),
                 xend = c(-0.5),
                 y    = c(-7),
                 yend = c(-7),
                 color = c("black"), 
                 linewidth = 1.5, linetype = c(1)) +
    panel_border(colour='grey85', size=1, linetype=1)
  
}

# partitioning of forcing 
{
  
  #______________________________________________________________________
  # FIGURE 5 TOp partitioning of NET forcing
  # forcing by VPD and Ysoil accounting for covariation (simple formula)
  
  
  ann_text1 <- data.frame(Ytf = -4, 
                          Forc_VPD_net = c(-7, -9),
                          lab = "Text",
                          Spp = factor('Pinus edulis',
                                       levels = c("Pinus edulis","Juniperus monosperma")),
                          Treatment = factor('Heat + Drought',
                                             levels = c("Ambient","Chamber Control","Drought",
                                                        'Heat', 'Heat + Drought'))
  )
  
  ann_text2 <- data.frame(Ytf = -4, 
                          Forc_VPD_net = c(-8),
                          lab = "Text",
                          Spp = factor('Pinus edulis',
                                       levels = c("Pinus edulis","Juniperus monosperma")),
                          Treatment = factor('Heat + Drought',
                                             levels = c("Ambient","Chamber Control","Drought",
                                                        'Heat', 'Heat + Drought'))
  )
  
  Los_Al %>%
    filter(Time=='Ymd') %>%
    mutate(Spp = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Spp = fct_recode(Spp, 'Juniperus monosperma' = 'J')) %>%
    ggplot(aes(x = Ytf, 
               y = Forc_VPD_net)) +
    geom_point(aes(x = Ytf, 
                   y = Forc_VPD_net),
               size = 3,
               shape = 16,
               color = 'red',
               alpha = 0.4) +
    geom_smooth(method='lm',
                formula= y~x,
                se = TRUE,
                color = 'red',
                size = 1.2) +
    geom_point(aes(x = Ytf,
                   y = Forc_Ysoil_net),
               size = 2,
               shape = 15,
               color = 'black',
               alpha = 0.4) +
    geom_smooth(aes(x = Ytf,
                    y = Forc_Ysoil_net),
                method='lm',
                formula= y~x,
                se = TRUE,
                color = 'black',
                size = 1.2) +
    geom_point(aes(x = Ytf,
                   y = Forc_cov_net),
               size = 2,
               shape = 15,
               color = 'blue',
               alpha = 0.4) +
    geom_smooth(aes(x = Ytf,
                    y = Forc_cov_net),
                method='lm',
                formula= y~x,
                se = TRUE,
                color = 'blue',
                size = 1.2) +
    geom_line(aes(x = Ytf,
                  y = Ytf),
              linetype = 'dashed', size = 1.5) +
    scale_color_viridis_c(option = "plasma") +
    theme_classic() +
    My_Theme +
    xlim(-10,0) +
    ylim(-10,0) +
    scale_x_continuous(breaks = c(0,-2,-4,-6,-8, -10), limits = c(-10,0)) +
    theme(plot.margin = rep(unit(0.5, "cm"), times = 4) ) +
    labs(title = "SUMO Experiment, Sevilleta, USA") +
    xlab('Total hydraulic forcing, MPa') +
    labs(y = expression('Net forcing by VPD and' ~ Psi[soil] ~ ', MPa') ) +
    facet_grid(Spp ~ Treatment) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    geom_text(data = ann_text1,label = c("Covar", 'VPD'),
              size = 6,
              fontface = "bold",
              color = c('blue', 'red')) +
    geom_text(data = ann_text2,label = c(expression(bold(Psi[soil]))),
              size = 6,
              fontface = "bold",
              color = c('black')) +
    panel_border(colour='grey85', size=1, linetype=1)
  
 }

# Stringency of regulation
{
  
  #____________________________________________________________________________
  # FIGURE 6
  # NET regulation against VPD
  Los_Al %>%
    mutate(Spp = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Spp = fct_recode(Spp, 'Juniperus monosperma' = 'J')) %>%
    filter(Time == 'Ymd') %>%
    ggplot(aes(x = Ytf, 
               y = Reg_VPD_net,
               group=Treatment)) +
    geom_point(aes(x = Ytf,
                   y = Reg_tot),
               size = 2,
               shape = 15,
               color = 'black',
               alpha = 1/2) +
    geom_point(size=2, 
               shape = 16,
               alpha = 1,
               aes(color = Treatment)
    ) +
    geom_hline(yintercept=0, size=1.5) +
    geom_abline(intercept=0, slope = 1, linetype = 'dotted', size=1.5) +
    geom_smooth(method='lm',
                formula= y ~ x,
                se = FALSE,
                aes(color = Treatment),
                size = 2) +
    geom_smooth(aes(x = Ytf,
                    y = Reg_tot),
                method='lm',
                formula= y ~ x,
                se = FALSE,
                size = 1,
                color = 'black') +
    theme_classic() +
    My_Theme +
    theme(plot.margin = rep(unit(0.5, "cm"), times = 4) ) +
    theme(legend.position="none") +
    scale_x_continuous(breaks = c(0,-2,-4,-6,-8, -10), limits = c(-10,0)) +
    labs(title = "Regulation against VPD, SUMO") +
    labs(x = expression(Psi[hf] ~ ', MPa'),
         y = expression(italic(S)[ws] ~ ',' ~ italic(S)[VPD] ~ ', MPa')) +
    facet_grid(Spp ~ Treatment) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") )   +
    panel_border(colour='grey85', size=1, linetype=1)
  


}

# Isohydricity plots (following Martinez-Vilalta et al)
{
  
  
  ann_text1 <- data.frame(Ypdreal = c(-5, -5), 
                          Ymdreal = c(-8, -9),
                          lab = "Text",
                          Spp = factor('Pinus edulis',
                                       levels = c("Pinus edulis","Juniperus monosperma")),
                          Treatment = factor('Heat + Drought',
                                             levels = c("Ambient","Chamber Control","Drought",
                                                        'Heat', 'Heat + Drought')) )

  ann_text2 <- data.frame(lab = rep("Text",2),
                          Spp = factor('Pinus edulis',
                                       levels = c("Pinus edulis","Juniperus monosperma")),
                          Treatment = factor('Heat + Drought',
                                             levels = c("Ambient","Chamber Control","Drought",
                                                        'Heat', 'Heat + Drought')) )
  
  Los_Al %>%
    filter(Time=='Ymd') %>%
    mutate(Spp = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Spp = fct_recode(Spp, 'Juniperus monosperma' = 'J')) %>%
    ggplot(aes(x = Ypdreal,
               y = Ymdreal,
               color = Treatment)) +
    geom_point(size=3) +
    geom_smooth(method = 'lm',
                 formula =  y ~ x,
                 se = FALSE,
                linetype = 'dashed', size = 1.5,
                 color = 'black') +
    geom_quantile(formula =  y ~ x,
                  method = 'rq',
                  quantiles = c(0.05), col = 'black', size = 1.5) +
    geom_abline(intercept=0, slope = 1, linetype = 'dotted', size=1.5) +
    # ylim(-2.7, 0) +
    scale_x_continuous(breaks = c(0,-2,-4,-6,-8, -10), limits = c(-10,0)) +
    theme_classic() +
    My_Theme +
    theme(legend.position="none") +
    labs(title = "Water status regulation") +
    xlab(expression(Psi[pd])) +
    ylab(expression(Psi[md])) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    facet_grid(Spp ~ Treatment) + 
    geom_text(data = ann_text1,
              label = c(expression(sigma), 'HS'),
              size = 6,
              fontface = "bold",
              color = 'black') +
    geom_segment(data = ann_text2,
             x    = c(-4, -4),
             xend = c(-0.5, -0.5),
             y    = c(-8, -9),
             yend = c(-8, -9),
             color = c("black","black"), 
             linewidth = 1.5, linetype = c(2,1)) +
    panel_border(colour='grey85', size=1, linetype=1)

  }

# STATS of regulation
{
  
# response of quantile GKmaxmd to treatments
  {

    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      group_by(Tree_ID, Treatment, Species) %>%
      select(Tree_ID, Treatment, Species, GKmaxmd) %>%
      summarise(GKmaxmd = mean(GKmaxmd)) %>%
      ungroup()
    
    
    summary(lm(GKmaxmd ~ Species * Treatment, data = pippo))
    anova(lm(GKmaxmd ~ Species * Treatment, data = pippo))

    m0 <-  lm(GKmaxmd ~ Species * Treatment, data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
    cld(emmeans(m0, ~ Species * Treatment, var = 'Ytf') )
    
    
    
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      filter(Species == 'P') %>% 
      group_by(Tree_ID, Treatment, Species) %>%
      select(Tree_ID, Treatment, Species, GKmaxmd ) %>%
      summarise(GKmaxmd = mean(GKmaxmd )) %>%
      ungroup()
    
    
    # test only valid if GKmax calculated by individual, not by species
    summary(lm(GKmaxmd ~ Treatment, data = pippo))
    anova(lm(GKmaxmd ~ Treatment, data = pippo))
    
    m0 <-  lm(GKmaxmd ~ Treatment, data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      filter(Species == 'J') %>% 
      group_by(Tree_ID, Treatment, Species) %>%
      select(Tree_ID, Treatment, Species, GKmaxmd ) %>%
      summarise(GKmaxmd = mean(GKmaxmd )) %>%
      ungroup()
    
    
    # test only valid if GKmax calculated by individual, not by species
    summary(lm(GKmaxmd ~ Treatment, data = pippo))
    anova(lm(GKmaxmd ~ Treatment, data = pippo))
    
    m0 <-  lm(GKmaxmd ~ Treatment, data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
  }
  
# response of predicted GKmaxmd (at fixed VPD and Ys) to treatments
  {

    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      group_by(Tree_ID, Treatment, Species) %>%
      select(Tree_ID, Treatment, Species, pred_GKmaxmd ) %>%
      summarise(pred_GKmaxmd = mean(pred_GKmaxmd )) %>%
      ungroup()
    
    # test only valid if GKmax calculated by individual, not by species
    summary(lm(pred_GKmaxmd ~ Species * Treatment, data = pippo))
    anova(lm(pred_GKmaxmd ~ Species * Treatment, data = pippo))
    
    m0 <-  lm(pred_GKmaxmd ~ Species * Treatment, data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
    cld(emmeans(m0, ~ Species * Treatment, var = 'Ytf') )
    
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      filter(Species == 'P') %>% 
      group_by(Tree_ID, Treatment, Species) %>%
      select(Tree_ID, Treatment, Species, pred_GKmaxmd ) %>%
      summarise(pred_GKmaxmd = mean(pred_GKmaxmd )) %>%
      ungroup()

    
    # test only valid if GKmax calculated by individual, not by species
    summary(lm(pred_GKmaxmd ~ Treatment, data = pippo))
    anova(lm(pred_GKmaxmd ~ Treatment, data = pippo))
    
    m0 <-  lm(pred_GKmaxmd ~ Treatment, data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )

    
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      filter(Species == 'J') %>% 
      group_by(Tree_ID, Treatment, Species) %>%
      select(Tree_ID, Treatment, Species, pred_GKmaxmd ) %>%
      summarise(pred_GKmaxmd = mean(pred_GKmaxmd )) %>%
      ungroup()
    
    
    # test only valid if GKmax calculated by individual, not by species
    summary(lm(pred_GKmaxmd ~ Treatment, data = pippo))
    anova(lm(pred_GKmaxmd ~ Treatment, data = pippo))
    
    m0 <-  lm(pred_GKmaxmd ~ Treatment, data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )

    
  }
  
# sensitivity of Ymd to Ypd alone
  {
    
    # separate by species
       pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd') %>%
      drop_na(Ypdreal)
    
    
    summary(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
    
    summary(lmer(Ymdreal ~ Ypdreal * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ypdreal * Treatment + (1|Tree_ID), data = pippo))
    
    
    pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd') %>%
      drop_na(Ypdreal)
    
    
    summary(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
    
    summary(lmer(Ymdreal ~ Ypdreal * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ypdreal * Treatment + (1|Tree_ID), data = pippo))
    
    # two species together
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      drop_na(Ypdreal, Ymdreal)
    
    summary(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
    
    summary(lmer(Ymdreal ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    check_model(m0)
    

    summary(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    
      }
  
  # sensitivity of DeltaY to Ypd alone
  {
    
    # separate by species
    pippo <- Los_Al %>%
      mutate(DeltaY = Ypdreal - Ymdreal) %>% 
      filter(Species == 'P') %>% 
      filter(Time=='Ymd') %>%
      drop_na(Ypdreal, DeltaY)
    
    
    summary(lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
    
    pippo <- Los_Al %>%
      mutate(DeltaY = Ypdreal - Ymdreal) %>% 
      filter(Species == 'J') %>% 
      filter(Time=='Ymd') %>%
      drop_na(Ypdreal, DeltaY)
    
    
    summary(lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(DeltaY ~ Ypdreal + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
    
    # two species together
    
    pippo <- Los_Al %>%
      mutate(DeltaY = Ypdreal - Ymdreal) %>% 
      filter(Time=='Ymd') %>%
      drop_na(Ypdreal, Ymdreal)
    
    summary(lmer(DeltaY ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(DeltaY ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(DeltaY ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(DeltaY ~ Species * Ypdreal + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Species', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
    CLD(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
    multcomp::cld(glht(m0, mcp(Species = "Tukey")))
    emmip(m0, Species ~ Treatment, cov.reduce = range)
    check_model(m0)
    
    
    
  }
  
  # sensitivity of Ymd to Ytf or Forc_VPD
  {
    
    pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd') %>%
      mutate(logYmd = log(-Ymdreal)) %>%
      mutate(logYtf = log(-Ytf)) %>%
      filter_at(vars(logYmd), all_vars(!is.infinite(.))) %>%
      drop_na(c(logYmd, logYtf))
    
    # it is concave up
    summary(lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    
    anova(lmer(Ymdreal ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    check_model(m0)
    
    summary(lmer(logYmd ~ logYtf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(logYmd ~ logYtf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(logYmd ~ logYtf + Treatment + (1|Tree_ID), data = pippo))
    
    
    pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd') %>%
      mutate(logYmd = log(-Ymdreal)) %>%
      mutate(logYtf = log(-Ytf)) %>%
      filter_at(vars(logYmd,logYtf), all_vars(!is.infinite(.))) %>%
      drop_na(logYmd, logYtf)
    
    # it is also concave up
    summary(lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    check_model(m0)
    
    m0 <-  lmer(Ymdreal ~ Ytf + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    check_model(m0)
    

    
    # two species together and plot
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      drop_na(Ytf, Ymdreal)
    
    summary(lmer(Ymdreal ~ Ytf * Species + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Ytf * Species + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ytf * Species + Treatment + (1|Tree_ID), data = pippo))
    
    anova(lmer(Ymdreal ~ Ytf * Species + Species * Treatment + (1|Tree_ID), data = pippo))

        m0 <-  lmer(Ymdreal ~ Ytf * Species + Treatment + (1|Tree_ID), data = pippo)
    check_model(m0)
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      drop_na(Forc_VPD, Ymdreal)
    
    
    summary(lmer(Ymdreal ~ Species * poly(Forc_VPD, degree = 1) * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Ymdreal ~ Species * poly(Forc_VPD, degree = 1) * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Species * poly(Forc_VPD, degree = 1) * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Species * poly(Forc_VPD, degree = 2) * Treatment + (1|Tree_ID), data = pippo)
    check_model(m0)
    
  
  }
  
  # partitioning of forcing
  {
    
    
    # separate by species
    pippo <- Los_Al %>%
      filter(Species == 'P') %>%
      filter(Time=='Ymd') %>%
      mutate(pForc_VPD = Forc_VPD_net / Ytf) %>%
      drop_na(Ypdreal, Ytf)
    
# proportion of NET forcing by VPD    
    summary(lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo))

    
    m0 <-  lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)

        
    # slope of Forc_VPD by treatment
    summary(lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    # differences in slope with Ytf
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    # absolute differences in VPD forcing by treatment adjusted by total forcing
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    
    # absolute differences in VPD forcing by treatment
    # boundary (singular) fit: see help('isSingular') (variance TreeID=0)
    m0 <-  lmer(Forc_VPD_net ~ Treatment + (1|Tree_ID), data = pippo)
    anova(lmer(Forc_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    m0 <-  lm(Forc_VPD_net ~ Treatment, data = pippo)
    anova(lm(Forc_VPD_net ~ Treatment, data = pippo))
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Forc_VPD_net') )

    plot(m0)

    
    
    # slope of Forc_Ysoil by treatment
    summary(lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    
    
    
    pippo <- Los_Al %>%
      filter(Species == 'J') %>%
      filter(Time=='Ymd') %>%
      mutate(pForc_VPD = Forc_VPD_net / Ytf) %>%
      drop_na(Ypdreal, Ytf)
    
    
    # proportion of forcing by VPD    
    summary(lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(pForc_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    # slope of NET Forc_VPD by treatment
    summary(lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Forc_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    
    # slope of Forc_Ysoil by treatment
    summary(lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Forc_Ysoil_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    
    
    
    
    
    
    # both species together
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      mutate(pForc_VPD = Forc_VPD_net / Ytf) %>%
      drop_na(Ypdreal, Ytf)
    
    
    summary(lmer(pForc_VPD ~ Ytf * Treatment * Spp + (1|Tree_ID), data = pippo))
    anova(lmer(pForc_VPD ~ Ytf * Treatment * Spp + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(pForc_VPD ~ Ytf * Treatment * Spp + (1|Tree_ID), data = pippo))
    
  
    }
  
  # sensitivity of Reg_tot
  {
    
    pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd') %>%
      drop_na(Reg_tot)
    
    summary(lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    
    summary(lmer(Reg_tot ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_tot ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_tot ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    check_model(m0)
    
    # whether slope varies by treatment
    m0 <-  lmer(Reg_tot ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    
    # whether absolute stringency varies by treatment
    m0 <-  lmer(Reg_tot ~ Treatment + (1|Tree_ID), data = pippo)
    anova(lmer(Reg_tot ~ Treatment + (1|Tree_ID), data = pippo))
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
    pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd') %>%
      mutate(logReg_tot = log(Reg_tot)) %>%
      filter_at(vars(logReg_tot), all_vars(!is.infinite(.))) %>%
      drop_na(logReg_tot)

    
    summary(lmer(log(Reg_tot) ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(log(Reg_tot) ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(log(Reg_tot) ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    
        
    pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd') %>%
      drop_na(Reg_tot)
    
    summary(lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo))
    
    anova(lmer(Reg_tot ~ Ytf * Treatment + (1|Tree_ID), data = pippo))

        # whether slope varies by treatment
    m0 <-  lmer(Reg_tot ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    
    # whether absolute stringency varies by treatment
    m0 <-  lmer(Reg_tot ~ Treatment + (1|Tree_ID), data = pippo)
    anova(lmer(Reg_tot ~ Treatment + (1|Tree_ID), data = pippo))
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
        
    pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd') %>%
      mutate(logReg_tot = log(Reg_tot)) %>%
      filter_at(vars(logReg_tot), all_vars(!is.infinite(.))) %>%
      drop_na(logReg_tot)

    
    summary(lmer(log(Reg_tot) ~ Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(log(Reg_tot) ~ Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(log(Reg_tot) ~ Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_tot ~ Ytf + Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    check_model(m0)
    
    
    # two species together
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      drop_na(Reg_tot)
    
    summary(lmer(Reg_tot ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_tot ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_tot ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_tot ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    cld(emmeans(m0, pairwise ~ Species, var = 'Ytf') )
    cld(emmeans(m0, ~ Species * Treatment, var = 'Ytf') )
    cld(emmeans(m0, ~ Species * Ytf * Treatment, var = 'Ytf') )
    
    
    
    pippo <- Los_Al %>%
      filter(Time=='Ymd') %>%
      mutate(logReg_tot = log(Reg_tot)) %>%
      filter_at(vars(logReg_tot), all_vars(!is.infinite(.))) %>%
      drop_na(logReg_tot)
    
    summary(lmer(log(Reg_tot) ~ Species * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(log(Reg_tot) ~ Species * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(log(Reg_tot) ~ Species * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_tot ~ Ytf * Species + Treatment + (1|Tree_ID), data = pippo)
    check_model(m0)
    
    
    
  }
  
  # regulation as proportion of Reg_tot
  {
    
    pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd') 
    
    summary(lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Reg_tot', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Reg_tot') )
    emmip(m0, Treatment ~ Reg_tot, cov.reduce = range)
    
    pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd') 
    
    summary(lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD ~ Reg_tot * Treatment + (1|Tree_ID), data = pippo))
    
    
  }
  
  # regulation against VPD
  {
    
    pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd')
    
    summary(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )

    # whether slope varies by treatment
    m0 <-  lmer(Reg_VPD ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    
    # whether absolute stringency varies by treatment
    m0 <-  lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo)
    anova(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
        pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd')
    
    summary(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD ~ Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Treatment') )

    # both species together
    pippo <- Los_Al %>%
      filter(Time=='Ymd')
    
    summary(lmer(Reg_VPD ~ Species + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD ~ Species + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD ~ Species + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD ~ Species + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Reg_VPD', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
    CLD(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
    emmip(m0, Species ~ Treatment, cov.reduce = range)
    check_model(m0)
    emmeans(m0, pairwise ~ Treatment | Species)
    
  }
  
  # NET regulation against VPD
  {
    
    pippo <- Los_Al %>%
      filter(Species == 'P') %>% 
      filter(Time=='Ymd')
    
    # whether slope varies by treatment
    anova(lmer(Reg_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    m0 <-  lmer(Reg_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
    # whether absolute stringency varies by treatment
    summary(lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )


        
    
    pippo <- Los_Al %>%
      filter(Species == 'J') %>% 
      filter(Time=='Ymd')
    

    # whether slope varies by treatment
    anova(lmer(Reg_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo))
    m0 <-  lmer(Reg_VPD_net ~ Ytf * Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
    # whether absolute stringency varies by treatment
    summary(lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD_net ~ Treatment + (1|Tree_ID), data = pippo)
    cld(emmeans(m0, pairwise ~ Treatment) )
    
    
    
    # both species together
    pippo <- Los_Al %>%
      filter(Time=='Ymd')
    
    summary(lmer(Reg_VPD_net ~ Species + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD_net ~ Species + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD_net ~ Species + Treatment + (1|Tree_ID), data = pippo))
    
    summary(lmer(Reg_VPD_net ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo))
    anova(lmer(Reg_VPD_net ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo))
    r.squaredGLMM(lmer(Reg_VPD_net ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo))
    
    m0 <-  lmer(Reg_VPD_net ~ Species * Ytf + Treatment + (1|Tree_ID), data = pippo)
    emtrends(m0, var = 'Reg_VPD_net', pairwise ~ Treatment)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    cld(emmeans(m0, Species ~ Ytf, var = 'Treatment') )
    cld(emmeans(m0, ~ Species * Treatment, var = 'Ytf') )
    check_model(m0)
    emmeans(m0, pairwise ~ Treatment | Species)
    
  }
  
}

# Hydroscape and Stringency areas
{
  
  # Total stringency area
  {
    
    # 1) convex hull from Reg-Yhf plot for all Sp/Treatment combinations (concept not correct)
    {
      
      
      hull <- data.frame(matrix(data = NA, nrow = 62, ncol = 5))
      colnames(hull) <- c('Species', 'Treatment', 'Tree_ID', 'A_above', 'A_below')
      
      m=1
      for(i in unique(Los_Al$Species)) {
        
        ifelse(m==1, n <- 1, n <- 33)
        
        for(j in unique(Los_Al$Treatment)) {
          
          pippo <- Los_Al %>%
            filter(Species == i) %>%
            filter(Treatment == j)
          
          for(l in unique(pippo$Tree_ID)) {
            
            pippo <- Los_Al %>%
              filter(Species == i) %>%
              filter(Treatment == j) %>%
              filter(Tree_ID == l) %>%
              drop_na(Ytf, Reg_tot)
            
            hull[n,1] <- i
            hull[n,2] <- j
            hull[n,3] <- l
            
            if(length(pippo$Species) > 0) {
              
              pippo1 <- pippo %>%
                select(Ytf, Reg_tot) %>% 
                drop_na(Ytf, Reg_tot) %>%
                filter(Reg_tot > 0) %>%
                full_join(.,
                          with(., data.frame(Ytf = c(minNA(Ytf), maxNA(Ytf)),
                                             Reg_tot = c(rep(0, 2))))
                )
              
              if(unique(is.finite(pippo1$Ytf))) {
                hull[n,4] <- pippo1 %>%
                  slice(chull(Ytf, Reg_tot)) %>%
                  as.matrix(.) %>%
                  areapl
              }
              
              pippo1 <- pippo %>%
                select(Ytf, Reg_tot) %>% 
                drop_na(Ytf, Reg_tot) %>%
                filter(Reg_tot < 0) %>%
                full_join(.,
                          with(., data.frame(Ytf = c(minNA(Ytf), maxNA(Ytf)),
                                             Reg_tot = c(rep(0, 2))))
                )
              
              if(unique(is.finite(pippo1$Ytf))) {
                
                hull[n,5] <- pippo1 %>%
                  slice(chull(Ytf, Reg_tot)) %>%
                  as.matrix(.) %>%
                  areapl
              }
            }
            n=n+1
          }
        }
        m=m+1
      }
    }
    
    # 2) convex hull from Ymd-Yhf plot for all Sp/Treatment combinations (concept not correct)
    {
      
      
      hull <- data.frame(matrix(data = NA, nrow = 62, ncol = 4))
      colnames(hull) <- c('Species', 'Treatment', 'Tree_ID', 'A_net')
      
      m=1
      for(i in unique(Los_Al$Species) ) {
        
        ifelse(m==1, n <- 1, n <- 33)
        
        for(j in unique(Los_Al$Treatment) ) {
          
          pippo <- Los_Al %>%
            filter(Species == i) %>%
            filter(Treatment == j)
          
          for(l in unique(pippo$Tree_ID) ) {
            
            pippo <- Los_Al %>%
              filter(Species == i) %>%
              filter(Treatment == j) %>%
              filter(Tree_ID == l) %>%
              drop_na(Ytf, Ymdreal)
            
            hull[n,1] <- i
            hull[n,2] <- j
            hull[n,3] <- l
            
            if(length(pippo$Species) > 0 ) {
              
              pippo1 <- pippo %>%
                select(Ytf, Ymdreal)
              
              if( unique(is.finite(pippo1$Ytf)) ) {
                hull[n,4] <- pippo1 %>%
                  slice(chull(x = Ytf, y = Ymdreal)) %>%
                  as.matrix(.) %>%
                  areapl
              }
            }
            n=n+1
          }
        }
        m=m+1
      }
    }
    
    # 3) figures and stats for hull and DeltaY
    {
      hull %>%
        drop_na() %>%
        mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
        mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
        group_by(Species, Treatment) %>%
        mutate(A_net = A_above - A_below) %>%
        summarise(Anet_mean = meanNA(A_net),
                  Anet_se   = sterr(A_net)) %>%
        ungroup() %>%
        ggplot(aes(y = Anet_mean,
                   x = Treatment),
               group = Species) + 
        geom_point(aes(color = Species),
                   position = position_dodge(width = 0.30), size = 3)  + 
        geom_errorbar(aes(ymin = Anet_mean - 1.96 * Anet_se,
                          ymax = Anet_mean + 1.96 * Anet_se,
                          group = Species,
                          color = Species),
                      position = position_dodge(width = 0.3),
                      width = 0.2, size = 1.2) +
        theme_classic() +
        My_Theme +
        xlab('Treatments') +
        ylab(expression('Total regulation stringency area,' ~MPa^2 ))
      
      
      Los_Al %>%
        filter(Time == 'Ymd') %>%
        ggplot(aes(x = Ytf,
                   y = Ymdreal-Ytf)) +
        geom_point() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_quantile(formula =  y ~ x,
                      method = 'rq',
                      quantiles = c(0.05, 0.95), col = 'blue') +
        facet_grid(Species ~ Treatment)
      
      hull <- hull %>%
        mutate(A_net = A_above - A_below)
        
      
      summary(lm(A_net  ~ Species + Treatment, hull))
      anova(lm(A_net  ~ Species + Treatment, hull))
      anova(lm(A_net  ~ Species * Treatment, hull))
      
      
      
      m0 <-  lm(A_net  ~ Species + Treatment, hull)
      cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
      cld(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
      
      pippo <-  hull %>%
        filter(Species == 'P') 
      
      anova(lm(A_net  ~ Treatment, pippo))
      m0 <- lm(A_net  ~ Treatment, pippo)
      cld(emmeans(m0, pairwise ~  'Treatment') )
      
      pippo <-  hull %>%
        filter(Species == 'J') 
      
      anova(lm(A_net  ~ Treatment, pippo))
      m0 <- lm(A_net  ~ Treatment, pippo)
      cld(emmeans(m0, pairwise ~  'Treatment') )
      
    }
    
    # 4) CORRECT ONE: quantile regression at tree level
    {
      # fits two regressions, one above (5%) and one below (95% quantile) cloud of points
      Strng_Area <-  Los_Al %>%
        filter(Time == 'Ymd') %>%
        group_by(Species, Treatment, Tree_ID) %>%  
        do(as.data.frame(lapply(
          summary(rq(Ymdreal ~ Ytf, 
                     data=.,
                     tau=c(0.05, 0.95)),se="rank"), coef))) %>%
        ungroup() %>%
        # retains only upper quantile
        select(!c('coefficients', 'lower.bd', 'upper.bd')) %>% 
        mutate(var = rep(c('int', 'slp'), times = 62)) %>%
        rename(coef_95 = 'coefficients.1',
               low_bd_95 = 'lower.bd.1',
               upp_bd_95 = 'upper.bd.1') %>% 
        pivot_wider(., names_from = 'var', 
                    values_from = c(coef_95, low_bd_95, upp_bd_95)) %>%
        left_join(Los_Al %>%
                    group_by(Species, Treatment, Tree_ID) %>%  
                    filter(Time == 'Ymd') %>% 
                    # calculates low boundary of Ytf and high boundary of Ymn distributions
                    summarise(minYtf95 = quantile(Ytf, probs = 0.05, na.rm=TRUE),
                              maxYmd95 = quantile(Ytf, probs = 0.95, na.rm=TRUE) ), 
                  ., 
                  by = c('Species', 'Treatment', 'Tree_ID')) %>%
        # estimates minimum Ymd on upper quantile regression
        mutate(minYmd_95 = coef_95_int + coef_95_slp * minYtf95) %>%
        # estimates maximum Ytf @ Ymd=0 on upper quantile regression (can be >ve or <ve)
        mutate(maxYtf_95 = (-1) * coef_95_int / coef_95_slp ) %>%
        # estimates Ymd on 1:1 line when 95% quantile regression intercepts 1:1 line
        mutate(Ymdintrcpt = coef_95_int / (1 - coef_95_slp)  ) %>%
        mutate(A_tot = minYtf95 * minYtf95 / 2 ) %>%
        # turns out this is valid independent of whether intercept is >ve or <ve
        mutate(A_above = (0 - minYmd_95) * (maxYtf_95 - minYtf95 ) / 2) %>%
        # the extra bit is either a scalene triangle or a trapezoide
        mutate(A_xtra_abv = 
                 ifelse( coef_95_int > 0, 
                         (0 - maxYmd95) * ((0-maxYtf_95)+(maxYmd95-(maxYmd95-coef_95_int)/coef_95_slp))/2, 
                         (0 - Ymdintrcpt) * maxYtf_95 / 2) 
        ) %>% 
        mutate(A_net = 
                 ifelse( coef_95_int > 0,
                         A_tot - (A_above + A_xtra_abv),
                         A_tot - (A_above - A_xtra_abv))
               )
      
      # predict 95% uncertainty from quantile regression
      {
        rqfits <- data.frame(matrix(data = NA, nrow = 10, ncol = 5))
        colnames(rqfits) <- c('Species', 'Treatment', 'fit', 
                              'minYmd_low', 'minYmd_upp')
        
        m = 1
        for(i in unique(Los_Al$Species) ) {
          
          ifelse(m == 1, n <- 1, n <- 6)
          
          for(j in unique(Los_Al$Treatment)) {
            
            pippo <- subset(Los_Al, subset= Species == i)
            pippo <- subset(pippo, subset= Treatment == j)
            pippo2 <- subset(Strng_Area, subset= Species == i)
            pippo2 <- subset(pippo2, subset= Treatment == j)
            
            rqfits[n,1] <- i
            rqfits[n,2] <- j
            # estimates uncertainty on 95% prediction quantile of minYmd_95
            rqfits[n,3:5]  <-    predict(rq(Ymdreal ~ Ytf,
                                            tau=0.95,
                                            data= pippo),
                                         interval = 'confidence',
                                         newdata = data.frame(Ytf = meanNA(pippo2$minYtf95)
                                         )
            ) / 1.96
            
            n=n+1
          }
          m=m+1
        }
        }
      
      # uncertainty in Ymd,min averaged by treatment
      Strng_Area <- left_join(Strng_Area, 
                              rqfits, by = c('Species', 'Treatment')) %>%
        mutate(A_above_upp = (0 - minYmd_upp) * (maxYtf_95 - minYtf95 ) / 2) %>%
        mutate(A_above_low = (0 - minYmd_low) * (maxYtf_95 - minYtf95 ) / 2) %>%
        mutate(Astrng = A_net) %>%
        mutate(Astrng_upp = ifelse( coef_95_int > 0,
                                    A_tot - (A_above_upp + A_xtra_abv),
                                    A_tot - (A_above_upp - A_xtra_abv)
                                    ) ) %>%
        mutate(Astrng_low = ifelse( coef_95_int > 0,
                                    A_tot - (A_above_low + A_xtra_abv),
                                    A_tot - (A_above_low - A_xtra_abv)
                                    ) ) %>%
        group_by(Species, Treatment, Tree_ID) %>% 
        summarise(Astrng = meanNA(Astrng),
                  Astrng_upp = meanNA(Astrng_upp),
                  Astrng_low = meanNA(Astrng_low)) 
      
      
      Strng_Area %>%
        mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
        mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J'))  %>%
        summarise(Astrng = mean(Astrng),
                  Astrng_upp = mean(Astrng_upp),
                  Astrng_low = mean(Astrng_low)) %>%
        ggplot(aes(y = Astrng,
                   x = Treatment),
               group = Species) + 
        geom_point(aes(color = Species),
                   position = position_dodge(width = 0.30), size = 3)  + 
        geom_errorbar(aes(ymin = Astrng_low,
                          ymax = Astrng_upp,
                          group = Species,
                          color = Species), 
                      position = position_dodge(width = 0.3),
                      width = 0.2, linewidth = 1.2) +
        theme_classic() +
        My_Theme +
        xlab('Treatments') +
        ylab('Regulation stringency area')
      
      
      summary(lm(Astrng  ~ Species + Treatment, Strng_Area))
      anova(lm(Astrng  ~ Species + Treatment, Strng_Area))
      anova(lm(Astrng  ~ Species * Treatment, Strng_Area))
      
      
      
      m0 <-  lm(Astrng  ~ Species + Treatment, Strng_Area)
      cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
      cld(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
      
      pippo <-  Strng_Area %>%
        filter(Species == 'P') 
      
      anova(lm(Astrng  ~ Treatment, pippo))
      summary(lm(Astrng  ~ Treatment, pippo))
      
      pippo <-  Strng_Area %>%
        filter(Species == 'J') 
      
      anova(lm(Astrng  ~ Treatment, pippo))
      summary(lm(Astrng  ~ Treatment, pippo))
      
    }
    
    
  }
  
    # Hydroscape area or mean DeltaY and stats
  {
    
    Hscape_Area <-  Los_Al %>%
      filter(Time == 'Ymd') %>%
      group_by(Species, Treatment, Tree_ID) %>%  
      do(as.data.frame(lapply(
        summary(rq(Ymdreal ~ Ypdreal, 
                   data=.,
                   tau=c(0.05, 0.95)),se="rank"), coef))) %>%
      ungroup() %>%
      select(!c(coefficients.1, lower.bd.1, upper.bd.1)) %>%
      mutate(var = rep(c('int', 'slp'), times = 62)) %>%
      rename(coef = 'coefficients',
             low_bd = 'lower.bd',
             upp_bd = 'upper.bd') %>% 
      pivot_wider(., names_from = 'var', values_from = c(coef, low_bd, upp_bd)) %>%
      left_join(Los_Al %>%
                  filter(Time == 'Ymd'), ., by = c('Species', 'Treatment', 'Tree_ID')) %>%
      group_by(Species, Treatment, Tree_ID) %>%  
      mutate(minYpd = quantile(Ypdreal, probs = 0.05, na.rm=TRUE))  %>%
      mutate(minYmd = coef_int + coef_slp * minYpd) %>%
      # area triangle below 1:1 line
      mutate(A_tot = minYpd^2 / 2) %>%
      # area triangle below lower quantile
      mutate(A_below = (minYpd-coef_int) * minYpd / 2) %>%
      # respective uncertainty bands
      mutate(A_below_upp = (minYpd-upp_bd_int) * minYpd / 2) %>%
      mutate(A_below_low = (minYpd-low_bd_int) * minYpd / 2) %>%
      # Hydroscape area of scalene triangle
      mutate(Astrng = A_tot - A_below) %>%
      # respective uncertainty bands
      mutate(Astrng_low = A_tot - A_below_upp) %>%
      mutate(Astrng_upp = A_tot - A_below_low) %>%
      mutate(DeltaY = Ypdreal - Ymdreal) %>%
      summarise(coef_int = mean(coef_int),
                minYpd = mean(minYpd),
                Astrng = mean(Astrng),
                A_tot  = mean(A_tot),
                A_below= mean(A_below),
                Astrng_upp = mean(Astrng_upp),
                Astrng_low = mean(Astrng_low),
                DeltaY = mean(DeltaY)) %>%
      ungroup()
    
    anova(lm(Astrng  ~ Species + Treatment, Hscape_Area))
    anova(lm(Astrng  ~ Species * Treatment, Hscape_Area))
    
    
    
    m0 <-  lm(Astrng  ~ Species + Treatment, Hscape_Area)
    cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
    cld(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
    
    pippo <-  Hscape_Area %>%
      filter(Species == 'P') 
    
    anova(lm(Astrng  ~ Treatment, pippo))
    
    pippo <-  Hscape_Area %>%
      filter(Species == 'J') 
    
    anova(lm(Astrng  ~ Treatment, pippo))

  }
  
  # comparison of Hydroscape and Stringency areas
  {
  
  AllAs <- full_join(Hscape_Area %>%
                       group_by(Species, Treatment) %>% 
                       summarise(A_Hscape = meanNA(Astrng),
                                 se_Hscape= sterr(Astrng)), 
                     Strng_Area %>% 
                       group_by(Species, Treatment) %>% 
                       summarise(A_strng = meanNA(Astrng),
                                 # from conf.int. to st.err.
                                 se_strng= meanNA(c(Astrng_upp-Astrng,Astrng-Astrng_low))/1.96,
                                 # across individuals instead.
                                 se_strng2= sterr(Astrng)
                                 )
                     )
  
  
  AllAs %>% 
    mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
    ggplot(aes(x=A_Hscape,
               y=A_strng,
               group=Species) ) +
    geom_point(aes(colour = Species,
                   shape = Treatment),
               size = 4 ) +
    geom_errorbar(aes(ymin = A_strng - se_strng2,
                      ymax = A_strng + se_strng2,
                      colour = Species) ) +
    geom_errorbarh(aes(xmin = A_Hscape - se_Hscape,
                       xmax = A_Hscape + se_Hscape,
                       colour = Species) ) +
    theme_classic() +
    My_Theme +
    xlab(expression('Hydroscape area, ' ~ A[HS] ~ ',' ~ MPa^2)) +
    ylab(expression('Stringency area,'  ~ A[S] ~ ',' ~ MPa^2))
  
 
  
  }
  
  
}

# sensitivity of forcing to sampling intensity
{
  
  # set at the beginning only
  finl_Los_Al <- NULL
  smpl_Los_Al <- NULL
  
  # sample size varies from 10 and 98 (max group size). 
  # samples are taken by species/treatment combinations
  # Each combination is sampled 1000 times
  n <- 600 # sample size/tree_ID (max is 450) (choose between n=5, 10, 25, 50, 100, 200, 400, 600)
  i <- 1
  
  prob_slct <- 0.95
  
  for(i in 1:1000) {
    
    
  smpl_Los_Al <- Los_Al %>%
    # sample size pools within and across individuals
    group_by(Species, Treatment) %>%
    sample_n(size= n, replace = TRUE) %>%
    select(!GKmaxmd) %>%
    ungroup()
  
  
  # calculates ref GKpd using quantiles over dates by chamber
  pippo <- smpl_Los_Al %>%
    filter(Time == 'Ymd') %>%
    group_by(Tree_ID) %>%
    # group_by(Species) %>%
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD,  + GKmaxpd * VPDpd
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    summarise(GKmaxmd = quantile(GKvary, probs = prob_slct, na.rm = TRUE)) %>%
    ungroup()  
  
  
  #__________________________________________________________________________  
  smpl_Los_Al <- left_join(smpl_Los_Al, pippo, by = 'Tree_ID')
  # smpl_Los_Al <- left_join(smpl_Los_Al, pippo, by = 'Species')
  #__________________________________________________________________________  
  
  
  pippo <- smpl_Los_Al %>%
    filter(Time == 'Ymd') %>%
    group_by(Tree_ID) %>%
    # group_by(Species) %>%
    mutate(Forc_VPD = ((-1) * GKmaxmd * VPD)) %>%
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    mutate(Ytf = Ysoil - GKmaxmd * VPDmd) %>% 
    mutate(Reg_tot = WP - Ytf) %>%
    mutate(pcForc_VPD = 100 * Forc_VPD /Ytf) %>%
    summarise(GKvary_smpl = quantile(GKvary, 
                                   probs = c(prob_slct), na.rm = TRUE),
              Ytf_smpl =  quantile(Ytf, 
                                probs = c(0.50), na.rm = TRUE),
              Forc_VPD_smpl =  quantile(Forc_VPD, 
                                     probs = c(0.50), na.rm = TRUE),
              PcForc_VPD_smpl =  quantile(pcForc_VPD, 
                                     probs = c(0.50), na.rm = TRUE)
    ) %>%
    ungroup()
  
  #__________________________________________________________________________  
  finl_Los_Al <- left_join(smpl_Los_Al, pippo, by = 'Tree_ID')
  # finl_Los_Al <- left_join(smpl_Los_Al, pippo, by = 'Species')
  #__________________________________________________________________________  
  
  print(i)
  
  i <- i +1
  }
  
pippo <-  finl_Los_Al %>%
    group_by(Species, Treatment) %>%
    summarise(GKvary025 = quantile(GKvary_smpl, 
                                   probs = c(0.025), na.rm = TRUE),
              GKvary050 = quantile(GKvary_smpl,
                                   probs = c(0.50), na.rm = TRUE),
              GKvary0975= quantile(GKvary_smpl,
                                   probs = c(0.975), na.rm = TRUE),
              Ytf025    = quantile(Ytf_smpl,
                                   probs = c(0.025), na.rm = TRUE),
              Ytf050    = quantile(Ytf_smpl,
                                   probs = c(0.50), na.rm = TRUE),
              Ytf0975   = quantile(Ytf_smpl,
                                   probs = c(0.975), na.rm = TRUE),
              Reg025    = quantile(Ymdreal-Ytf_smpl,
                                   probs = c(0.025), na.rm = TRUE),
              Reg050    = quantile(Ymdreal-Ytf_smpl,
                                   probs = c(0.50), na.rm = TRUE),
              Reg0975   = quantile(Ymdreal-Ytf_smpl,
                                   probs = c(0.975), na.rm = TRUE),
              Forc_VPD025=quantile(Forc_VPD_smpl,
                                   probs = c(0.025), na.rm = TRUE),
              Forc_VPD50 =quantile(Forc_VPD_smpl,
                                   probs = c(0.50), na.rm = TRUE),
              Forc_VPD0975=quantile(Forc_VPD_smpl,
                                   probs = c(0.975), na.rm = TRUE),
              PcForc_VPD025=quantile(PcForc_VPD_smpl,
                                   probs = c(0.025), na.rm = TRUE),
              PcForc_VPD50=quantile(PcForc_VPD_smpl,
                                   probs = c(0.50), na.rm = TRUE),
              PcForc_VPD0975=quantile(PcForc_VPD_smpl,
                                   probs = c(0.975), na.rm = TRUE)
    ) %>%
  mutate(sample_size = n)

# the first sample size/tree ID
# sensitivity <- pippo

# subsequent ones
sensitivity <- full_join(sensitivity, pippo)

sensitivity %>%
  mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
  mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
  ggplot(aes(x = sample_size,
             y = GKvary050,
             group=Treatment)) +
  geom_point(aes(colour = Treatment), size = 2) +
  geom_errorbar(aes(ymin = GKvary025,
                    ymax = GKvary0975,
                    colour = Treatment)) +
  theme_classic() +
  My_Theme +
  theme(legend.position="none") + 
  scale_x_continuous(trans = 'log10') +
  xlab(expression('Sample size,' ~ log[10] ~ 'n')) +
  ylab(expression(G[tot] ~ '/' ~ K[L])) +
  ylim(-2,2) +
  facet_grid(Species ~ Treatment) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold.italic"),
    strip.text.y = element_text(size = 14, face = "bold.italic") )  +
  panel_border(colour='grey85', size=1, linetype=1)


sensitivity %>%
  mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
  mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
  ggplot(aes(x = sample_size,
             y = Ytf050,
             group=Treatment)) +
  geom_point(aes(colour = Treatment), size = 2) +
  geom_errorbar(aes(ymin = Ytf025,
                    ymax = Ytf0975,
                    colour = Treatment)) +
  theme_classic() +
  My_Theme +
  theme(legend.position="none") +
  scale_x_continuous(trans = 'log10') +
  xlab(expression('Sample size,' ~ log[10] ~ 'n')) +
  ylab(expression('Total hydraulic forcing')) +
  ylim(-10,0) +
  facet_grid(Species ~ Treatment) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold.italic"),
    strip.text.y = element_text(size = 14, face = "bold.italic") )  +
  panel_border(colour='grey85', size=1, linetype=1)


sensitivity %>%
  mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
  mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
  ggplot(aes(x = sample_size,
             y = Reg050,
             group=Treatment)) +
  geom_point(aes(colour = Treatment), size = 2) +
  geom_errorbar(aes(ymin = Reg025,
                    ymax = Reg0975,
                    colour = Treatment)) +
  theme_classic() +
  My_Theme +
  theme(legend.position="none") +
  scale_x_continuous(trans = 'log10') +
  xlab(expression('Sample size,' ~ log[10] ~ 'n')) +
  ylab(expression('Total regulation,' ~ italic(S)^ws)) +
  ylim(-5, 5) +
  facet_grid(Species ~ Treatment) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold.italic"),
    strip.text.y = element_text(size = 14, face = "bold.italic") )  +
  panel_border(colour='grey85', size=1, linetype=1)


sensitivity %>%
  mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
  mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
  ggplot(aes(x = sample_size,
             y = Forc_VPD50,
             group=Treatment)) +
  geom_point(aes(colour = Treatment), size = 2) +
  geom_errorbar(aes(ymin =Forc_VPD025,
                    ymax = Forc_VPD0975,
                    colour = Treatment)) +
  theme_classic() +
  My_Theme +
  theme(legend.position="none") +
  scale_x_continuous(trans = 'log10') +
  xlab(expression('Sample size,' ~ log[10] ~ 'n')) +
  ylab(expression('VPD forcing')) +
  facet_grid(Species ~ Treatment) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold.italic"),
    strip.text.y = element_text(size = 14, face = "bold.italic") )  +
  panel_border(colour='grey85', size=1, linetype=1)


sensitivity %>%
  mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
  mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
  ggplot(aes(x = sample_size,
             y = PcForc_VPD50,
             group=Treatment)) +
  geom_point(aes(colour = Treatment), size = 2) +
  geom_errorbar(aes(ymin =PcForc_VPD025,
                    ymax = PcForc_VPD0975,
                    colour = Treatment)) +
  theme_classic() +
  My_Theme +
  theme(legend.position="none") +
  scale_x_continuous(trans = 'log10') +
  xlab(expression('Sample size,' ~ log[10] ~ 'n')) +
  ylab(expression('Percent of VPD forcing')) +
  facet_grid(Species ~ Treatment) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold.italic"),
    strip.text.y = element_text(size = 14, face = "bold.italic") )  +
  panel_border(colour='grey85', size=1, linetype=1)



}

# changes in conclusions if one estimates S at reference Ysoil and VPD
{
  
  Ypdref = -0.1
  VPDref = 0.1
  
  pippo1 <-  Los_Al %>%
    filter(Time == 'Ymd',
           Species == 'J') %>%
    group_by(Tree_ID) %>%  
    mutate(Ysoil = Ypdreal) %>%  
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    drop_na(GKvary) %>%
    split(list(.$Species)) %>%
    lapply(., function(x) { predict(lm(GKvary ~ Ypdreal * VPDmd,
                                       x[,]), 
                                    newdata = data.frame(Ypdreal = Ypdref, VPDmd = VPDref)) } ) %>%
    do.call(rbind, .) %>%
    data.frame(.) %>%
    mutate(Species = rownames(.)) %>%
    ungroup() %>%
    rename(pred_GKmaxmd = X1) %>%
    left_join(Los_Al %>%
                filter(Time == 'Ymd') %>% 
                select(!pred_GKmaxmd), 
              ., 
              by = c('Species')) %>%
    group_by(Species, Treatment, Tree_ID) %>%  
    mutate(pred_Ytf = Ysoil - pred_GKmaxmd * VPDmd) %>% 
    mutate(pred_Forc_VPD = ((-1) * pred_GKmaxmd * VPD)) %>%
    mutate(pred_Reg_tot = WP - pred_Ytf) %>%
    mutate(pcForc_VPD = 100* Forc_VPD /Ytf) %>%
    mutate(pred_pcForc_VPD = 100* pred_Forc_VPD /pred_Ytf) %>%
    ungroup()
  
  pippo2 <-  Los_Al %>%
    filter(Time == 'Ymd',
           Species == 'P') %>%
    group_by(Tree_ID) %>%  
    mutate(Ysoil = Ypdreal) %>%  
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    drop_na(GKvary) %>%
    split(list(.$Species)) %>%
    lapply(., function(x) { predict(lm(GKvary ~ Ypdreal + VPDmd,
                                       x[,]), 
                                    newdata = data.frame(Ypdreal = Ypdref, VPDmd = VPDref)) } ) %>%
    do.call(rbind, .) %>%
    data.frame(.) %>%
    mutate(Species = rownames(.)) %>%
    ungroup() %>%
    rename(pred_GKmaxmd = X1) %>%
    left_join(Los_Al %>%
                filter(Time == 'Ymd') %>% 
                select(!pred_GKmaxmd), 
              ., 
              by = c('Species')) %>%
    group_by(Species, Treatment, Tree_ID) %>%  
    mutate(pred_Ytf = Ysoil - pred_GKmaxmd * VPDmd) %>% 
    mutate(pred_Forc_VPD = ((-1) * pred_GKmaxmd * VPD)) %>%
    mutate(pred_Reg_tot = WP - pred_Ytf) %>%
    mutate(pcForc_VPD = 100* Forc_VPD /Ytf) %>%
    mutate(pred_pcForc_VPD = 100* pred_Forc_VPD /pred_Ytf) %>%
    ungroup()
  
  pippo <- full_join(pippo1, pippo2)
  
  pippo %>%
    mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
    ggplot(aes(x = GKmaxmd,
               y = pred_GKmaxmd,
               group=Treatment)) +
    geom_point(aes(colour = Treatment), size = 2) +
    geom_abline(intercept=0, slope =1, linetype = 'dashed', size=1.5) +
    theme_classic() +
    My_Theme +
    theme(legend.position="none") + 
    xlab(expression('GKmax from quantile, ' ~ G/K[max])) +
    ylab(expression('GKmax from regression, ' ~ G/K[max])) +
    facet_grid(Species ~ Treatment) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    panel_border(colour='grey85', size=1, linetype=1)  
  
  pippo %>%
    mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
    ggplot(aes(x = Ytf,
               y = pred_Ytf,
               group=Treatment)) +
    geom_point(aes(colour = Treatment), size = 2) +
    geom_abline(intercept=0, slope =1, linetype = 'dashed', size=1.5) +
    theme_classic() +
    My_Theme +
    theme(legend.position="none") + 
    xlab(expression('Forcing from quantile' ~ Psi[hf] ~ ', MPa')) +
    ylab(expression('Forcing from regression' ~ Psi[hf] ~ ', MPa')) +
    ylim(-16, 0) +
    facet_grid(Species ~ Treatment) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    panel_border(colour='grey85', size=1, linetype=1)  
  
  
  pippo %>%
    mutate(Species = fct_recode(Species, 'Pinus edulis' = 'P')) %>%
    mutate(Species = fct_recode(Species, 'Juniperus monosperma' = 'J')) %>%
    ggplot(aes(x = Reg_tot,
               y = pred_Reg_tot,
               group=Treatment)) +
    geom_point(aes(colour = Treatment), size = 2) +
    geom_abline(intercept=0, slope =1, linetype = 'dashed', size=1.5) +
    theme_classic() +
    My_Theme +
    theme(legend.position="none") + 
    xlab(expression('Regulation from quantile' ~ italic(S)[tot] ~ ', MPa')) +
    ylab(expression('Regulation from regression' ~ italic(S)[tot] ~ ', MPa')) +
    facet_grid(Species ~ Treatment) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    panel_border(colour='grey85', size=1, linetype=1)  
  
  
}

# comparison sigma Jordi versus sigma hydraulic forcing
{
  
  
  # Jordi's slopes and intercepts
  sigma_Jordi <- Los_Al %>%
    filter(Time == 'Ymd') %>%
    # group_by(Tree_ID) %>%
    group_by(Species, Treatment) %>%
    do(as.data.frame(cbind(rbind(summary(lm(Ymdreal ~ Ypdreal, data=.))$coefficients[-1,]),
                           rbind(summary(lm(Ymdreal ~ Ypdreal, data=.))$coefficients[-2,])))
    ) %>%
    select(!starts_with(c('Std. Error', 't value', 'Pr(>|t|)'))) %>%
    rename(intercept_Jordi = 'Estimate...7',
           sigma_Jordi = 'Estimate...3')
  
  
  Los_Al <- full_join(Los_Al %>%
                        filter(Time == 'Ymd'),
                      sigma_Jordi   ) 
  
  
  
  # thermodynamic forcing slopes and intercepts
  sigma_tf <- Los_Al %>%
    filter(Time == 'Ymd') %>%
    # group_by(Tree_ID) %>%
    group_by(Species, Treatment) %>%
    do(as.data.frame(cbind(rbind(summary(lm(Ymdreal ~ Ytf, data=.))$coefficients[-1,]),
                           rbind(summary(lm(Ymdreal ~ Ytf, data=.))$coefficients[-2,])))
    ) %>%
    select(!starts_with(c('Std. Error', 't value', 'Pr(>|t|)'))) %>%
    rename(intercept_tf = 'Estimate...7',
           sigma_tf = 'Estimate...3')
  
  
  Los_Al <- full_join(Los_Al %>%
                        filter(Time == 'Ymd'),
                    sigma_tf   ) 
  
  
  # figure

  tbl <- Los_Al %>% 
    group_by(Species, Treatment) %>% 
    summarise(sigma_J_mn = meanNA(sigma_Jordi),
              sigma_J_se =  sterr(sigma_Jordi),
              incpt_J_mn = meanNA(intercept_Jordi),
              incpt_J_se =  sterr(intercept_Jordi),
              sigma_tf_mn =  meanNA(sigma_tf),
              sigma_tf_se =  sterr(sigma_tf),
              incpt_tf_mn =  meanNA(intercept_tf),
              incpt_tf_se =  sterr(intercept_tf)
    ) %>%
    ungroup()
  
  tbl %>%
    ggplot(aes(x = sigma_J_mn,
               y = sigma_tf_mn),
           group = Species) + 
    geom_point(aes(color = Treatment),
               size = 3)  + 
    geom_errorbarh(aes(xmin = sigma_J_mn - 1.96 * sigma_J_se,
                       xmax = sigma_J_mn + 1.96 * sigma_J_se,
                       group = Species,
                       color = Treatment),
                  width = 0.2, size = 1.2) +
    geom_errorbar(aes(ymin = sigma_tf_mn - 1.96 * sigma_tf_se,
                      ymax = sigma_tf_mn + 1.96 * sigma_tf_se,
                      group = Species,
                      color = Treatment),
                  width = 0.01, size = 1.2) +
    geom_abline(slope=1, intercept=0) +
    theme_classic() +
    My_Theme +
    xlim(0,0.8) +
    ylim(0,0.65) +
    xlab('sigma according to Martinez-Vilalta et al. (2014)') +
    ylab(expression('sigma hydraulic forcing')) +
    annotate("text",
             x = c(0.40, 0.65),
             y = c(0.05, 0.50),
             color = "black",
             size = 8,
             label = c("P. edulis",
                       'J. monosperma') )
    
  
  
  
  
}

# comparison of 95% quantiles in forcing and regulation cross spp and treatments
{
  
   
  
  # figure
  
  tbl <- Los_Al %>% 
    filter(Time == 'Ymd') %>%
    group_by(Species, Treatment, Tree_ID) %>% 
    summarise(Forc_max = quantile(Ytf, probs = 0.05, na.rm = TRUE),
              Reg_max  = quantile(Reg_tot, probs = 0.95, na.rm = TRUE)) %>%
    ungroup()
  
  pippo <-  lm(Forc_max  ~ Species * Treatment, data = tbl)
  summary(pippo)
  anova(pippo)
  m0 <-  lm(Forc_max  ~ Species * Treatment, data = tbl) 
  cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
  cld(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
  cld(emmeans(m0, ~ Species * Treatment) )
  
  pippo <-  lm(Reg_max  ~ Species * Treatment, data = tbl)
  summary(pippo)
  anova(pippo)
  m0 <-  lm(Reg_max  ~ Species * Treatment, data = tbl) 
  cld(emmeans(m0, pairwise ~ Treatment, var = 'Species') )
  cld(emmeans(m0, pairwise ~ Species, var = 'Treatment') )
  cld(emmeans(m0, ~ Species * Treatment) )
  
  
  tbl %>%
    ggplot(aes(x = sigma_J_mn,
               y = sigma_tf_mn),
           group = Species) + 
    geom_point(aes(color = Treatment),
               size = 3)  + 
    geom_errorbarh(aes(xmin = sigma_J_mn - 1.96 * sigma_J_se,
                       xmax = sigma_J_mn + 1.96 * sigma_J_se,
                       group = Species,
                       color = Treatment),
                   width = 0.2, size = 1.2) +
    geom_errorbar(aes(ymin = sigma_tf_mn - 1.96 * sigma_tf_se,
                      ymax = sigma_tf_mn + 1.96 * sigma_tf_se,
                      group = Species,
                      color = Treatment),
                  width = 0.01, size = 1.2) +
    geom_abline(slope=1, intercept=0) +
    theme_classic() +
    My_Theme +
    xlim(0,0.8) +
    ylim(0,0.65) +
    xlab('sigma according to Martinez-Vilalta et al. (2014)') +
    ylab(expression('sigma hydraulic forcing')) +
    annotate("text",
             x = c(0.40, 0.65),
             y = c(0.05, 0.50),
             color = "black",
             size = 8,
             label = c("P. edulis",
                       'J. monosperma') )
  
  
  
  
  
}



#### Drake E tereticornis dataset ####

  
# this is the E tereticornis experiment of 2014 (combined warming * drought)

# met data
{
  
  setwd("D:/Dropbox/old stuff/old research/old grants/Australia UWS phloem project/WTC results")
  
  # this can be improved by using actual chamber-specific data
  file <- 'WTC3_Enviro_Master.csv'
  met_Drake <- read.csv(file, sep = ",", header = TRUE)
  
  
  # calculate saturated VP using formulas in Yuan et al 2019
  Pmsl <- 1013.25 # in hPa (pressure mean sea level)
  Z <- 95 # altitude (m)
  
  # Saturation vapour pressure
  met_Drake$Pmst <- with(met_Drake,
                         Pmsl * ( (AirTC_Avg + 273.16) / ( (AirTC_Avg + 273.16) + 0.0065 * Z) ) ^ 5.625 )
  met_Drake$fw <- with(met_Drake,
                       (1 + 7 * 10^(-4) + 3.46 * 10 ^(-6) * Pmst ) )
  
  met_Drake$SatVap <- with(met_Drake,
                           6.112 * fw * exp( (17.67 * AirTC_Avg) / (AirTC_Avg + 243.5)) )
  
  # vpd = esat(100-RH)/100
  met_Drake$VPD <- with(met_Drake, 
                        SatVap /10 *(100- RH ) / 100 ) # divide by 10 to go from hPa to in kPa
  
  # first split DateTime into date and hours
  met_Drake <-   met_Drake %>% 
    mutate(
      date = strsplit(TIMESTAMP, " ")  %>% 
        as.data.frame() %>% 
        t %>% 
        data.frame(stringsAsFactors = F) %>% 
        pull(1),
      Hour = strsplit(TIMESTAMP, " ")  %>% 
        as.data.frame() %>% 
        t %>% 
        data.frame(stringsAsFactors = F) %>% 
        pull(2)
    ) %>%
    dplyr::select(c(date, Hour, AirTC_Avg, RH, PPFD_Avg, VPD)) %>%
    mutate(date = parse_date_time(date, orders = 'dmy'))
  
  # then split hours only into numeric Hour and average data for the fractions
  met_Drake <-   met_Drake %>% 
    mutate(
      Hour = strsplit(Hour, ":")  %>%
        as.data.frame() %>%
        t %>%
        data.frame(stringsAsFactors = F) %>%
        pull(1)
    ) %>%
    group_by(date, Hour) %>%
    reframe(AirTC_Avg, RH, PPFD_Avg, VPD) %>%
    filter(Hour == '04' | Hour == '13') %>%
    mutate(date = parse_date_time(date, orders = 'ymd'))
  
  met_Drake <- met_Drake %>%
    select( ! c(AirTC_Avg, RH, PPFD_Avg)) %>%
    mutate(chamber = rep.int(1:12, times = length(met_Drake$Hour)/12)) %>%
    mutate(chamber = as.character(chamber)) %>%
    mutate(chamber = ifelse(chamber=='10' | chamber=='11' | chamber=='12',
                            paste('C', chamber, sep =''),
                            paste('C0', chamber, sep =''))
           ) %>%
    group_by(date, chamber) %>%
    pivot_wider(., names_from = 'Hour', values_from = 'VPD') %>%
    mutate(VPDpd = `04`,
           VPDmd = `13`) %>%
    pivot_longer(., cols = c(`04`, `13`), names_to = 'Hour', values_to = 'VPD')
  
  
}

# add file with actual VPD data by chamber
{

  setwd("D:/Dropbox/My Documents/science/papers/my papers/Homoiohydry paper/data/Drake tereticornis")

  # this improves by using actual chamber-specific data
  file <- 'WTC_daily_VPD_by_chamber_20130914-20140526.csv'
  VPD_corr_Drake <- read.csv(file, sep = ",", header = TRUE)
  
  
  # first split DateTime into date and hours
  VPD_corr_Drake <-   VPD_corr_Drake %>% 
    mutate(
      date = strsplit(DateTime, " ")  %>% 
        as.data.frame() %>% 
        t %>% 
        data.frame(stringsAsFactors = F) %>% 
        pull(1),
      Hour = strsplit(DateTime, " ")  %>% 
        as.data.frame() %>% 
        t %>% 
        data.frame(stringsAsFactors = F) %>% 
        pull(2)
    ) %>%
    dplyr::select(!Tair_al) %>%
    mutate(date = parse_date_time(date, orders = 'dmy'))
  
  # then split hours only into numeric Hour and average data for the fractions
  VPD_corr_Drake <-   VPD_corr_Drake %>% 
    mutate(
      Hour = strsplit(Hour, ":")  %>%
        as.data.frame() %>%
        t %>%
        data.frame(stringsAsFactors = F) %>%
        pull(1)
    ) %>%
    group_by(chamber, T_treatment, Water_treatment, date, Hour) %>%
    rename(VPD = VPDair) %>%
    reframe(VPD) %>%
    ungroup() %>%
    filter(Hour == '04' | Hour == '13') %>%
    mutate(date = parse_date_time(date, orders = 'ymd'))
  
  VPD_corr_Drake <- VPD_corr_Drake %>%
    group_by(chamber, T_treatment, Water_treatment, date, Hour) %>%
    pivot_wider(., names_from = 'Hour', values_from = 'VPD') %>%
    mutate(VPDpd = `04`,
           VPDmd = `13`) %>%
    pivot_longer(., cols = c(`04`, `13`), names_to = 'Hour', values_to = 'VPD') %>%
    mutate(Hour = case_when(Hour == '04' ~ 'predawn',
                            Hour == '13' ~ 'midday')) %>%
    rename(Time = Hour)
    
  

  
}

# WP data
{
  
  
  setwd("D:/Dropbox/My Documents/science/papers/my papers/Homoiohydry paper/data/Drake tereticornis")
  file <- 'WTC_TEMP_CM_WATERPOTENTIAL-PREDAWN-MIDDAY_20130515-20140424_L2.csv'
  WP_Drake <- read.csv(file, sep = ",", header = TRUE)
  
  
  
  
  # assumes pre-dawn at 4:00 and midday at 13:00
  WP_Drake_long <- WP_Drake %>%
    mutate(Ypdreal = predawn / (-10) ) %>%
    mutate(Ymdreal = midday / (-10) ) %>%
    gather(., Time, WP, c("predawn", "midday")) %>%
    mutate(Hour = case_when(Time == 'predawn' ~ '04',
                            Time == 'midday' ~ '13')) %>%
    drop_na(Ypdreal, Ymdreal) %>%
    select(! stem ) %>%
    mutate(Hour = as.character(Hour)) %>%
    mutate(date = parse_date_time(date, orders = 'dmy')) %>%
    mutate(WP = WP / (-10) )
  
  
  # # this shows that controls became treatments through the study
  # WP_Drake_long  %>%
  #   filter(Time == 'midday') %>%
  #   ggplot(aes(x = date,
  #              y = Ymdreal)) +
  #   geom_point(size=3) +
  #   theme_classic() +
  #   My_Theme +
  #   labs(title = "WP regulation Eucalyptus Australia") +
  #   xlab('cate') +
  #   ylab('WP') +
  #   facet_wrap(vars(chamber))

  
}

# merge all data
{
  
  
  Drake <- left_join(WP_Drake_long, met_Drake, by= c('date', 'Hour', 'chamber')) %>%
    drop_na(VPD) %>%
    group_by(date, chamber, Time) %>%
    summarise_if(., is.numeric, mean, na.rm = TRUE) %>%
    ungroup()
  
  
  
  # recover treatment info
  Drake <- left_join(Drake, WP_Drake_long %>%
                       group_by(date, chamber, Time) %>%
                       summarise_if(., is.character, unique) %>%
                       select(! c(Hour)) %>%
                       ungroup(), 
                     by= c('date', 'Time', 'chamber'))
  
  
  # final VPD correction for heated chambers

  
  pippo <- VPD_corr_Drake %>%
    filter(Time == 'midday')
  
  # anova(lmer(VPDpd ~ T_treatment + Water_treatment + (1|chamber) + (1|date), data = VPD_corr_Drake))
  # summary(lmer(VPDpd ~ T_treatment + Water_treatment + (1|chamber) + (1|date), data = VPD_corr_Drake))
  m0 <-  lmer(VPD ~ T_treatment + (1|chamber), data = pippo)
  
Drake <-  Drake %>%
    rowwise() %>%
    mutate(new_var = 
             ifelse(Time == "predawn", VPDpd,
                    ifelse(T_treatment == "elevated",
                           VPDmd + fixef(m0)['T_treatmentelevated'], VPDmd)
                    )) %>% 
  ungroup()

Drake <- Drake %>%
  select(!VPD) %>%
  rename(VPD = new_var) %>%
  select(!c(id, leaf))

predawn <- Drake %>%
  filter(Time == 'predawn') %>%
  mutate(VPDpd = VPD)

midday <- Drake %>%
  filter(Time == 'midday') %>%
  mutate(VPDmd = VPD)

Drake <- full_join(predawn, midday)
  
   
  # calculates ref GKpd between first and second data and averages over chambers
  # GKmaxpd  <- Drake %>% 
  #   filter(Time == 'predawn') %>% 
  #   arrange(date) %>%
  #   # Group(chamber) gets one set of values for all IDs for first date at Ypd 
  #   group_by(chamber) %>%
  #   mutate(GKpd = (first(Ypdreal)- nth(Ypdreal, n=2)) /
  #            (nth(VPDpd, n=2) - first(VPDpd))) %>% 
  #   ungroup() %>%
  #   # uses same criterion as for midday GK (slightly more conservative to reduce outliers)
  #   summarise(GKpdmean = quantile(GKpd, probs = 0.50, na.rm = TRUE)) %>%
  #   pull(GKpdmean) %>%
  #   as.numeric(.)
  
  
  # Drake <- cbind(Drake, GKmaxpd)
  

  # calculates ref GKpd using quantiles over dates by chamber
  pippo <- Drake %>%
    filter(Time == 'midday') %>%
    #___________________________________________________________________________
    #_CRUCIAL LINE
    # group_by(chamber) %>%
    # Group(chamber) gets one set of values for all IDs
    #___________________________________________________________________________
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD,  + GKmaxpd * VPDpd
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    summarise(GKmaxmd = quantile(GKvary, probs = 0.95, na.rm = TRUE)) %>%
    ungroup()  
  
  
  #___________________________________________________________________________
  #_CRUCIAL LINE

  # Drake <- left_join(Drake, pippo, by = 'chamber')
  Drake <- cbind(Drake, pippo)
  #___________________________________________________________________________
  
  
  pippo <- Drake %>%
    #___________________________________________________________________________
    #_CRUCIAL LINE
    group_by(chamber) %>%
    #___________________________________________________________________________
    select(c(chamber, T_treatment, Water_treatment, Ymdreal, Ypdreal, VPDmd, VPDpd, GKmaxmd)) %>%
    mutate(Ysoil = Ypdreal) %>%  
    mutate(GKvary = (Ysoil - Ymdreal) /VPDmd) %>%
    # rather liberal range around GKmaxmd to get better stats
    filter(GKvary < 1.2 * GKmaxmd &
           GKvary > 0.8 * GKmaxmd) %>%
    summarise(VPDmdref = quantile(VPDmd, probs = 0.05, na.rm = TRUE),
              Ypdref = quantile(Ypdreal, probs = 0.95, na.rm = TRUE)) %>%
    ungroup()
  
  #___________________________________________________________________________
  #_CRUCIAL LINE
  # pippo <- pippo %>%
    # add_row(chamber = 'C05', VPDmdref = mean(c(0.158, 0.970)), Ypdref = mean(c(-0.474, -0.24))) %>%
    # add_row(chamber = 'C09', VPDmdref = mean(c(0.158, 0.970)), Ypdref = mean(c(-0.474, -0.24)))
  #___________________________________________________________________________

    
  #___________________________________________________________________________
  #_CRUCIAL LINE
  Drake <- left_join(Drake, pippo, by = c('chamber'))
  # Drake <- cbind(Drake, pippo)
  #___________________________________________________________________________

  
  
  predawn <- Drake %>%
    filter(Time == 'predawn') %>%
    group_by(chamber) %>%
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD,  + GKmaxpd * VPDpd
    mutate(Ytf = Ysoil) %>% 
    mutate(Reg_tot = WP - Ytf)%>% 
    mutate(sequence = 1) %>%
    ungroup() 
  
  midday <- Drake %>%
    filter(Time == 'midday') %>%
    group_by(chamber) %>%
    mutate(Ysoil = Ypdreal) %>%  # correction for nighttime VPD,  + GKmaxpd * VPDpd
    mutate(Ytf = Ysoil- GKmaxmd * VPDmd) %>% 
    mutate(Reg_tot = WP - Ytf) %>%
    mutate(sequence = 2) %>%
    ungroup()  
  
  Drake <- full_join(predawn, midday) %>%
    arrange(sequence, .by_group=TRUE)
  
  
  
  Drake <- Drake %>%
    mutate(Forc_Ysoil = Ysoil) %>%
    mutate(Forc_VPD = ((-1) * GKmaxmd * VPD)) %>%
    mutate(Reg_Ysoil = Ymdreal - Ysoil) %>%
    mutate(Reg_VPD = Ymdreal - ((-1) * GKmaxmd * VPD)) %>%
    ungroup()
  
  pippo <- Drake %>%
    filter(Time=='midday') %>%
    split(list(.$T_treatment, .$Water_treatment)) %>%
    lapply(., function(x) { lm(VPDmd ~ Ypdreal,
                               x[,])$coefficients } ) %>%
    do.call(rbind, .) %>%
    data.frame(.) %>%
    mutate(W_T_Treat = rownames(.)) %>%
    mutate(    T_treatment = unlist(strsplit(W_T_Treat, split = '[.]'))[c(1,3,5,7)],
           Water_treatment = unlist(strsplit(W_T_Treat, split = '[.]'))[c(2,4,6,8)]) %>%
    select(!c(X.Intercept., W_T_Treat)) %>%
    relocate(c(T_treatment, Water_treatment)) %>%
    rename(dVPDdYpd = Ypdreal) %>%
    ungroup()
  
  # reset slopes to zero for all eucs since P always >> 0.05
  pippo <- pippo %>% 
    mutate(dVPDdYpd = 0)
  
  
  
  Drake <- left_join(Drake, pippo, by = c('T_treatment', 'Water_treatment')) %>%
    mutate(Forc_VPD_net = Forc_VPD + (+1) * GKmaxmd * (Ysoil-Ypdref) * dVPDdYpd) %>%
    mutate(Reg_VPD_net = Ymdreal - Forc_VPD_net)
  
  
  pippo <- Drake %>%
    filter(Time=='midday') %>%
    split(list(.$T_treatment, .$Water_treatment)) %>%
    lapply(., function(x) { lm(Ypdreal ~ VPDmd,
                               x[,])$coefficients } ) %>%
    do.call(rbind, .) %>%
    data.frame(.) %>%
    mutate(W_T_Treat = rownames(.)) %>%
    mutate(    T_treatment = unlist(strsplit(W_T_Treat, split = '[.]'))[c(1,3,5,7)],
           Water_treatment = unlist(strsplit(W_T_Treat, split = '[.]'))[c(2,4,6,8)]) %>%
    select(!c(X.Intercept., W_T_Treat)) %>%
    relocate(c(T_treatment, Water_treatment)) %>%
    rename(dYpddVPD = VPDmd) %>%
    ungroup()
  
  # reset slopes to zero for all eucs since P always >> 0.05
  pippo <- pippo %>%
    mutate(dYpddVPD = 0)
  
  
  
  Drake <- left_join(Drake, pippo, by = c('T_treatment', 'Water_treatment')) %>%
    mutate(Forc_Ysoil_net = Ysoil - (+1) * (VPDmd-VPDmdref) * dYpddVPD) %>%
    mutate(Forc_cov_net = (VPDmd-VPDmdref) * dYpddVPD - (GKmaxmd * (Ysoil-Ypdref) * dVPDdYpd)) %>%
    mutate(Reg_Ysoil_net = Ymdreal - Forc_Ysoil_net)  %>%
    mutate(Reg_cov_Ysoil_VPD = Ymdreal - Forc_cov_net)
  
  
  
  }

# total forcing
{

  #_______________________________________________________________________________
  # FIGURE 3B
  
  ann_text1 <- data.frame(Ytf = c(-8), 
                          WP  = c(-2.6),
                          lab = "Text",
                          T_treatment = factor('elevated T',
                                               levels = c("ambient T","elevated T")),
                          Water_treatment = factor('drydown',
                                                   levels = c('control', 'drydown')) )
  
  ann_text2 <- data.frame(lab = "Text",
                          T_treatment = factor('elevated T',
                                               levels = c("ambient T","elevated T")),
                          Water_treatment = factor('drydown',
                                                   levels = c('control', 'drydown')) )
  
  Drake %>%
    mutate(T_treatment = fct_recode(T_treatment, "ambient T" = "ambient")) %>%
    mutate(T_treatment = fct_recode(T_treatment, "elevated T" = "elevated")) %>%
    filter(Time == 'midday') %>%
    mutate(InvYtf = (-1) * Ytf) %>%
    mutate(InvWP = (-1) * WP) %>%
    mutate(T_water_Treat = as.factor(paste(T_treatment, Water_treatment))) %>% 
    group_by(T_treatment, Water_treatment) %>%
    mutate(InvRegrWP = predict(lm(log(InvWP) ~ log(InvYtf)))) %>%
    ungroup() %>%
    ggplot(aes(x = Ytf,
               y = WP) ) +
    geom_point(size=3, aes(color = chamber) ) +
    geom_point(size=2, color = 'black',
               data = function(x) {x[x$Time %in% c("predawn"), ]} ) +
    geom_abline(intercept=0, slope = 1, linetype = 'dashed', size=1.0) +
    geom_quantile(formula =  y ~ x,
                  method = 'rq',
                  quantiles = c(0.95), col = 'black', size = 1.5) +
    geom_abline(intercept=0, slope = 1, linetype = 'dotted', size=1.5) +
    geom_line(aes(x=Ytf,
                  y=(-1)*exp(InvRegrWP),
                  color = T_water_Treat),
              linetype = 'dashed', size = 1.5) +
    theme_classic() +
    My_Theme +
    labs(title = "Water status regulation E. tereticornis, Australia") +
    labs(x = expression(Psi[hf] ~ ', MPa'),
         y = expression(Psi[md] ~ ', MPa')) +
    xlim(-10,0) +
    ylim(-3,0) +
    theme(legend.position="none") +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    facet_grid(T_treatment ~ Water_treatment) +
    geom_text(data = ann_text1,
              label = expression(A[s]),
              size = 6,
              fontface = "italic",
              color = 'black') +
    geom_segment(data = ann_text2,
                 x    = c(-7.5),
                 xend = c(-6),
                 y    = c(-2.6),
                 yend = c(-2.6),
                 color = c("black"), 
                 linewidth = 1.5, linetype = c(1)) +
    panel_border(colour='grey85', size=1, linetype=1)
  
  

    }

# Stringency of regulation
{
  
  Drake %>%
    mutate(T_treatment = fct_recode(T_treatment, "ambient T" = "ambient")) %>%
    mutate(T_treatment = fct_recode(T_treatment, "elevated T" = "elevated")) %>%
    filter(Time=='midday') %>%
    mutate(Year = ifelse(date < "2014-02-17 UTC", "yr_1", "yr_2") ) %>%
    filter(Year == 'yr_2') %>%
    ggplot(aes(x = Ytf, 
               y = Reg_VPD)) +
    geom_point(aes(x = Ytf,
                   y = Reg_tot),
               size = 2,
               shape = 15,
               color = 'black',
               alpha = 1/1) +
    geom_point(size = 3,
               shape = 16,
               color = 'red',
               alpha = 1/1) +
    scale_color_viridis_c(option = "plasma") +
    geom_abline(intercept=0, slope = 0, linetype = 1, size=1.5) +
    geom_smooth(method='lm',
                formula= y~x,
                # se = FALSE,
                aes(x = Ytf,
                    y = Reg_tot),
                color = 'black',
                size = 1) +
    geom_smooth(method='lm',
                formula= y~x,
                se = FALSE,
                color = 'red',
                size = 1) +
    theme_classic() +
    My_Theme +
    theme(plot.margin = rep(unit(0.5, "cm"), times = 4) ) +
    labs(title = "Regulation against VPD, WTC Australia") +
    xlab(expression(Psi[soil] ~ ', MPa')) +
    labs(y = expression(S[ws] ~ ',' ~ S[VPD] ~ ', MPa') ) +
    facet_grid(T_treatment ~ Water_treatment) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    panel_border(colour='grey85', size=1, linetype=1)
}

# partitioning of forcing 
{
  
#______________________________________________________________________
  # FIGURE 5B partitioning of NET forcing
  # forcing by VPD and Ysoil accounting for covariation (simple formula)
  
  ann_text1 <- data.frame(Ytf = -2, 
                          Forc_VPD_net = c(-3, -4.2),
                          lab = "Text",
                          T_treatment = factor('elevated T',
                                       levels = c("ambient T","elevated T")),
                          Water_treatment = factor('drydown',
                                             levels = c('control', 'drydown'))
  )
  
  ann_text2 <- data.frame(Ytf = -2, 
                          Forc_VPD_net = c(-3.6),
                          lab = "Text",
                          T_treatment = factor('elevated T',
                                               levels = c("ambient T","elevated T")),
                          Water_treatment = factor('drydown',
                                                   levels = c('control', 'drydown'))
  )
  
  Drake %>%
    mutate(T_treatment = fct_recode(T_treatment, "ambient T" = "ambient")) %>%
    mutate(T_treatment = fct_recode(T_treatment, "elevated T" = "elevated")) %>%
    filter(Time=='midday') %>%
    ggplot(aes(x = Ytf, 
               y = Forc_VPD_net)) +
    scale_color_viridis_c(option = "plasma") +
    geom_point(aes(x = Ytf, 
                   y = Forc_VPD_net),
               size = 3,
               shape = 16,
               color = 'red',
               alpha = 1/2) +
    geom_smooth(method='lm',
                formula= y~x,
                se = TRUE,
                color = 'red',
                size = 1) +
    geom_point(aes(x = Ytf,
                   y = Forc_Ysoil_net),
               size = 2,
               shape = 15,
               color = 'black',
               alpha = 1/2) +
    geom_smooth(aes(x = Ytf,
                    y = Forc_Ysoil_net),
                method='lm',
                formula= y~x,
                se = TRUE,
                color = 'black',
                size = 1) +
    geom_point(aes(x = Ytf,
                   y = Forc_cov_net),
               size = 2,
               shape = 15,
               color = 'blue',
               alpha = 1/2) +
    geom_smooth(aes(x = Ytf,
                    y = Forc_cov_net),
                method='lm',
                formula= y~x,
                se = TRUE,
                color = 'blue',
                size = 1) +
    geom_line(aes(x = Ytf,
                  y = Ytf),
              linetype = 'dashed', size = 1.5) +
    theme_classic() +
    My_Theme +
    theme(plot.margin = rep(unit(0.5, "cm"), times = 4) ) +
    labs(title = "Eucalyptus tereticornis, Australia") +
    xlab('Total forcing, MPa') +
    labs(y = expression('Net forcing by VPD and' ~ Psi[soil] ~ ', MPa') ) +
    facet_grid(T_treatment ~ Water_treatment) +
    # xlim(-1,3) +
    ylim(-6,0) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    geom_text(data = ann_text1,label = c("Covar", 'VPD'),
              size = 6,
              fontface = "bold",
              color = c('blue', 'red')) +
    geom_text(data = ann_text2,label = c(expression(bold(Psi[soil]))),
              size = 6,
              fontface = "bold",
              color = c('black')) +
    panel_border(colour='grey85', size=1, linetype=1)
  
}

# Isohydricity plots (MV style)
{
  
  ann_text1 <- data.frame(Ypdreal = c(-2, -2), 
                          Ymdreal = c(-2.2, -2.5),
                          lab = "Text",
                          T_treatment = factor('elevated T',
                                       levels = c("ambient T","elevated T")),
                          Water_treatment = factor('drydown',
                                             levels = c("control","drydown")) )
  
  ann_text2 <- data.frame(lab = rep("Text",2),
                          T_treatment = factor('elevated T',
                                       levels = c("ambient T","elevated T")),
                          Water_treatment = factor('drydown',
                                             levels = c("control","drydown")) )
  
  
  
  Drake  %>%
    mutate(T_treatment = fct_recode(T_treatment, "ambient T" = "ambient")) %>%
    mutate(T_treatment = fct_recode(T_treatment, "elevated T" = "elevated")) %>%
    ggplot(aes(x = Ypdreal,
               y = Ymdreal)) +
    geom_point(size=3, color = 'red') +
    geom_smooth(method = 'lm',
                formula =  y ~ x,
                se = FALSE,
                linetype = 'dashed', size = 1.5,
                color = 'black') +
    geom_quantile(formula =  y ~ x,
                  method = 'rq',
                  quantiles = c(0.05), col = 'black', size = 1.5) +
    geom_abline(intercept=0, slope = 1, linetype = 'dotted', size=1.5) +
    xlim(-3, 0) +
    ylim(-3, 0) +
    theme_classic() +
    My_Theme +
    labs(title = "Water status regulation, E. tereticornis, Australia") +
    xlab(expression(Psi[pd])) +
    ylab(expression(Psi[md])) +
    facet_grid(T_treatment ~ Water_treatment) +
    geom_text(data = ann_text1,
              label = c(expression(sigma), 'HS'),
              size = 6,
              fontface = "bold",
              color = 'black') +
    geom_segment(data = ann_text2,
                 x    = c(-1.8, -1.8),
                 xend = c(-1.2, -1.2),
                 y    = c(-2.2, -2.5),
                 yend = c(-2.2, -2.5),
                 color = c("black","black"),
                 linewidth = 1.5, linetype = c(2,1)) +
    theme(
      strip.text.x = element_text(size = 14, face = "bold.italic"),
      strip.text.y = element_text(size = 14, face = "bold.italic") ) +
    panel_border(colour='grey85', size=1, linetype=1)
  
  

}

# stats of regulation
{

  pippo <- Drake %>%
    mutate(Year = ifelse(date < "2014-02-17 UTC", "yr_1", "yr_2") ) %>%
    filter(Time=='midday') %>%
    filter(Year == 'yr_2') %>%
    mutate(pForc_VPD = Forc_VPD / Ytf) %>%
    mutate(Treatment = paste(T_treatment, Water_treatment, sep='_')) 
  # %>%
    # mutate(logReg_tot = log(Reg_tot)) %>%
    # filter_at(vars(logReg_tot), all_vars(!is.infinite(.))) %>%
    # mutate(logYmd = log(-Ymdreal)) %>%
    # mutate(logYtf = log(-Ytf)) %>%
    # filter_at(vars(logYmd), all_vars(!is.infinite(.))) %>%
    # drop_na(logYmd, logYtf, logReg_tot)
  
  # partitioning of forcing
  {
    
    pippo <- Drake %>%
      filter(Time=='midday') %>%
      mutate(Year = ifelse(date < "2014-02-17 UTC", "yr_1", "yr_2") ) %>%
      filter(Year == 'yr_2') %>%
      mutate(pForc_VPD = Forc_VPD / Ytf) %>%
      mutate(Treatment = paste(T_treatment, Water_treatment, sep='_')) %>%
      drop_na(Ypdreal, Ytf)
    
    # proportion of forcing by VPD    
    summary(lmer(pForc_VPD ~ Ytf + Treatment + (1|chamber), data = pippo))
    anova(lmer(pForc_VPD ~ Ytf + Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(pForc_VPD ~ Ytf + Treatment + (1|chamber), data = pippo))
    
    
    m0 <-  lmer(pForc_VPD ~ Ytf + Treatment + (1|chamber), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ytf') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    
                 
    summary(lmer(Ytf ~ Forc_VPD * Treatment + (1|chamber), data = pippo))
    anova(lmer(Ytf ~ Forc_VPD * Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(Ytf ~ Forc_VPD * Treatment + (1|chamber), data = pippo))
    
    
    
    }
  
  # response of GKratio to treatments
  {
    
    
    pippo <- Drake %>%
      filter(Time=='midday') %>%
      mutate(Year = ifelse(date < "2014-02-17 UTC", "yr_1", "yr_2") ) %>%
      filter(Year == 'yr_2') %>%
      mutate(Treatment = paste(T_treatment, Water_treatment, sep='_')) %>%
      select(chamber, Treatment, GKmaxmd) %>%
      group_by(chamber, Treatment) %>%
      summarise(GKmaxmd = mean(GKmaxmd)) %>%
      ungroup()
    
    summary(lm(GKmaxmd ~ Treatment, data = pippo))
    anova(lm(GKmaxmd ~ Treatment, data = pippo))
    
    
    }
  
  # sensitivity of Ymd to Ypd alone
  {
    
    summary(lmer(Ymdreal ~ Ypdreal + Treatment + (1|chamber), data = pippo))
    anova(lmer(Ymdreal ~ Ypdreal + Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ypdreal + Treatment + (1|chamber), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Ypdreal + Treatment + (1|chamber), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
 
    
    summary(lmer(log(-Ymdreal) ~ log(-Ypdreal) + Treatment + (1|chamber), data = pippo))
    anova(lmer(log(-Ymdreal) ~ log(-Ypdreal) + Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(log(-Ymdreal) ~ log(-Ypdreal) + Treatment + (1|chamber), data = pippo))
    
    
  }
  
  # sensitivity of DeltaY to Ypd alone
  {
    
    # separate by species
    pippo <- Drake %>%
      mutate(DeltaY = Ypdreal - Ymdreal) %>% 
      filter(Time=='midday') %>%
      mutate(Year = ifelse(date < "2014-02-17 UTC", "yr_1", "yr_2") ) %>%
      filter(Year == 'yr_2') %>%
      mutate(Treatment = paste(T_treatment, Water_treatment, sep='_')) %>%
      drop_na(Ypdreal, DeltaY)
    
    
    summary(lmer(DeltaY ~ Treatment + (1|chamber), data = pippo))
    anova(lmer(DeltaY ~ Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(DeltaY ~ Treatment + (1|chamber), data = pippo))
    
    m0 <-  lmer(DeltaY ~ Ypdreal + Treatment + (1|chamber), data = pippo)
    emtrends(m0, var = 'Ypdreal', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ypdreal') )
    emmip(m0, Treatment ~ Ypdreal, cov.reduce = range)
    check_model(m0)
    
    summary(lmer(log(DeltaY) ~ Treatment + (1|chamber), data = pippo))
    anova(lmer(log(DeltaY) ~ Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(log(DeltaY) ~ Treatment + (1|chamber), data = pippo))
    
    
  }
  
  
  # sensitivity of Ymd to Yhf
  {
    
    summary(lmer(Ymdreal ~ Ytf + Treatment + (1|chamber), data = pippo))
    anova(lmer(Ymdreal ~ Ytf + Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(Ymdreal ~ Ytf + Treatment + (1|chamber), data = pippo))
    
    m0 <-  lmer(Ymdreal ~ Ytf * Treatment + (1|chamber), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ymd') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    check_model(m0)
    
    
    summary(lmer(log(-Ymdreal) ~ log(-Ytf) + Treatment + (1|chamber), data = pippo))
    anova(lmer(log(-Ymdreal) ~ log(-Ytf) + Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(log(-Ymdreal) ~ log(-Ytf) + Treatment + (1|chamber), data = pippo))
    
    
  }
  
  # sensitivity of Reg_tot
  
  {
    
    summary(lmer(log(Reg_tot) ~ T_treatment + Water_treatment + (1|chamber), data = pippo))
    anova(lmer(log(Reg_tot) ~ T_treatment + Water_treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(log(Reg_tot) ~ T_treatment + Water_treatment + (1|chamber), data = pippo))
    
    
    summary(lmer(log(Reg_tot) ~ Treatment + (1|chamber), data = pippo))
    anova(lmer(log(Reg_tot) ~ Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(log(Reg_tot) ~ Treatment + (1|chamber), data = pippo))
    
    summary(lmer(Reg_tot ~ Ytf + Treatment + (1|chamber), data = pippo))
    anova(lmer(Reg_tot ~ Ytf + Treatment + (1|chamber), data = pippo))
    r.squaredGLMM(lmer(Reg_tot ~ Ytf + Treatment + (1|chamber), data = pippo))
    
    
    m0 <-  lmer(Reg_tot ~ Ytf * Treatment + (1|chamber), data = pippo)
    emtrends(m0, var = 'Ytf', pairwise ~ Treatment)
    CLD(emmeans(m0, pairwise ~ Treatment, var = 'Ymd') )
    emmip(m0, Treatment ~ Ytf, cov.reduce = range)
    check_model(m0)
    
    
    
  }

  
}


#### theoretical figure 1 ####
Ytf_theor1 <- seq(from=-0.1, to=-4, length.out=40)
Ytf_theor2 <- seq(from=-0.1, to=-0.1, length.out=40)
Ylf_theor1 <- seq(from=-0.1, to=-4, length.out=40)
Ylf_theor2 <- seq(from=-0.1, to=--0.1, length.out=40)
Ylf_theor3 <- seq(from=-0.1, to=-4, length.out=40)

theor <- data.frame(Ytf_theor1, Ytf_theor2, Ylf_theor1, Ylf_theor2, Ylf_theor3)

p1 <-ggplot(theor, aes(x = Ytf_theor1,
                  y = Ylf_theor1)) +
  geom_line(linewidth = 1, linetype = 2) +
  annotate(geom = "segment", 
           x    = c(-0.5, -1.5, -2.3, -3.5, 0),
           xend = c(-0.5, -1.5, -2.3, -3.5, 0),
           y    = c(-4.0, -4.0, -4.0, -4.0, 0),
           yend = c(-0.5, -1.5, -2.5, -3.5, -4),
           color = c("black","black","black","black","black"), 
           linewidth = 1, linetype = 2) +
  annotate(geom = "segment", 
           x    = c(-2.3, -0.5, -1.5, -3.5),
           xend = c(-2.3, -2.3, -3.5, -3.5),
           y    = c(-2.3, -0.5, -1.6, -3.6),
           yend = c(-1.9, -2.3, -3.6, -1.9),
           color = c("blue","blue","red","red"), 
           linetype = 1, size = 2, 
           arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", 
           x = c(-1.2, -2.2, -3.56, -0.85, -3.42, -2.22, -0.08),
           y = c(-1.6, -2.6, -3.7, -0.5, -1.7, -1.7, -0.05),
           color = "black", 
           size = 8, 
           label = c('B (pd)', 'C','E', 'A (pd)', 'F (md)', 'D(md)', 'O') ) +
  annotate("text", 
           x = c(-1.4, -2.2, -0.4, -3.4, -0.1),
           y = c(-3.9, -3.9, -3.9, -3.9, -3.9),
           color = "black", 
           size = 8, 
           label = c('B*minute', 'C*minute', 'A*minute', 'E*minute', 'O*minute'), 
           parse = TRUE ) +
  annotate("text", 
           x = c(0.05),
           y = c(0.1),
           color = "white", 
           size = 8, 
           label = c(''), 
           parse = TRUE ) +
  ylim(-4, 0.05) +
  theme_classic() +
  My_Theme +
  labs(title = "Perfecly isohydric at seasonal scale") +
  xlab('Hydraulic forcing, MPa') +
  ylab('Leaf water potential, MPa') + 
  coord_cartesian(expand = FALSE) # eliminates padding of X and Y axes

p2 <-ggplot(theor, aes(x = Ytf_theor1,
                  y = Ylf_theor1)) +
  geom_line(linewidth = 1, linetype = 2) +
  annotate(geom = "segment", 
           x    = c(-0.5, -1.5, -2.3, -3.5, 0),
           xend = c(-0.5, -1.5, -2.3, -3.5, 0),
           y    = c(-4.0, -4.0, -4.0, -4.0, 0),
           yend = c(-0.5, -1.5, -2.5, -3.5, -4),
           color = c("black","black","black","black","black"), 
           linewidth = 1, linetype = 2) +
  annotate(geom = "segment", 
           x    = c(-2.3, -0.5, -1.5, -3.5),
           xend = c(-2.3, -2.3, -3.5, -3.5),
           y    = c(-2.3, -0.5, -1.6, -3.6),
           yend = c(-1.9, -2.3, -3.6, -3.2),
           color = c("blue","blue","red","red"), 
           linetype = 1, size = 2, 
           arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", 
           x = c(-1.2, -2.2, -3.56, -0.85, -3.42, -2.22, -0.08),
           y = c(-1.6, -2.6, -3.7, -0.5, -3.1, -1.7, -0.05),
           color = "black", 
           size = 8, 
           label = c('B (pd)', 'C','E', 'A (pd)', 'F (md)', 'D(md)', 'O') ) +
  annotate("text", 
           x = c(-1.4, -2.2, -0.4, -3.4, -0.1),
           y = c(-3.9, -3.9, -3.9, -3.9, -3.9),
           color = "black", 
           size = 8, 
           label = c('B*minute', 'C*minute', 'A*minute', 'E*minute', 'O*minute'), 
           parse = TRUE ) +
  annotate("text", 
           x = c(0.05),
           y = c(0.1),
           color = "white", 
           size = 8, 
           label = c(''), 
           parse = TRUE ) +
  ylim(-4, 0.05) +
  theme_classic() +
  My_Theme +
  labs(title = "Perfecly anisohydric at seasonal scale") +
  xlab('Hydraulic forcing, MPa') +
  ylab('') +
  coord_cartesian(expand = FALSE) # eliminates padding of X and Y axes


grid.arrange(p1, p2, nrow=1)


