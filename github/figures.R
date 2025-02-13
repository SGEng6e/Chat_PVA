## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, cache = TRUE, include = FALSE, fig.width = 12, fig.height = 8,
                      fig.align='center', out.width="90%", purl = knitr::hook_purl)
library(tidyverse)
library(ggthemes)
library(MCMCvis)
library(coda)
library(leaflet)
library(figpatch)
library(patchwork)
library(prettyunits)

### LOAD DATA ##############################################################################################

load("./pva_chain1.rda")
load("./pva_chain2.rda")
load("./pva_chain3.rda")
load("./pva_chain4.rda")

pva_chains <- mcmc.list(mget(ls(pattern = "pva_chain"))) ; rm(list = ls(pattern = "pva_chain\\d"))
load("./ybch_con.rda")
load("./ybch_dat.rda")

### PROCESS DATA ###########################################################################################

pva_results <- MCMCsummary(pva_chains, probs = c(0.025,0.05,0.5,0.95,0.975)) %>% 
  rownames_to_column("Term") %>% 
  mutate(Parameter = str_extract(Term, ".*(?=(\\[.{1,5}\\])|((?<!\\])$))"),
         Idx = str_extract(Term, "(?<!((tem\\[\\d, )|(det\\[\\d, )|(_d\\[)|(_t\\[)))\\d{1,2}(?=\\]$)"),
         Scenario = str_extract(Term, "(?<!((tem\\[)|(det\\[)))(?<=\\[)\\d(?=,)"),
         Year = if_else(Parameter %in% c("Pr_ext","Pr_qext","Pr_K"),
                        as.character(as.numeric(Idx)+2021),
                        as.character(as.numeric(Idx)+2001)))

pva_smry <- pva_results %>% 
  select(Term, mean, sd, `5%`, `95%`, Parameter, Idx, Scenario, Year) 

N <- pva_smry %>% 
  filter(str_detect(Parameter, "^N")) %>% 
  mutate(across(c(`5%`, `95%`), \(x) if_else(x > 1000, 1000, x)),
         mean = if_else(mean > 1100, NA_real_, mean),
         across(where(is.numeric), \(x) pretty_signif(x, digits = 2))) %>% 
  split(.$Parameter) %>% 
  map(\(x) split(x, x$Scenario))

phi_a <- pva_smry %>% 
  filter(str_detect(Parameter, "phi_a")) %>% 
  mutate(across(where(is.numeric), \(x) pretty_signif(x, digits = 2)),
         context = if_else(Year > 2020, "proj", "hist")) %>% 
  split(.$context)

phi_j <- pva_smry %>% 
  filter(str_detect(Parameter, "phi_j")) %>% 
  mutate(across(where(is.numeric), \(x) pretty_signif(x, digits = 2)),
         context = if_else(Year > 2020, "proj", "hist")) %>% 
  split(.$context)

pi <- pva_smry %>% 
  filter(str_detect(Parameter, "^pi")) %>% 
  mutate(across(where(is.numeric), \(x) pretty_signif(x, digits = 2)),
         context = if_else(Year > 2022, "proj", "hist")) %>% 
  split(.$context)

ela_sen <- pva_smry %>% 
  select(Term, mean, sd, `5%`, `95%`) %>% 
  filter(str_detect(Term, "(ela)|(sen)")) %>% 
  mutate(quantity = if_else(str_detect(Term, "ela"), "Elasticity", "Sensitivity"),
         Term = factor(case_match(Term, c("sen_phij","ela_phij") ~ 2,c("sen_phia","ela_phia") ~ 1,
                                  c("sen_fer","ela_fer") ~ 3,c("sen_omega","ela_omega") ~ 5,
                                  c("sen_alpha","ela_alpha") ~ 4),
                       levels = c(5:1), labels = c(
                         "Age at last\nbreeding","Age at first\nbreeding",
                         "Fertility","Juvenile\nsurvival","Adult\nsurvival")))

thm_obj <- theme(
  legend.text = element_text(size = 13),
  legend.key = element_rect(colour = "transparent",fill = "transparent"),
  legend.title = element_text(size = 11.5),
  panel.spacing = unit(5, "mm"),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 13),
  axis.title = element_text(size = 14),
  strip.background = element_blank(),
  strip.text = element_text(size = 12),
  panel.background = element_blank(),
  panel.grid = element_blank())

# ## ----fig1, fig.cap="Representative satellite image of the Okanagan river where channelization,
# agricultural expansion, and urban development during the 1900s led to declines of riparian-
# associated species in the south Okanagan Valley of Canada. Southern Mountain chat (*Icteria virens
# auricollis*) territories were monitored along the river between 2002 and 2021. Cluster points
# indicate approximate locations of clusters of nests, while blue pins indicate the approximate
# locations of single nests.",
# include = T, out.width = '100%', fig.dpi=240

# Data are not shown to protect the privacy of private land owners.

## ----fig2, fig.cap="Abundance of breeding pairs of Southern Mountain chat (*Icteria virens auricollis*)
# in the south Okanagan Valley of Canada between 2002 and 2021, and projected under 5 scenarios to 2041.
# Under scenario 1, demographic rates of fertility, adult survival, and juvenile survival were projected
# according to their 2002-2021 means. In scenarios 2 and 3, rates of parasitism by the brown-headed
# cowbird (*Molothrus ater*) increased and decreased by 20% respectively. Under scenarios 4 and 5, we
# simulated a habitat degradation event that suppressed the breeding population by 10% and 20% annually,
# respectively. Solid lines and shaded ribbons indicate posterior means and 90% credible intervals of
# population estimates. Black horizontal dashed line at *y=1000* indicates the approximate carrying
# capacity of our study region.",
# include = T, fig.dim=c(7.5,6), fig.dpi=240, cache=FALSE, warning=FALSE----

N$N_tot %>% 
  list_rbind() %>% 
  mutate(Year = as.numeric(Year)) %>% 
  filter(!(Scenario > 1 & Year < 2022)) %>%
  ggplot(aes(x = Year, y = as.numeric(mean), colour = Scenario, fill = Scenario, group = Scenario))+
  geom_line(lwd = 0.8)+
  geom_ribbon(aes(x=Year, ymax=as.numeric(`95%`), ymin=as.numeric(`5%`),fill=Scenario),colour=NA, alpha=0.2)+
  scale_fill_colorblind()+
  scale_colour_colorblind()+
  geom_hline(yintercept = 1000, linetype = 2, linewidth = 0.75, colour = "black")+
  geom_vline(xintercept = 2022, linetype = 2, linewidth = 0.75, colour = "firebrick2")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1000,150), limits = c(0,1050))+
  labs(y = "Breeding pair abundance")+
  thm_obj+
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.8),
        axis.title.y = element_text(vjust = 2.5))


## ----fig3, fig.cap="The adult survival rate of the Southern Mountain chat (*Icteria virens auricollis*)
# has the largest proportional effect (elasticity) on the population asymptotic growth rates, followed by
# juvenile survival, and fertility. Asymptotic growth rate is least sensitive to structural population 
# parameters of age at first and last breeding. Elasticity is computed as sensitivity rescaled to account
# for the magnitude of both the population growth rate and each respective matrix model element. Point 
# estimates and line ranges are posterior means and 90% Bayesian credible intervals.", 
# include = T, fig.dim=c(6.5,5), fig.dpi=240, cache=FALSE----

ggplot(ela_sen, aes(x = mean, y = Term))+
  facet_wrap(~quantity, scales = "free_x")+
  geom_point()+
  geom_linerange(aes(x = mean, xmax = `95%`, xmin = `5%`, y = Term))+
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "darkgrey")+
  labs(y = NULL, x = "Quantity estimate (dimensionless)")+
  theme(panel.spacing.x = unit(25, "pt"))+
  thm_obj
