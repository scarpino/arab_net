library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(cowplot)

pfish <- function( comm, treat_in_comm, treat_n, total_n){
  x1 = treat_in_comm
  x2 = comm - treat_in_comm
  x3 = treat_n - treat_in_comm
  x4 = total_n - treat_n

  return(fisher.test(matrix(c(x1,x2,x3,x4), ncol=2))$p.value)
}

comms_raw <- read.csv("../data/comm_props_raw_21Jan16.csv")

n_genes <- sum(comms_raw$size)
n_dry <- sum(comms_raw$Drought)
n_cold <- sum(comms_raw$Cold)

comms <- comms_raw %>%
  gather(treatment, treat_in, Drought, Cold)%>%
  mutate(treat_total = ifelse(treatment =="Drought", n_dry, n_cold),
         expect = ifelse(treatment =="Drought", n_dry*prop, n_cold*prop))%>%
  rowwise()%>%
  mutate(p.val = pfish(size, treat_in, treat_total, n_genes))%>%
  ungroup()%>%mutate(p.adj = p.adjust(p.val, method = "BH"))

commsFilter <- comms %>%
  filter(size > 9)%>%
  mutate(p.adj = p.adjust(p.val, method = "BH"),
         relDiff = (treat_in-expect)/treat_total,
         colCat = ifelse(p.adj < 0.05, ifelse(relDiff<0, "down", "up"), "ns"))

