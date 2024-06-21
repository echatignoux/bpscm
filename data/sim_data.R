##-----------------
##  SIM_DATA  
##-----------------
## Author: Edouard Chatignoux <edouard.chatignoux@gmail.com>
## Created: 03 may 2024
##-----------------
### Commentary:
##'
##'  
##' 
##-----------------

###' Init
##'==================================== 
####' Packages
##'---------------------------------------- 

library(magrittr)
library(tidyverse)
library(INLA)
library(sf)
theme_set(theme_bw())
bis <- geom_abline(intercept=0, slope=1)
source("./R/prep_data.R")


seed <- 24
set.seed(seed)

###'Read shape and age rate ==================================== 
shap_geo <- readRDS("./data/shap_geo.rds")
lop <- readRDS("./data/lop.rds")

###' Simulate data ==================================== 
####'Generate icar random effect  ---------------------------------------- 
## Scaled adjency matrix 
n <- nrow(shap_geo)
cp <- scale_shap(shap_geo)
spdep::nb2INLA(spdep::poly2nb(shap_geo), file="graph.inla")

A_cstr <- map(cp$scales$comp,~{
  w<-which(cp$dt_comp$comp==.x)
  (cp$adj_mat[w,]  %>% as.matrix %>% colSums()>0) %>%
    as.numeric() %>%
    matrix(ncol=n)
}) %>% 
  do.call(rbind, .)
A_cstr <- A_cstr[rowSums(A_cstr)!=0,] 
if (class(A_cstr)[1]=="numeric") 
  A_cstr%<>% matrix(nrow=1)
e_cstr <- rep(0,nrow(A_cstr))

Q = INLA:::inla.pc.bym.Q("graph.inla")
Q = INLA:::inla.scale.model(Q,  constr=list(A=A_cstr, e=e_cstr))
QQ = Q
diag(QQ) = diag(QQ) + 1e-6

## icar and iid random effects for b_ij and c_j 
dt_re <- shap_geo %>%
  as_tibble() %>%
  select(geo) %>% 
  mutate(c_icar = inla.qsample(1, Q=QQ, constr = list(A=A_cstr, e=e_cstr), seed = seed) %>% as.numeric,
         b1_icar = inla.qsample(1, Q=QQ, constr = list(A=A_cstr, e=e_cstr), seed = seed + 1) %>% as.numeric,
         b2_icar = inla.qsample(1, Q=QQ, constr = list(A=A_cstr, e=e_cstr), seed = seed + 2) %>% as.numeric)  %>%
  mutate(c_iid = rnorm(n),
         b1_iid = rnorm(n),
         b2_iid = rnorm(n)) 

####'Parameters for the lambda ---------------------------------------- 
param <- tibble(prox = c("X1","X2","X3"),
                N = c(5,7,10),
                sd_c = c(0.15,0.1,0.3),
                sd_b = c(0.25,0.15,0),
                phi=c(.5,.7,.8))

####'Generate BP distributed data according to the parameters ---------------------------------------- 
## Scale rates to match N levels

dt <- lop %>%
  left_join(param, by  = c("prox")) %>% 
  group_by(prox) %>%
  mutate(ar = ar * N / sum(ar*py) * length(unique(geo))
         )
  
## Set c_j, b_ij and simulate the X_ijk
dt %<>%
  left_join(dt_re, by = "geo") %>%
  mutate( c = sqrt(phi) * c_icar + sqrt(1-phi) * c_iid,
         b1 = sqrt(phi) * b1_icar + sqrt(1-phi) * b1_iid,
         b2 = sqrt(phi) * b2_icar + sqrt(1-phi) * b2_iid
         ) %>%
  mutate(
    re_tm = case_when(prox=="X1"~b1,
                      prox=="X2"~b2,
                      prox=="X3"~c),
    re =sd_c*c+ifelse(prox=="X1",sd_b*b1,sd_b*b2),
    lbd=py*ar*exp(re)) %>% 
  rowwise() %>% 
  mutate(Y=rpois(n=1, lambda = lbd)) %>%
  ungroup()

## Sum the X_ijk to get Y_ijk
dt %<>% 
  group_by(age,geo, py) %>%
  mutate_at( vars(Y,lbd,ar), ~ .x + sum(.x * (prox=="X3"))) %>%
  filter(prox != "X3") %>%
  mutate(prox = ifelse(prox=="X1","Y1","Y2")) %>%
  bind_rows(dt) %>%
  select(geo,age,prox,py,Y,lbd,ar,re = re_tm) 

saveRDS(dt,"./data/data_sim.rds")

## sim_data.R ends here
