##-----------------
##  POST  
##-----------------
## Created: 03 may 2024
##-----------------
### Commentary:
##'
##' Function to extract posterior means
##' of quantities of interest for BP-SCM
##'
##-----------------

##' Proxys labels
lib_prox <- . %>%
  mutate(prox = factor(prox,
                       levels=1:6,
                       labels = c("Y1","Y2","X1","X2","X3","sh"))
  )

##' Summary of parameters
##' 
##' @title get_pars
##' @param .x : the BP-SCM stanfit
##' @param pars : params to summarise
##' @return Tibble with posterior summary for model parameters
##' 
##' @author 
get_pars <- function(.x, pars = NULL){
  if (is.null(pars)) {
    pars<-intersect(.x@model_pars,c("Int","sd_c","sd_b","tau","tau_sh","tau_sp","w","phi","p_sh","delta"))
  } 
  
  samp <- rstan::extract(.x, pars = pars)
  res <- tibble(param = names(samp)) %>%
    mutate(su = map(param,~ {
      as.matrix(samp[[.x]]) %>%
        tidybayes::summarise_draws("mean", "median","sd", "rhat","ess_bulk","ess_tail",
                                   ~quantile(.x, probs = c(0.025, 0.975)))
    }
    )) %>%
    unnest(su) %>%
    mutate(prox = str_extract(variable,"\\d+") %>% as.numeric) %>%
    select(param,  prox, mean, median, sd, low = `2.5%`, up = `97.5%`, rhat, contains("ess"))
  
  res%<>%
    mutate(prox = ifelse(str_detect(param,"sd_|Int|phi"),prox+2,prox))%>%
    bind_rows(res %>% filter(param=="sd_c", prox==3) %>%
                mutate(param="tau", prox=5))
  
  if (.x@model_name == "pscm")
    res %<>%
      mutate(prox = ifelse(param=="sd_c",6,prox))  
  
  res %>% 
    lib_prox %>%
    arrange(param,prox) %>%
    rename(pm=mean)
}

##' Age rates rho_ik and r_ik
##'
##' @title get_rates
##' @param .x : the BP-SCM stanfit
##' @param .y : data list used for the fit (obtained with prep_data)
##' @return Tibble with posterior summary for age rates
##' @author 
get_rates <- function(.x,.y){
  res <- spread_draws(.x,r[id_data,prox]) %>%
    summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                    ~quantile(.x, probs = c(0.025, 0.975))) %>%
    ungroup()
  res%<>% 
    bind_rows(spread_draws(.x,rho[id_data,prox]) %>%
                summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                                ~quantile(.x, probs = c(0.025, 0.975))) %>%
                mutate(prox = prox+2))%>%
    ungroup()
  res  %>%
    lib_prox %>%
    left_join(.y$id_dt %>% 
                select(age,id_data) %>%
                unique(),by="id_data") %>%
    select(prox,age,pm = mean,sd,low=`2.5%`,up=`97.5%`,rhat,ess_bulk,ess_tail)
}

##' Total number of cases mu_ij. and lambda_ij.
##'
##' @title get_lm
##' @param .x : the BP-SCM stanfit
##' @param .y : data list used for the fit (obtained with prep_data)
##' @return Tibble with posterior summary for mean total number of cases per geographical area and proxy
##' @author 
get_lm <- function(.x,.y){
  res <- spread_draws(.x,mu[id_geo,prox]) %>%
    summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                    ~quantile(.x, probs = c(0.025, 0.975)))
  res %<>% bind_rows(
    spread_draws(.x,lbd[id_geo,prox]) %>%
      summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                      ~quantile(.x, probs = c(0.025, 0.975)))%>%
      mutate(prox=prox+2)
  )
  
  res %>%
    lib_prox%>%
    ungroup() %>%
    left_join(.y$id_dt %>% select(id_geo, geo) %>%  unique(), by="id_geo") %>%
    ungroup() %>%
    select(prox,geo,pm = mean,sd,low=`2.5%`,up=`97.5%`,rhat,ess_bulk,ess_tail)
}

##' Relative risks nu_ij and theta_ij
##'
##' @title get_rr
##' @param .x : the BP-SCM stanfit
##' @param .y : data list used for the fit (obtained with prep_data)
##' @return Tibble with posterior summary for relative risks  per geographical area and proxy
##' @author 
get_rr <- function(.x,.y){
  res <- spread_draws(.x,theta[id_geo,prox]) %>%
    summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                    ~quantile(.x, probs = c(0.025, 0.975)))
  
  res %<>%  bind_rows(
    spread_draws(.x,nu[id_geo,prox]) %>%
      summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                      ~quantile(.x, probs = c(0.025, 0.975))) %>%
      mutate(prox = prox+2)
  )
  res %<>%
    lib_prox%>%
    ungroup() %>%
    left_join(.y$id_dt %>% select(id_geo, geo) %>%  unique(), by="id_geo") %>%
    ungroup() %>%
    select(prox,geo,pm = mean,sd,low=`2.5%`,up=`97.5%`,rhat,ess_bulk,ess_tail) 
}

##' Random effects b_ij and c_j
##'
##' @title get_re
##' @param .x : the BP-SCM stanfit
##' @param .y : data list used for the fit (obtained with prep_data)
##' @return Tibble with posterior summary for rnadom effects per geographical area
##' @details b_ij random effects are scaled by their standard devriation sigma_i
##' @author 
get_re <- function(.x,.y){
  res <- spread_draws(.x,bc[id_geo,prox]) %>%
    summarise_draws("mean", "sd", "rhat","ess_bulk","ess_tail",
                    ~quantile(.x, probs = c(0.025, 0.975))) 
  cat("End RE\n")
  res %>%
    lib_prox%>%
    ungroup() %>%
    left_join(.y$id_dt %>% select(id_geo, geo) %>%  unique(), by="id_geo") %>%
    ungroup() %>%
    select(prox,geo,pm = mean,sd,low=`2.5%`,up=`97.5%`,rhat,ess_bulk,ess_tail) 
}


