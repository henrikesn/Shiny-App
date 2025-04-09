# helpers.R for Simulation App

library(multcomp)
library(dplyr)
library(sandwich)

simulate_data <- function(n_group, n_samps, effect, std_list, n_sim){
  set.seed(123)
  
  mean_H0 <- 3
  means_H1 <- mean_H0 + effect * seq(0, n_group-1) #!!!!!
  #print(means_H1)
  
  simulation_H0 <- as.vector(replicate(n_sim,
                                       {rnorm(n_samps*n_group, 
                                              mean = mean_H0, 
                                              sd = rep(std_list, each = n_samps))}))

  

  simulation_H1 <- as.vector(replicate(n_sim, {
    unlist(lapply(1:n_group, function(g) {
      rnorm(n_samps, mean = means_H1[g], sd = std_list[g])
    }))
  }))
  
  
  group_labels <- rep(1:n_group, each = n_samps)
  
  group_labels_H0 <- rep(group_labels, times = n_sim)
  group_labels_H1 <- rep(group_labels, times = n_sim)
  

  
  data_H0 <- data.frame(
    values = simulation_H0,
    group = group_labels_H0,
    hypothesis = "H0",
    sim_id = rep(1:n_sim, each = n_group * n_samps)
  )
  
  
  data_H1 <- data.frame(
    values = simulation_H1,
    group = group_labels_H1,
    hypothesis = "H1",
    sim_id = rep(1:n_sim, each = n_group * n_samps)
  )
  
  
  data_H0$group <- as.factor(data_H0$group)
  data_H1$group <- as.factor(data_H1$group)
  
  data_sim <- rbind(data_H0, data_H1)
   
  
  return(data_sim)
  
}


fit_models <- function(data) {
  
  data_H0 <- data %>% filter(hypothesis == 'H0')
  data_H1 <- data %>% filter(hypothesis == 'H1') 
  
  sim_H0 <- split(data_H0, data_H0$sim_id)
  sim_H1 <- split(data_H1, data_H1$sim_id)
  
  models_H0 <- lapply(sim_H0, function(d) aov(values ~ group, data = d))
  models_H1 <- lapply(sim_H1, function(d) aov(values ~ group, data = d))
  
  return(list(H0 = models_H0, H1 = models_H1))
}


multiple_comp_H0 <- function(models, comparison) {
  
  lapply(models$H0, function(model) {
    switch(comparison,
           "many-to-one comparisons" = glht(model, linfct = mcp(group = "Dunnett"), alternative = 'two.sided'),
           "all pairwise comparisons" = glht(model, linfct = mcp(group = "Tukey"), alternative = 'two.sided'))
  })
}


multiple_comp_H1 <- function(models, comparison) {
  
  lapply(models$H1, function(model) {
    switch(comparison,
         "many-to-one comparisons" = glht(model, linfct = mcp(group = "Dunnett"), alternative = 'two.sided'),
         "all pairwise comparisons" = glht(model, linfct = mcp(group = "Tukey"), alternative = 'two.sided'))
  })
}
  
multiple_comp_H0_sandwich <- function(models, comparison) {
  
  lapply(models$H0, function(model) {
    switch(comparison,
           "many-to-one comparisons" = glht(model, linfct = mcp(group = "Dunnett"), vcov = sandwich, alternative = 'two.sided'),
           "all pairwise comparisons" = glht(model, linfct = mcp(group = "Tukey"), vcov = sandwich, alternative = 'two.sided'))
  })
}


multiple_comp_H1_sandwich <- function(models, comparison) {
  
  lapply(models$H1, function(model) {
    switch(comparison,
           "many-to-one comparisons" = glht(model, linfct = mcp(group = "Dunnett"), vcov = sandwich, alternative = 'two.sided'),
           "all pairwise comparisons" = glht(model, linfct = mcp(group = "Tukey"), vcov = sandwich, alternative = 'two.sided'))
  })
}




# n_group <- 5
# n_samps <- 5
# effect <- 0.2
# std_list <- c(1,2,3,4,5)
# test <- simulate_data(n_group, n_samps, effect, std_list)
# 
# models <- fit_models(test)
# modelH0 <- models$H0
# model1H0 <- modelH0[[1]]
# 
# #comp <- multiple_comp_H0(model1H0, "many-to-one comparisons")
# 
# versuch <- glht(model1H0, linfct = mcp(group = "Dunnett"))
# p <- summary(versuch, test = adjusted(NULL))
# p1 <- summary(versuch, test = adjusted("single-step"))
