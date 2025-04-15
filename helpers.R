# helpers.R for Simulation App

library(multcomp)
library(dplyr)
library(sandwich)
library(effectsize)

simulate_data <- function(n_group, n_samps_list, effect, std_list, n_sim){
  set.seed(123)
  
  mean_H0 <- 3
  means_H1 <- mean_H0 + effect * seq(0, n_group-1) #!!!!!
  #print(means_H1)

  simulation_H0 <- unlist(replicate(n_sim, {
    unlist(lapply(1:n_group, function(g) {
      
      rnorm(n = n_samps_list[g], mean = mean_H0, sd = std_list[g])
      
    }))
  }, simplify = FALSE))
  
  simulation_H1 <- unlist(replicate(n_sim, {
    unlist(lapply(1:n_group, function(g) {
      
      rnorm(n = n_samps_list[g], mean = means_H1[g], sd = std_list[g])
      
    }))
  }, simplify = FALSE))
  
  
  group_labels <- unlist(lapply(1:n_group, function(g) rep(g, each = n_samps_list[g])))
  
  group_labels_H0 <- rep(group_labels, times = n_sim)
  group_labels_H1 <- rep(group_labels, times = n_sim)
  
  sim_ids <- rep(1:n_sim, each = sum(n_samps_list))
  
  data_H0 <- data.frame(
    values = simulation_H0,
    group = group_labels_H0,
    hypothesis = "H0",
    sim_id = sim_ids
  )
  
  
  data_H1 <- data.frame(
    values = simulation_H1,
    group = group_labels_H1,
    hypothesis = "H1",
    sim_id = sim_ids
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


effectsizes_mto <- function(data, ref, effect, n_samps_list, std_list) {
  group_ids <- unique(data$group)
  combinations <- setdiff(group_ids, ref)
  
  effectlist <- lapply(unique(data$sim_id), function(sim) {
    
    data_sub <- data%>%filter(hypothesis == 'H1', sim_id == sim)
    
    lapply(combinations, function(comb) {
      
      pair <- unique(c(as.character(comb), as.character(ref)))
      
      data_sub_pair <- data_sub %>%
        filter(group %in% pair) %>%
        mutate(group = factor(as.character(group), levels = as.character(pair)))

      
      cohen <- cohens_d(values ~ group, ref.group = as.character(ref), data = data_sub_pair)
      
      means_ref <- data_sub_pair %>%
        filter(group == as.character(ref)) %>%
        summarise(mean_ref = mean(values))
      
      means_comp <- data_sub_pair %>%
        filter(as.character(group) == as.character(comb)) %>%
        summarise(mean_comp = mean(values))
      
      mean_diff <- means_comp$mean_comp - means_ref$mean_ref
      
      true_mean_diff <- (as.integer(pair[1]) - as.integer(pair[2])) * effect
      
      n_ref <- n_samps_list[as.integer(pair[2])]
      
      n_comp <- n_samps_list[as.integer(pair[1])]
      
      std_ref <- std_list[as.numeric(pair[2])]
      
      std_comp <- std_list[as.numeric(pair[1])]
      
      pooled_std <- sqrt((std_ref^2 * (n_ref-1) + std_comp^2 * (n_comp-1))/(n_ref+n_comp-2))
      
      true_cohen <- abs(true_mean_diff) / pooled_std 
      
      
      
      effects <- data.frame(
        sim = sim,
        ref = pair[2],
        comp = pair[1],
        mean_ref = means_ref$mean_ref,
        mean_comp = means_comp$mean_comp,
        mean_diff = mean_diff,
        cohens_d = abs(cohen$Cohens_d),
        true_mean_diff = true_mean_diff,
        true_cohen = true_cohen,
        comparison = paste0(pair[1], " vs ", pair[2])
        )
    })
   }) 

  
  temp <- unlist(effectlist, recursive = FALSE)
  result <- bind_rows(temp)
  
  return(result)
}

effectsizes_pw <- function(data, effect, n_samps_list, std_list) {
  
  group_ids <- unique(data$group)
  combinations <- combn(group_ids, 2, simplify = FALSE)
  
  effectlist <- lapply(unique(data$sim_id), function(sim) {  # list of list of data.frames
    
    data_sub <- data %>% filter(hypothesis == 'H1', sim_id == sim)
    
    lapply(combinations, function(pair) {
      
      data_sub_pair <- data_sub %>%
        filter(group %in% pair) %>%
        mutate(group = factor(as.character(group), levels = as.character(pair)))
      
      cohen <- cohens_d(values ~ group, data = data_sub_pair)
      
      means_ref <- data_sub_pair %>%
        filter(as.character(group) == as.character(pair[1])) %>%
        summarise(mean_ref = mean(values))
      
      means_comp <- data_sub_pair %>%
        filter(as.character(group) == as.character(pair[2])) %>%
        summarise(mean_comp = mean(values))
      
      mean_diff = means_comp$mean_comp - means_ref$mean_ref
      
      true_mean_diff = (as.numeric(pair[2]) - as.numeric(pair[1])) * effect
      
      n_ref <- n_samps_list[pair[1]]
      
      n_comp <- n_samps_list[pair[2]]
      
      std_ref <- std_list[pair[1]]
      
      std_comp <- std_list[pair[2]]
      
      pooled_std <- sqrt((std_ref^2 * (n_ref-1) + std_comp^2 * (n_comp-1))/(n_ref+n_comp-2))
      
      true_cohen <- abs(true_mean_diff) / pooled_std 
      

      effects <- data.frame(
        sim = sim,
        ref = pair[1],
        comp = pair[2],
        mean_ref = means_ref$mean_ref,
        mean_comp = means_comp$mean_comp,
        mean_diff = mean_diff,
        cohens_d = abs(cohen$Cohens_d),
        true_mean_diff = true_mean_diff,
        true_cohen = true_cohen,
        comparison = paste0(pair[2], " vs ", pair[1])
      )
 
    })
  })
  
  temp <- unlist(effectlist, recursive = FALSE)
  result <- bind_rows(temp)
  
  return(result)
}
  



##########################################################################


n_group <- 5
n_samps_list <- c(2,3,4,5,6)
effect <- 0.2
std_list <- c(1,2,3,4,5)
n_sim <- 100
ref <- 3
test <- simulate_data(n_group, n_samps_list, effect, std_list, n_sim)

data_H1 <- test%>%filter(hypothesis == 'H1')


versuch <- effectsizes_pw(test, 0.5, n_samps_list, std_list)
versuchmto <- effectsizes_mto(test,2, 0.5, n_samps_list, std_list)

gefiltert <- versuch %>% distinct(comparison, true_cohen)
gfiltertmto <- versuchmto %>% distinct(comparison, true_cohen)

