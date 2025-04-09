library(shiny)
library(bslib)
library(ggplot2)
library(shinyBS) # for bsModal

# define user interface
#
source("helpers.R")

ui <- page_sidebar(
  title = "Multiple Comparison Procedures",
  sidebar = sidebar("", position = "left",
    
    helpText(
      "Compare means of several groups with different multiple comparison procedures: 
      choose the number of samples and the sample sizes, the standard deviations, the type of comparisons 
      as well as the procedures to be compared. Afterwards, press 'run simulation'."
    ),
    
    numericInput(
      "nsim",
      "Number of simulations: ",
      value = 1000,
      min = 0,
      max = 10000,
      step = 10
    ),
  
    numericInput(
      "num_groups", 
      "Number of groups: ", 
      value = 3, 
      min = 3, 
      max = 15, 
      step = 1
    ),
  
    uiOutput("sd_inputs"),
  
    sliderInput(
     "samplesize",
     "Set sample size per group",
     min = 0,
     max = 200,
     value = 32
    ),
    
    sliderInput(
      "effectsize",
      "Set effect size",
      min = 0,
      max = 5,
      value = 0.2,
      step = 0.1
    ),
    
    radioButtons(
      "sandwich",
      "Should a heteroscedasticity-consistent estimate of the covariance matrix be used?",
      choices = list(
        "Yes",
        "No"
      )
    ),
  
    selectInput(
     "comparisons",
     label = "Choose the type of comparisons",
     choices = 
       list(
         "many-to-one comparisons", 
         "all pairwise comparisons"), 
    selected = "all pairwise comparisons"
  ),
  
 
  uiOutput("procedure_choices"),
  withMathJax(),

  actionButton(
    "run",
    label = "run simulation"
  )),
  
  mainPanel(
    plotOutput("powerplot"),
    plotOutput("fwerplot"),
    plotOutput("pcerrorplot")
    
  )
)
  
  

server <- function(input, output) {
  
  output$sd_inputs <- renderUI({
    num_groups <- input$num_groups
    
    input_list <- lapply(1:num_groups, function(i) {
      numericInput(paste0("sd_",i),
                   paste("standard deviation for group", i, ":"),
                   value = 1, min = 0, step = 0.1)
    })
    
      
    do.call(tagList, input_list) 
    })
  
  output$procedure_choices <- renderUI({
    req(input$comparisons)
    
    if (input$comparisons == "many-to-one comparisons") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Dunnett", "Step-Down Dunnett")
    } else {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Tukey", "Shaffer", "Westfall")
    }
  
    checkboxes_with_info <- lapply(choices, function(proc) {
      fluidRow(
        column(1, 
               checkboxInput(
                 inputId = paste0("check_", proc),
                 label = NULL,
                 value = FALSE
               )
        ),
        column(7, tags$label(`for` = paste0("check_", proc), proc)),
        column(1, actionLink(
          inputId = paste0("info_", proc),
          label = icon("info-circle")
        ))
      )
    })
    
    tagList(
      tags$label("Choose the procedures to be compared:"),
      tags$div(style = "margin-bottom: 10px;", checkboxes_with_info)
    )
  })
  
  selected_methods <- reactive({
    if (input$comparisons == "many-to-one comparisons") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Dunnett", "Step-Down Dunnett")
    } else {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Tukey", "Shaffer", "Westfall")
    }
    
    Filter(function(proc) {
      isTRUE(input[[paste0("check_", proc)]])
    }, choices)
  })


  
  procedure_info_text <- function(procedure) {
    switch(procedure,
           "Bonferroni" = withMathJax("Single-step procedure that controls the FWER. 
         The unadjusted p-values \\(p_1, p_2, ..., p_m\\) are compared with the threshold \\(\\frac{\\alpha}{m}\\)."),
           "Holm" = withMathJax("Step-down version of Bonferroni that controls the FWER and is at least as powerful as Bonferroni. 
                                The null hypothesis \\(H_{(1)}\\) associated with the smallest p-value \\(p_1\\)
                                is tested first: if \\(p_1 > \\frac{\\alpha}{m}\\), the procedure stops without any rejections.
                                Otherwise, \\(H_{(1)}\\) is rejected and \\(H_{(2)}\\) is tested at the larger significance level \\(\\frac{\\alpha}{(m-1)}\\).
                                These steps are repeated until either the first non-rejection occurs or all hypotheses \\(H_{(1)},...,H_{(m)}\\) are rejected."),
           "Hochberg" = withMathJax("Step-up extension of the Simes test that can be considered a reversed Holm procedure. It is at least as powerful as Holm.
                                    The null hypothesis \\(H_{(m)}\\) associated with the largest p-value \\(p_m\\) is tested first:
                                    if \\(p_m \\leq {\\alpha}\\ \\), all hypotheses \\(H_{(1)},...,H_{(m)}\\) are rejected. 
                                    Otherwise, \\(H_{(m)}\\) is retained and \\(H_{(m-1)}\\) is tested at the smaller significance level \\(\\frac{\\alpha}{2}\\).
                                    If \\(p_{(m-1)} \\leq \\frac{\\alpha}{2}\\), the procedure stops and all hypotheses \\(H_{(1)},...,H_{(m-1)}\\) are rejected.
                                    Theses steps are repeated until either the first rejection occurs or all null hypotheses \\(H_{(1)},...H_{(m)}\\) are retained."),
           "Hommel" = "Improved version of the Simes test that is uniformly more powerful than Hochberg.
           The Simes test is applied to each intersection hypothesis of a closed test procedure.",
           "Dunnett" = "A procedure for comparing mutliple treatment groups against a single control group.
           It accounts for correlations between the comparisons by using the multivariate t-distribution to 
           jointly evaluate all test statistics. The FWER is controlled and the procedure achieves higher power than Bonferroni.",
           "Step-Down Dunnett" = withMathJax(
             "Step-wise extension of Dunnett's test for comparing multiple treatment groups against a single control. 
              The null hypotheses \\( H_{(1)}, \\ldots, H_{(m)} \\) (corresponding to the ordered p-values 
              \\( p_{(1)} \\leq \\ldots \\leq p_{(m)} \\)) are tested sequentially, 
              and the adjusted p-values \\( q_{(i)} = \\max_{I : i \\in I} p_I \\) are obtained, 
              stopping as soon as \\( q_{(i)} > \\alpha \\). 
              (The index set \\( I \\) refers to intersections of elementary hypotheses in the closed testing procedure that include \\( i \\)) 
              The use of max-\\(t\\) tests allows for shortcuts to the full closure. 
              Therefore, power is higher than for Dunnett's test."),
           "Tukey" = "A procedure for all pairwise comparisons of several means. It takes the maximum 
           over the absolute values of all pairwise test statistics that are univariate t distributed.
           Adjusted p-values are computed from the underlying multivarate t distribution thus accounting for correlations between the test statistics.",
           "Shaffer" = "A stepwise extension of the Tukey test using Bonferroni tests.",
           "Westfall" = "A stepwise extension of the Tukey test using truncated closure. Relevant subsets are found by using
           rank conditions involving contrast matrices for the various subsets. The set of subsets considered
           in computing the adjusted p-values is often smaller (than for e.g. the step-down Dunnett test) due to the restricted combination condition.
           The procedure is uniformly more powerful than Tukey and Shaffer."
    )
  }
  
  observe({
    all_choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Dunnett", 
                     "Step-Down Dunnett", "Tukey", "Shaffer", "Westfall")
    
    lapply(all_choices, function(proc) {
      observeEvent(input[[paste0("info_", proc)]], {
        showModal(modalDialog(
          title = paste(proc, "Information"),
          procedure_info_text(proc),
          easyClose = TRUE,
          footer = NULL
        ))
      })
    })
  })
  
  # UI names -> R names 
  
  method_mapping <- c(
    "Bonferroni" = "bonferroni",
    "Holm" = "holm",
    "Hochberg" = "hochberg",
    "Hommel" = "hommel",
    "Dunnett" = "single-step", 
    "Tukey" = "single-step",
    "Shaffer" = "Shaffer",
    "Westfall" = "Westfall",
    "Benjamini-Hochberg" = "BH",
    "Step-Down Dunnett" = "free" 
  )
  
  # generation of data
  simulation <- eventReactive(input$run,{
    std_values <- sapply(1:input$num_groups, function(i) input[[paste0("sd_", i)]])
    
    simulated_data <- simulate_data(
     n_group =  input$num_groups,
     n_samps = input$samplesize,
     effect = input$effectsize,
     std_list = std_values,
     n_sim = input$nsim
      )
    
    
   return(simulated_data)
  })
  
  #########################################################################################
  
  # generation of aov models
  models <- reactive({
    fit_models(simulation())
  })
  
  # generation of multiple comparisons (H0) (many-to-one or all pairwise)
  results_H0 <- reactive({ 
    req(models())
    if (input$sandwich == "Yes") {
      multiple_comp_H0_sandwich(models(), input$comparisons)
    } else {
      multiple_comp_H0(models(), input$comparisons)
    }
  })
  
  # generation of multiple comparisons (H1) (many-to-one or all pairwise)
  results_H1 <- reactive({
    req(models())
    if (input$sandwich == "Yes") {
      multiple_comp_H1_sandwich(models(), input$comparisons)
    } else {
      multiple_comp_H1(models(), input$comparisons)
    }
  })
  
  
  ######################################################################################

  
  procedure_results <- reactive({
    req(selected_methods()) # at least one procedure has to be chosen
    
    results_list_H0 <- list()
    results_list_H1 <- list()
    
    for (method in selected_methods()) {
      r_method <- method_mapping[method]
      
      summaries_H0 <- lapply(results_H0(), function(x) summary(x, test = adjusted(type = r_method)))
      summaries_H1 <- lapply(results_H1(), function(x) summary(x, test = adjusted(type = r_method)))
      
      results_list_H0[[method]] <- summaries_H0
      results_list_H1[[method]] <- summaries_H1
    }
      
    results <- list(
      H0 = results_list_H0,
      H1 = results_list_H1
    )
    
    #print(results)
    return(results)
    
  })
  
  #########################################################################################
  
  # average power
  power_results <- eventReactive(input$run, {
    
    results <- procedure_results()
    results_H1 <- results$H1
    
    power_list <- list()
    
    for (method in selected_methods()) {
      p <- sapply(results_H1[[method]], function(y) {
        mean(y$test$pvalues < 0.05)})
      #print(p)
      power_list[[method]] <- p
      
      
    }
    #print(power_list)
    
    average_power <- lapply(power_list, function(l) mean(unlist(l)))
    #print(average_power)
    
    return(average_power)
    
  })
  
  
  
  # per-comparison error rate (expected proportion of type I errors)
  pc_error_results <- eventReactive(input$run, {
    
    results <- procedure_results()
    results_H0 <- results$H0
    
    error_list <- list()
    
    for (method in selected_methods()) {
      pc <- sapply(results_H0[[method]], function(y) {
        mean(y$test$pvalues < 0.05)})
      error_list[[method]] <- pc
      
    }
    
    #print(error_list)
    pc_error <- lapply(error_list, function(vec) mean(unlist(vec)))
    #print(pc_error)
    
    return(pc_error)
    
  })
  
  
  
  # family-wise error rate
  fwer_results <- eventReactive(input$run, {
    
    results <- procedure_results()
    results_H0 <- results$H0
    
    
    fwer_list <- list()
    
    for (method in selected_methods()) {
      method_results <- results_H0[[method]]
      reject_any <- sapply(method_results, function(r) {
        any(r$test$pvalues < 0.05)})
      fwer_list[[method]] <- mean(reject_any)  
    
    }
    
    return(fwer_list)
  })
  
  #########################################################################################
    
  # show average power results
    
    output$powerplot <- renderPlot({
      power <- power_results()
      req(power)
      
      power_df <- data.frame(
        Procedure = names(power),
       Average_Power = unlist(power)
      )
      
      #print(power_df)
      
      ggplot(power_df, aes(x = Procedure, y = Average_Power)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Multiple Comparisons Procedures: Average Power",
             x = "Multiple Comparison Procedures",
             y = "Average Power") +
        theme_minimal() 
    })
  
  # show per-comparison error results
  
   output$pcerrorplot <- renderPlot({
     pc_error <- pc_error_results()
     req(pc_error)
     
     pc_error_df <- data.frame(
       Procedure = names(pc_error),
       PC_Error = unlist(pc_error)
     )
     
     ggplot(pc_error_df, aes(x = Procedure, y = PC_Error)) +
       geom_bar(stat = "identity", position = "dodge") +
       labs(title = "Multiple Comparisons Procedures: Per-comparison error rate",
            x = "Multiple Comparison Procedures",
            y = "Per-comparison error rate") +
       theme_minimal()
   })
   
   # show family wise error results
   
   output$fwerplot <- renderPlot({
     fwer <- fwer_results()
     req(fwer)
     
     fwer_df <- data.frame(
       Procedure = names(fwer),
       FWER = unlist(fwer)
     )
     
     ggplot(fwer_df, aes(x = Procedure, y = FWER)) +
       geom_bar(stat = "identity", position = "dodge") +
       labs(title = "Multiple Comparisons Procedures: Family-wise error rate",
            x = "Multiple Comparison Procedures",
            y = "Family-wise error rate") +
       theme_minimal()
   })
    
}


# Run the app ----
shinyApp(ui = ui, server = server)
    