# BETA REGRESSION MASTER FUNCTION
BetaRegression_RareVariants <- function(dataset, adjustment){ 
  columns_ignore <- c("percentage", "marker", "BRIDGES_ID", "ER_statusIndex", "AgeDiagIndex", "study")
  
  if(adjustment == "baseline"){ # Age and Study
    
    matrix_to_fill <- matrix(nrow = length(which(!names(dataset) %in% columns_ignore)), ncol = 15)
    
    dataset <- dataset %>%
      mutate(ER_statusIndex = case_when(ER_statusIndex == 888 ~ NA_real_, TRUE ~ ER_statusIndex)) 
    
    colnames(matrix_to_fill) <- c("gene", "estimate", "ci.lb", "ci.ub", "std.error", "p.value", "vif", "pseudo.r.squared", "loglik",
                                  "estimate.age", "ci.lb.age", "ci.ub.age", "std.error.age", "p.value.age", "vif.age")
    possible_error_beta_regression <- function(index){
      out <- tryCatch(exp = {
        suppressWarnings(
          Pr <- betareg(percentage ~ dataset[[index]] + AgeDiagIndex + study, data = dataset, warning = FALSE, link = 'logit')
        )
        return(Pr)
      },
      error = function(e){ #message("error")
        return(c(colnames(dataset[i]), rep(888, 14)))
      }, 
      warning = function(w){ #message("warning")
        suppressWarnings(
          return(c(colnames(dataset[i]), rep(888, 14)))
        )}
      )
      return(out)
    }
    
    for(i in which(!names(dataset) %in% columns_ignore)){
      output <- possible_error_beta_regression(i)
      
      
      if (class(output) == "character") { # When output = c(colnames(dosage_icogs[x]), rep(NA, 8)) meaning an error occured, then add it to the matrix right away.
        matrix_to_fill[i-length(columns_ignore),] <- output
        
      } else if (class(output) == "betareg") {# if a true answer, then put it in the confint line
        PR_ci <- confint(output, level = 0.95)
        Pr_sum <- summary(output)
        Pr_vif <- vif(output)
        intermediate_SNP_data <- c(colnames(dataset[i]), # SNP name
                                   Pr_sum$coefficients$mean[2,1], # Estimate gene
                                   PR_ci[2,1], # ci.lb gene
                                   PR_ci[2,2], # ci.ub gene
                                   Pr_sum$coefficients$mean[2,2], # std.error gene
                                   Pr_sum$coefficients$mean[2,4], # p.value gene
                                   Pr_vif[1], # vif gene
                                   
                                   Pr_sum$pseudo.r.squared, # pseudo.r.squared
                                   Pr_sum$loglik, # loglik
                                   
                                   Pr_sum$coefficients$mean[3,1], # estimate.age
                                   PR_ci[3,1], # ci.lb.age
                                   PR_ci[3,2], # ci.ub.age
                                   Pr_sum$coefficients$mean[3,2], # std.error.age
                                   Pr_sum$coefficients$mean[3,4], # p.value.age
                                   Pr_vif[2] # vif.age
        )
        matrix_to_fill[i-length(columns_ignore),] <- intermediate_SNP_data
      }
    }
    
  } 
  else if(adjustment == "adjusted"){ # ER, Age and Study as covariates
    
    matrix_to_fill <- matrix(nrow = length(which(!names(dataset) %in% columns_ignore)), ncol = 21)
    
    dataset <- dataset %>%
      mutate(ER_statusIndex = case_when(ER_statusIndex == 888 ~ NA_real_, TRUE ~ ER_statusIndex)) 
    
    
    colnames(matrix_to_fill) <- c("gene", "estimate", "ci.lb", "ci.ub", "std.error", "p.value", "vif", "pseudo.r.squared", "loglik",
                                  "estimate.age", "ci.lb.age", "ci.ub.age", "std.error.age", "p.value.age", "vif.age",
                                  "estimate.ER", "ci.lb.ER", "ci.ub.ER", "std.error.ER", "p.value.ER", "vif.ER")
    possible_error_beta_regression <- function(index){
      out <- tryCatch(exp = {
        suppressWarnings(
          Pr <- betareg(percentage ~ dataset[[index]] + ER_statusIndex + AgeDiagIndex + study, data = dataset, warning = FALSE, link = 'logit')
        )
        return(Pr)
      },
      error = function(e){
        # message("error")
        return(c(colnames(dataset[i]), rep(888, 20)))
      }, 
      warning = function(w){
        # message("warning")
        suppressWarnings(
          return(c(colnames(dataset[i]), rep(888, 20)))
        )}
      )
      return(out)
    }
    
    
    for(i in which(!names(dataset) %in% columns_ignore)){
      output <- possible_error_beta_regression(i)
      
      if (class(output) == "character") { # When output = c(colnames(dosage_icogs[x]), rep(NA, 8)) meaning an error occured, then add it to the matrix right away.
        matrix_to_fill[i-length(columns_ignore),] <- output
        
      } else if (class(output) == "betareg") {# if a true answer, then put it in the confint line
        PR_ci <- confint(output, level = 0.95)
        Pr_sum <- summary(output)
        Pr_vif <- vif(output)
        intermediate_SNP_data <- c(colnames(dataset[i]), # SNP name
                                   Pr_sum$coefficients$mean[2,1], # Estimate gene
                                   PR_ci[2,1], # ci.lb gene
                                   PR_ci[2,2], # ci.ub gene
                                   Pr_sum$coefficients$mean[2,2], # std.error gene
                                   Pr_sum$coefficients$mean[2,4], # p.value gene
                                   Pr_vif[1], # vif gene 
                                   
                                   Pr_sum$pseudo.r.squared, # pseudo.r.squared
                                   Pr_sum$loglik, # loglik
                                   
                                   Pr_sum$coefficients$mean[4,1], # estimate.age
                                   PR_ci[4,1], # ci.lb.age
                                   PR_ci[4,2], # ci.ub.age
                                   Pr_sum$coefficients$mean[4,2], # std.error.age
                                   Pr_sum$coefficients$mean[4,4], # p.value.age
                                   Pr_vif[3], # vif.age
                                   
                                   Pr_sum$coefficients$mean[3,1], # estimate.ER
                                   PR_ci[3,1], # ci.lb.ER
                                   PR_ci[3,2], # ci.ub.ER
                                   Pr_sum$coefficients$mean[3,2], # std.error.ER
                                   Pr_sum$coefficients$mean[3,4], # p.value.ER
                                   Pr_vif[2] # vif.ER
        )
        matrix_to_fill[i-length(columns_ignore),] <- intermediate_SNP_data
      }
    }
    
  }
  
  return(matrix_to_fill)
}


# SENSITIVITY ANALYSIS
BetaRegression_RareVariants.sensitivity <- function(dataset, carrier_dataset){ 
  
  columns_ignore.tmp <- c("percentage", "marker", "BRIDGES_ID", "ER_statusIndex", "AgeDiagIndex", "study", "to.be.included")
  
  matrix_to_fill.tmp <- matrix(nrow = length(which(!names(dataset) %in% columns_ignore.tmp)), ncol = 15)
  
  dataset.tmp <- dataset %>%
    mutate(ER_statusIndex = case_when(ER_statusIndex == 888 ~ NA_real_, TRUE ~ ER_statusIndex)) 
  
  colnames(matrix_to_fill.tmp) <- c("gene", "estimate", "ci.lb", "ci.ub", "std.error", "p.value", "vif", "pseudo.r.squared", "loglik",
                                    "estimate.age", "ci.lb.age", "ci.ub.age", "std.error.age", "p.value.age", "vif.age")
  
  
  
  for(i in which(!names(dataset.tmp) %in% columns_ignore.tmp)){
    
    patterns <- exprs(!!sym(colnames(dataset.tmp[i])) == 1 & BRIDGES_ID %in% carrier_dataset$BRIDGES_ID ~ 1,
                      !!sym(colnames(dataset.tmp[i])) == 0 & BRIDGES_ID %in% carrier_dataset$BRIDGES_ID ~ 2,
                      !!sym(colnames(dataset.tmp[i])) == 0 & BRIDGES_ID %!in% carrier_dataset$BRIDGES_ID ~ 0)
    
    dataset.tmp.for.func <- dataset.tmp  %>%
      mutate(to.be.included = case_when(!!!patterns)) %>% # "non_carrier"
      filter(to.be.included != 2) # Here we filter out all the non-carriers that have a variant in any other gene
    
    
    possible_error_beta_regression.sensitivity <- function(dat, i){
      out <- tryCatch(exp = {
        suppressWarnings(
          Pr <- betareg(percentage ~ to.be.included + AgeDiagIndex + study, data = dat, warning = FALSE, link = "logit")
        )
        return(Pr)
      },
      error = function(e){ #message("error")
        return(c(colnames(dat[i]), rep(888, 14)))
      }, 
      warning = function(w){ #message("warning")
        suppressWarnings(
          return(c(colnames(dat[i]), rep(888, 14)))
        )}
      )
      return(out)
    }
    
    output.tmp <- possible_error_beta_regression.sensitivity(dat = dataset.tmp.for.func, i)
    
    if (class(output.tmp) == "character") { # When output = c(colnames(dosage_icogs[x]), rep(NA, 8)) meaning an error occured, then add it to the matrix right away.
      matrix_to_fill.tmp[(i-length(columns_ignore.tmp))+1,] <- output.tmp
      
    } else if (class(output.tmp) == "betareg") {# if a true answer, then put it in the confint line
      PR_ci <- confint(output.tmp, level = 0.95)
      Pr_sum <- summary(output.tmp)
      Pr_vif <- vif(output.tmp)
      intermediate_SNP_data <- c(colnames(dataset.tmp[i]), # SNP name
                                 Pr_sum$coefficients$mean[2,1], # Estimate gene
                                 PR_ci[2,1], # ci.lb gene
                                 PR_ci[2,2], # ci.ub gene
                                 Pr_sum$coefficients$mean[2,2], # std.error gene
                                 Pr_sum$coefficients$mean[2,4], # p.value gene
                                 Pr_vif[1], # vif gene
                                 
                                 Pr_sum$pseudo.r.squared, # pseudo.r.squared
                                 Pr_sum$loglik, # loglik
                                 
                                 Pr_sum$coefficients$mean[3,1], # estimate.age
                                 PR_ci[3,1], # ci.lb.age
                                 PR_ci[3,2], # ci.ub.age
                                 Pr_sum$coefficients$mean[3,2], # std.error.age
                                 Pr_sum$coefficients$mean[3,4], # p.value.age
                                 Pr_vif[2] # vif.age
      )
      matrix_to_fill.tmp[(i-length(columns_ignore.tmp))+1,] <- intermediate_SNP_data
    }
  }
  return(matrix_to_fill.tmp)
}





# RESULTS OF THE MODEL
transform_output_function <- function(dataset){
  dataset %>%
    as_tibble() %>% 
    select(gene, estimate, ci.lb, ci.ub, p.value) %>%
    mutate_at(vars(-("gene")), as.numeric) %>%
    mutate(estimate = case_when(estimate == 888 ~ NA_real_, TRUE ~ estimate), 
           ci.lb = case_when(ci.lb == 888 ~ NA_real_, TRUE ~ ci.lb), 
           ci.ub = case_when(ci.ub == 888 ~ NA_real_, TRUE ~ ci.ub), 
           p.value = case_when(p.value == 888 ~ NA_real_, TRUE ~ p.value)) %>%
    # THIS SECTIONS CREATES THE ODDS RATIO
    mutate(estimate = exp(estimate),
           ci.lb = exp(ci.lb),
           ci.ub = exp(ci.ub)) %>% 
    mutate(estimate = round(as.numeric(estimate), digits = 2),
           ci.lb = round(as.numeric(ci.lb), digits = 2),
           ci.ub = round(as.numeric(ci.ub), digits = 2),
           p.value = round(as.numeric(p.value), digits = 3))
}
