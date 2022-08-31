# Alpha Diversity Functions











# Set Names Beta ----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

set_names_beta <- function(x) {
  setNames(x, methods.beta)
}



# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 


# Beta Anova --------------------------------------------------------------

beta_anova <- function(beta.model, methods, names, num.cores = num.cores){
  
  # Assigns number of cores/threads for parallel computing (aka %dopar%) 
  # num.cores = 4
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # Anova
  res <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    mod <- beta.model[[beta]]
    anova(mod, by = "term") %>% 
      # tidy() %>% 
      as_tibble(rownames = "term") %>%
      rename("p.value" = "Pr(>F)" ,
             "statistic" = "F") %>%
      mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
      mutate(metric = beta, .before = 1) %>%
      arrange(desc(statistic))  # highest to lowest effect size 
  } %>%
    bind_rows()
  
  stopCluster(cl)
  return(res)
}




#Ordination data table list ---------------------------------------------
#   Description: 
#   Input: 
#   Output: 

ord_dt_list <- function(model, physeq){
  return(
    lapply(model, function(model) {
      get.biplot.data(smpls = physeq, ord = model)})
  )
  
}


# Sample Coordinate Datatable ---------------------------------------------
#   Description: 
#   Input: 
#   Output: 

samp_coord_dt <- function(ord.datatable.list, names){
  
  print(paste0("names: ", names))  # TEST
  
  sample.coord.dt <- lapply(names(ord.datatable.list), function(beta) {  # CHANGE
    ord.datatable.list[[beta]]$sample.coords[, Dist := beta]
    return(ord.datatable.list[[beta]]$sample.coords)
  }) %>% rbindlist()
  
  sample.coord.dt[
    , Dist := factor(Dist, levels = names(names)[c(1:length(names))]) ]
  
  return(sample.coord.dt)
}


# Axis Labels -------------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

axis_labels <- function(ord.datatable.list){
  return(axis.labels <- lapply(names(ord.datatable.list), function(beta) {
    data.table(
      Dist = beta,
      X.lab = ord.datatable.list[[beta]]$axes.labs[1],
      Y.lab = ord.datatable.list[[beta]]$axes.labs[2]
    )
  })  %>% rbindlist())
}



# Homogeneity of Dispersion  -------------------------------------------------
#   Description: 
#   Input: 
#   Output: Statistical results

beta_hom_disper <- function(data, physeq, vars, methods = "bray", plot = F, factor = F){
  
  tmp.list <- list()
  
  # Temp workaround in the case that data needs to be "factorized"
  if(isTRUE(factor)){
    data$vars <- factor(data[[vars]])
  }
  
  for (m in methods) {
    
    set.seed(1)
    # print(m)
    tmp.dist <- phyloseq::distance(physeq, method = m)
    # print(tmp.dist)
    tmp.mod <- betadisper(tmp.dist, data[[vars]])
    tmp.anova <- anova(tmp.mod)
    tmp.permTest <- permutest(tmp.mod, pairwise = TRUE, permutations = 999) 
    tmp.mod.HSD <- TukeyHSD(tmp.mod)
    
    cat("\n")
    print(paste0("Method: ", m))
    # print(tmp.anova)
    print(tmp.permTest) 
    # tmp.permTest %>% pander()
    # tmp.list[m] <- c(tmp.permTest)
    
    if(isTRUE(plot)){
      # plot(tmp.mod)
      boxplot(tmp.mod)
      # plot(tmp.mod.HSD)
    }
    
  }
  # tmp.dist <- phyloseq::distance(physeq, method = methods)
  # tmp.mod <- betadisper(tmp.dist, data[[vars]])
  # tmp.permTest <- permutest(tmp.mod, pairwise = TRUE, permutations = 99)
  
  # return()
}



# Beta Full dbRDA  --------------------------------------------------------------
#   Description: Distance based redundancy analysis beta diversity metrics
#   Input: variables, distance list, beta methods, method names, phyloseq obj
#   Output: statistical model


full_dbrda_model  <- function(vars, distance, methods, names, physeq, data, terms = 1, num.cores = num.cores){  
  
  # Starts parallel computing
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # Build full model
  beta.model.full <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    # progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    dist.mat <- distance[[beta]]
    form <- if(length(vars) == 1){
      paste0("dist.mat ~ ", vars)
    } else{paste0("dist.mat ~ (", paste0(vars, collapse = "+") ,")^", terms)}
    print(form)  #  TEST
    capscale(as.formula(form),
             data = data,
             na.action = na.omit # na.omit only non-missing site scores are shown
             #comm = otu.matrix(physeq)  # Error might originate here, I think
    )
  }
  
  stopCluster(cl)
  return(beta.model.full)
  
}

# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 



# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 

