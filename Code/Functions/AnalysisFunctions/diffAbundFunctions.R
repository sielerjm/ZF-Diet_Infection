# Differential Abundance Functions






# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 

run_ancombc2 <- function(tmp.pseq, 
                         tmp.taxLevel = "Genus",
                         tmp.fixFormula = NULL, 
                         tmp.randFormula = NULL,
                         tmp.group = NULL,
                         tmp.pairwise = ifelse(!is.null(tmp.group), T, F),
                         tmp.mdfdr = ifelse(!is.null(tmp.group), list(fwer_ctrl_method = "BH", B = 100), NULL)
                         ){
  
  set.seed(42)
 
  return(ancombc2(data = tmp.pseq, tax_level = tmp.taxLevel,
                  fix_formula = tmp.fixFormula,
                  rand_formula = tmp.randFormula,
                  p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = tmp.group, #struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 8, verbose = TRUE,
                  global = ifelse(!is.null(tmp.group), T, F), pairwise = tmp.pairwise, 
                  dunnet = F, trend = F,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  mdfdr_control = tmp.mdfdr
        )
  ) 
  
}
