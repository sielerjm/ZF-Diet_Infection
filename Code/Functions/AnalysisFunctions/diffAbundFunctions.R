# Alpha Diversity Functions















# ANCOM-BC  -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 

run_ANCOMBC_funcs <- function(
    physeq,
    tmp.datatable,
    tmp.params
    
    ){
  
  # Genus level data
  Genus.data = aggregate_taxa(physeq, "Genus")
  
  
  # Family level data
  Family.data = aggregate_taxa(physeq, "Family")
  
  
  # Phylum level data
  Phylum.data = aggregate_taxa(physeq, "Phylum")
  
  
  
}






# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 

run_ancombc <- function(){
  out = ancombc(phyloseq = tax.level, 
                formula = "Diet", 
                p_adj_method = "BH", 
                # lib_cut = 0, 
                # prv_cut = 0,
                group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
                struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
                neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
                tol = 1e-05, 
                max_iter = 100,
                conserve = T, 
                alpha = 0.05, 
                global = T
                #assay_name = "counts"
  )
  
  res = out$res
  res_global = out$res_global
  
  return(c(res, res_global))
}





# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 

diffAbund_statsTable <- function(tmp.W, 
                                 tmp.P, 
                                 tmp.Q, 
                                 tmp.Diff, 
                                 caption.text = ""){
  
  tmp.tab <- tab_w %>%
    tibble::rownames_to_column("Taxon") %>%
    left_join(., tab_p %>% tibble::rownames_to_column("Taxon"), by = "Taxon") %>%
    left_join(., tab_q %>% tibble::rownames_to_column("Taxon"), by = "Taxon") %>%
    left_join(., tab_diff %>% tibble::rownames_to_column("Taxon"), by = "Taxon") %>%
    filter(Taxon %in% sig_taxa) %>%
    arrange(-W) %>%
    flextable() %>%
    set_caption(caption = paste0("ANCOMBC global test result (", caption.text,")")) %>%
    align(j = 2:5, align = "right") %>%
    colformat_double(j = 2, digits = 2) %>%
    colformat_double(j = 3:4, digits = 3) %>%
    merge_v(j = 1) %>%
    set_formatter(values = list("p_val" = p_val_format,
                                "q_val" = p_val_format) ) %>%
    autofit()
  
  return(tmp.tab)
  
}






# -------------------------------------------------------------------------
# Description: 
# Input: 
# Output: 
