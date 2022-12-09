
# Analysis Helper ---------------------------------------------------------
# Description: Instantiates important variables and settings for analyses


# Naming Conventions ------------------------------------------------------
# Description: format for naming variables/functions consistently
#   - If multiple ords per `<>`, Capitalize subsequent words (eg. <moreThanOneWord>)

# VARIABLES: 
#   <highLevelVarType>.<subtype>.<subset>.<variables>.<Modifications>
#
#   dt.<subtype>.<subset>.<modifications>  # datatables
#   df.<subtype>.<subset>.<modifications> # dataframes
#   plot.<subset>.<y-var>.<x-vars...>  # plots/figures
#   table.<subset>.<y-var>.<x-vars...>  # tables
#   mod.<subtype>.<subset>.<y>.<x-vars...>  # models (lm, glm, etc.)
#   ps.<subtype>.<subset>  # phyloseq objects

# FUNCTIONS:
#   Should be descriptive, but short enough to know the main task

# Set Environmental Variables ---------------------------------------------

# Analysis ID
analysis.ID <- paste0(
  "ZF-DietStudy_",  # Data subsetted? If so, how? "<name>_"
  Sys.Date(),  # Date of analysis
  "_rf"
)

# Number of cores
# - Automated
#   - detectCores()  # Tells you how many cores your comp has available
#   - If cores greater than 4, use ~90%. If 4 or less, use ~50%
num.cores = ifelse(detectCores() > 4, round(detectCores()*.9), round(detectCores()*.5) )
# num.cores = 1


# Import Data -------------------------------------------------------------

# Load Phyloseq Object
ps.all <- readRDS(paste0(path.input, "/phyloseq_rarefied_2022-08-17.rds"))

# Update dataframes and datatables ----------------------------------------


## Update Phyloseq object --------------------------------------------------

sample_data(ps.all) <- sample.data.frame(ps.all) %>%
  mutate(Timepoint = case_when(
    Timepoint == "3mpf" ~ "4mpf",
    Timepoint == "6mpf" ~ "7mpf"
  ))


# Load Sample Data
df.all <- sample.data.frame(ps.all)
dt.all <- sample.data.table(ps.all) 



# Update dataframes and datatables ----------------------------------------


## Update Phyloseq object --------------------------------------------------


## Update Dataframe object --------------------------------------------------


## Update Datadable object --------------------------------------------------


# Subset Data -------------------------------------------------------------
# - Different analyses will required looking at different subsets of the data
#   - All at T0: Technically all controls. 
#       - How do diets differ?
#   - Controls, T0-T1: 
#       - How do diets differ across development?
#   - T1: 
#       - How do diets differ across treatments?
#   - Exposed/Final: 
#       - How does Diet impact X, depending on Treatment/Infection?
#       - How does Infection.Status impact X?


## Initial ----------------------------------------------------------

df.T0 <- df.all %>%
  filter(Timepoint == "4mpf")

# Datatable
dt.T0 <- dt.all %>%
  filter(Timepoint == "4mpf")

# Phyloseq Object
ps.T0 <- ps.all %>%
  subset_samples(Timepoint == "4mpf")

## Initial and Final ----------------------------------------------------

# Dataframe
df.conT0T1 <- df.all %>%
  filter(PrePostExp != "Exposed" )

# Datatable
dt.conT0T1 <- dt.all %>%
  filter(PrePostExp != "Exposed" )

# Phyloseq Object
ps.conT0T1 <- ps.all %>%
  subset_samples(PrePostExp != "Exposed" )

# df.all %>%
#   filter(PrePostExp != "Exposed" ) %>%
#   group_by(PrePostExp) %>%
#   count()



## Final ----------------------------------------------------------------
# Dataframe
df.T1 <- df.all %>%
  filter(Timepoint == "7mpf")

# Datatable
dt.T1 <- dt.all %>%
  filter(Timepoint == "7mpf")

# Phyloseq Object
ps.T1 <- ps.all %>%
  subset_samples(Timepoint == "7mpf")

## Final and Control ----------------------------------------------------------------
# Dataframe
df.conT1 <- df.all %>%
  filter(Timepoint == "7mpf" & Exposure == "Unexposed")

# Datatable
dt.conT1 <- dt.all %>%
  filter(Timepoint == "7mpf" & Exposure == "Unexposed")

# Phyloseq Object
ps.conT1 <- ps.all %>%
  subset_samples(Timepoint == "7mpf" & Exposure == "Unexposed")

## Final and Exposed ------------------------------------------------
# - See what the effect of exposure had on final time point 

# Dataframe
df.expFin <- df.all %>%
  filter(Timepoint == "7mpf" & Exposure == "Exposed")

# Datatable
dt.expFin <- dt.all %>%
  filter(Timepoint == "7mpf" & Exposure == "Exposed")

# Phyloseq Object
ps.expFin <- ps.all %>%
  subset_samples(Timepoint == "7mpf" & Exposure == "Exposed")



# Diets -------------------------------------------------------------------

## Gemma -------------------------------------------------------------------

ps.Gemma <- ps.all %>%
  subset_samples(Diet == "Gemma")

ps.Gemma.con <- ps.Gemma %>%
  subset_samples(PrePostExp != "Exposed")


## Watts -------------------------------------------------------------------

ps.Watts <- ps.all %>%
  subset_samples(Diet == "Watts")

ps.Watts.con <- ps.Watts %>%
  subset_samples(PrePostExp != "Exposed")


## ZIRC --------------------------------------------------------------------

ps.ZIRC <- ps.all %>%
  subset_samples(Diet == "ZIRC")

ps.ZIRC.con <- ps.ZIRC %>%
  subset_samples(PrePostExp != "Exposed")



# Alpha-Diversity ---------------------------------------------------------

methods.alpha <- c("Observed", "Shannon", "Simpson") %>% 
  purrr::set_names()


## Calculate Alpha Scores -------------------------------------------------

# Creates a datatable of alpha diversity scores for each sample

dt.alphaScores.all <- alpha_base(ps.all,  # Phyloseq object
                                 methods.alpha,  # list of alpha methods
                                 "Sample",  # Column name for sample IDs
                                 F  # Set T if you have phylogenetic data
                                 ) 

                                 
dt.alphaScores.T0 <- alpha_base(ps.T0,  # Phyloseq object
                                methods.alpha,  # list of alpha methods
                                "Sample",  # Column name for sample IDs
                                F  # Set T if you have phylogenetic data
                                ) 

dt.alphaScores.conT0T1 <- alpha_base(ps.conT0T1,  # Phyloseq object
                                     methods.alpha,  # list of alpha methods
                                     "Sample",  # Column name for sample IDs
                                     F  # Set T if you have phylogenetic data
) 


dt.alphaScores.T1 <- alpha_base(ps.T1,  # Phyloseq object
                                methods.alpha,  # list of alpha methods
                                "Sample",  # Column name for sample IDs
                                F  # Set T if you have phylogenetic data
) 

dt.alphaScores.conT1 <- alpha_base(ps.conT1,  # Phyloseq object
                                   methods.alpha,  # list of alpha methods
                                   "Sample",  # Column name for sample IDs
                                   F  # Set T if you have phylogenetic data
) 

dt.alphaScores.expFin <- alpha_base(ps.expFin,  # Phyloseq object
                                    methods.alpha,  # list of alpha methods
                                    "Sample",  # Column name for sample IDs
                                    F  # Set T if you have phylogenetic data
) 


## Normalize Alpha Scores -------------------------------------------------

# Normalize alpha scores from 0 to 1
dt.alphaScores.norm.all <- norm_alpha_score(dt.alphaScores.all, df.all, methods.alpha)
dt.alphaScores.norm.T0 <- norm_alpha_score(dt.alphaScores.T0, df.T0, methods.alpha)
dt.alphaScores.norm.conT0T1 <- norm_alpha_score(dt.alphaScores.conT0T1, df.conT0T1, methods.alpha)
dt.alphaScores.norm.T1 <- norm_alpha_score(dt.alphaScores.T1, df.T1, methods.alpha)
dt.alphaScores.norm.conT1 <- norm_alpha_score(dt.alphaScores.conT1, df.conT1, methods.alpha)
dt.alphaScores.norm.expFin <- norm_alpha_score(dt.alphaScores.expFin, df.expFin, methods.alpha)


## Alpha Datatable -------------------------------------------------------

# Make a datatabe containing sample data and alpha diversity
dt.alphaPlus.all <- alpha_dataTable(dt.all, dt.alphaScores.norm.all)
dt.alphaPlus.T0 <- alpha_dataTable(dt.T0, dt.alphaScores.norm.T0)
dt.alphaPlus.conT0T1 <- alpha_dataTable(dt.conT0T1, dt.alphaScores.norm.conT0T1)
dt.alphaPlus.T1 <- alpha_dataTable(dt.T1, dt.alphaScores.norm.T1)
dt.alphaPlus.conT1 <- alpha_dataTable(dt.conT1, dt.alphaScores.norm.conT1)
dt.alphaPlus.expFin <- alpha_dataTable(dt.expFin, dt.alphaScores.norm.expFin)


# Melt (pivot_longer) data table for easy plotting and statistical analysis
dt.alphaPlus.all.melt <- melt_alphaDataTable(dt.alphaPlus.all)
dt.alphaPlus.T0.melt <- melt_alphaDataTable(dt.alphaPlus.T0)
dt.alphaPlus.conT0T1.melt <- melt_alphaDataTable(dt.alphaPlus.conT0T1)
dt.alphaPlus.T1.melt <- melt_alphaDataTable(dt.alphaPlus.T1)
dt.alphaPlus.conT1.melt <- melt_alphaDataTable(dt.alphaPlus.conT1)
dt.alphaPlus.expFin.melt <- melt_alphaDataTable(dt.alphaPlus.expFin)



# Beta-Diversity ----------------------------------------------------------

# Distance lists
distList.all <- gen.dist.matrices(ps.all, methods = "taxonomic", cores = num.cores)
methods.beta <- names(distList.all) %>% set_names(., .)

distList.T0 <- gen.dist.matrices(ps.T0, methods = "taxonomic", cores = num.cores)
distList.conT0T1 <- gen.dist.matrices(ps.conT0T1, methods = "taxonomic", cores = num.cores)
distList.T1 <- gen.dist.matrices(ps.T1, methods = "taxonomic", cores = num.cores)
distList.conT1 <- gen.dist.matrices(ps.conT1, methods = "taxonomic", cores = num.cores)
distList.expFin <- gen.dist.matrices(ps.expFin, methods = "taxonomic", cores = num.cores)

distList.conT0T1.Gemma <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "Gemma"), methods = "taxonomic", cores = num.cores)
distList.conT0T1.Watts <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "Watts"), methods = "taxonomic", cores = num.cores)
distList.conT0T1.ZIRC <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "ZIRC"), methods = "taxonomic", cores = num.cores)

distList.all.Gemma <- gen.dist.matrices(subset_samples(ps.all, Diet == "Gemma"), methods = "taxonomic", cores = num.cores)
distList.all.Watts <- gen.dist.matrices(subset_samples(ps.all, Diet == "Watts"), methods = "taxonomic", cores = num.cores)
distList.all.ZIRC <- gen.dist.matrices(subset_samples(ps.all, Diet == "ZIRC"), methods = "taxonomic", cores = num.cores)



# Differential Abundance --------------------------------------------------

tax.level = "Genus"


## Diet (3mo, controls) ----------------------------------------------------

diffAb.conT0 = run_ancombc2(tmp.pseq = ps.T0, 
                         tmp.taxLevel = tax.level,
                         tmp.fixFormula = paste0(c("Diet"), collapse = " + "), 
                         tmp.randFormula = NULL,
                         tmp.group = "Diet"
                          )


# output_T0 = ancombc2(data = pseq, tax_level = "Genus",
#                   fix_formula = "Diet",
#                   p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
#                   prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                   group = "Diet", #struc_zero = TRUE, neg_lb = TRUE,
#                   alpha = 0.05, n_cl = 8, verbose = TRUE,
#                   global = TRUE, pairwise = TRUE, dunnet = F, trend = F,
#                   iter_control = list(tol = 1e-2, max_iter = 20, 
#                                       verbose = TRUE),
#                   em_control = list(tol = 1e-5, max_iter = 100),
#                   # lme_control = lme4::lmerControl(),
#                   mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
# )




## Diet (6mo,  controls) ---------------------------------------------------

diffAb.conT1 = run_ancombc2(tmp.pseq = ps.conT1, 
                         tmp.taxLevel = tax.level,
                         tmp.fixFormula = paste0(c("Diet"), collapse = " + "), 
                         tmp.randFormula = NULL,
                         tmp.group = "Diet"
                         )



## Development (controls) --------------------------------------------------

# ADD INTERACTION
tmp.pseq <- ps.conT0T1
sample_data(tmp.pseq) <- microbiome::meta(tmp.pseq) %>%
  mutate(Diet.Time = paste0(Diet,".",Timepoint), .after = Age) 

# Sanity check counts are correct

microbiome::meta(tmp.pseq) %>%
  group_by(Diet.Time) %>%
  count()


diffAb.conT0T1_DietTime = run_ancombc2(tmp.pseq = tmp.pseq, 
                        tmp.taxLevel = tax.level,
                        tmp.fixFormula = "Diet.Time", 
                        tmp.randFormula = NULL,
                        tmp.group = "Diet.Time", 
                        tmp.pairwise = F,
                        tmp.mdfdr = NULL
)


diffAb.conT0T1_Time = run_ancombc2(tmp.pseq = ps.conT0T1, 
                                       tmp.taxLevel = tax.level,
                                       tmp.fixFormula = "Timepoint", 
                                       tmp.randFormula = NULL,
                                       tmp.group = NULL
)

diffAb.conT0T1_Diet = run_ancombc2(tmp.pseq = ps.conT0T1, 
                                   tmp.taxLevel = tax.level,
                                   tmp.fixFormula = paste0(c("Diet"), collapse = " + "), 
                                   tmp.randFormula = NULL,
                                   tmp.group = "Diet"
)




diffAb.conT0T1_Gemma = run_ancombc2(tmp.pseq = ps.Gemma.con, 
                                   tmp.taxLevel = tax.level,
                                   tmp.fixFormula = paste0(c("Timepoint"), collapse = " + "), 
                                   tmp.randFormula = NULL,
                                   tmp.group = NULL
)

diffAb.conT0T1_Watts = run_ancombc2(tmp.pseq = ps.Watts.con, 
                                    tmp.taxLevel = tax.level,
                                    tmp.fixFormula = paste0(c("Timepoint"), collapse = " + "), 
                                    tmp.randFormula = NULL,
                                    tmp.group = NULL
)

diffAb.conT0T1_ZIRC = run_ancombc2(tmp.pseq = ps.ZIRC.con, 
                                    tmp.taxLevel = tax.level,
                                    tmp.fixFormula = paste0(c("Timepoint"), collapse = " + "), 
                                    tmp.randFormula = NULL,
                                    tmp.group = NULL
)

# Body Condition Score

diffAb.conT0T1_DietBCS = run_ancombc2(tmp.pseq = ps.conT0T1, 
                                       tmp.taxLevel = tax.level,
                                       tmp.fixFormula = "Diet + Body.Condition.Score", 
                                       tmp.randFormula = NULL,
                                       tmp.group = "Diet", 
                                       tmp.pairwise = F,
                                       tmp.mdfdr = NULL
)

diffAb.conT0T1_ZIRCBCS = run_ancombc2(tmp.pseq = ps.ZIRC.con, 
                                       tmp.taxLevel = tax.level,
                                       tmp.fixFormula = "Body.Condition.Score", 
                                       tmp.randFormula = NULL,
                                       tmp.group = NULL, 
                                       tmp.pairwise = F,
                                       tmp.mdfdr = NULL
)

diffAb.conT0T1_ZIRCBCS_Phylum = run_ancombc2(tmp.pseq = ps.ZIRC.con, 
                                      tmp.taxLevel = "Phylum",
                                      tmp.fixFormula = "Body.Condition.Score", 
                                      tmp.randFormula = NULL,
                                      tmp.group = NULL, 
                                      tmp.pairwise = F,
                                      tmp.mdfdr = NULL
)





## Exposure ----------------------------------------------------------------

tmp.pseq <- ps.all
sample_data(tmp.pseq)$PrePostExp = factor(sample_data(tmp.pseq)$PrePostExp, levels = c("Unexposed", "Exposed", "Pre-exposure"))


diffAb.all_DietPrePostExp = run_ancombc2(tmp.pseq = tmp.pseq, 
                         tmp.taxLevel = tax.level,
                         tmp.fixFormula = paste0(c("Diet", "PrePostExp"), collapse = " + "), 
                         tmp.randFormula = NULL,
                         tmp.group = "Diet"
)

diffAb.all_PrePostExp = run_ancombc2(tmp.pseq = tmp.pseq, 
                                     tmp.taxLevel = tax.level,
                                     tmp.fixFormula = paste0(c("PrePostExp"), collapse = " + "), 
                                     tmp.randFormula = NULL,
                                     tmp.group = "PrePostExp"
)

tmp.pseq <- ps.T1
sample_data(tmp.pseq)$PrePostExp = factor(sample_data(tmp.pseq)$Exposure, levels = c("Unexposed", "Exposed"))

diffAb.T1_GemmaExp = run_ancombc2(tmp.pseq = subset_samples(tmp.pseq, Diet == "Gemma"), 
                                    tmp.taxLevel = tax.level,
                                    tmp.fixFormula = paste0(c("Exposure"), collapse = " + "), 
                                    tmp.randFormula = NULL,
                                    tmp.group = NULL
)

diffAb.T1_WattsExp = run_ancombc2(tmp.pseq = subset_samples(tmp.pseq, Diet == "Watts"), 
                                    tmp.taxLevel = tax.level,
                                    tmp.fixFormula = paste0(c("Exposure"), collapse = " + "), 
                                    tmp.randFormula = NULL,
                                    tmp.group = NULL
)

diffAb.T1_ZIRCExp = run_ancombc2(tmp.pseq = subset_samples(tmp.pseq, Diet == "ZIRC"), 
                                   tmp.taxLevel = tax.level,
                                   tmp.fixFormula = paste0(c("Exposure"), collapse = " + "), 
                                   tmp.randFormula = NULL,
                                   tmp.group = NULL
)

diffAb.T1_Exp = run_ancombc2(tmp.pseq = tmp.pseq, 
                            tmp.taxLevel = tax.level,
                            tmp.fixFormula = paste0(c("Exposure"), collapse = " + "), 
                            tmp.randFormula = NULL,
                            tmp.group = NULL
)

# ADD INTERACTION
tmp.pseq <- ps.all
sample_data(tmp.pseq) <- microbiome::meta(tmp.pseq) %>%
  mutate(Diet.Exp = paste0(Diet,".",PrePostExp), .after = Age) 

sample_data(tmp.pseq)$PrePostExp = factor(sample_data(tmp.pseq)$PrePostExp, levels = c("Unexposed", "Exposed", "Pre-exposure"))

# Sanity check counts are correct

microbiome::meta(tmp.pseq) %>%
  group_by(Diet.Exp) %>%
  count()


diffAb.all_Diet.Exp = run_ancombc2(tmp.pseq = tmp.pseq, 
                                       tmp.taxLevel = tax.level,
                                       tmp.fixFormula = "Diet.Exp", 
                                       tmp.randFormula = NULL,
                                       tmp.group = "Diet.Exp", 
                                       tmp.pairwise = F,
                                       tmp.mdfdr = NULL
)

tmp.pseq <- ps.all
sample_data(tmp.pseq)$PrePostExp = factor(sample_data(tmp.pseq)$Exposure, levels = c("Unexposed", "Exposed"))


diffAb.all_GemmaExp = run_ancombc2(tmp.pseq = subset_samples(tmp.pseq, Diet == "Gemma"), 
                                  tmp.taxLevel = tax.level,
                                  tmp.fixFormula = paste0(c("PrePostExp"), collapse = " + "), 
                                  tmp.randFormula = NULL,
                                  tmp.group = NULL
)

diffAb.all_WattsExp = run_ancombc2(tmp.pseq = subset_samples(tmp.pseq, Diet == "Watts"), 
                                  tmp.taxLevel = tax.level,
                                  tmp.fixFormula = paste0(c("PrePostExp"), collapse = " + "), 
                                  tmp.randFormula = NULL,
                                  tmp.group = NULL
)

diffAb.all_ZIRCExp = run_ancombc2(tmp.pseq = subset_samples(tmp.pseq, Diet == "ZIRC"), 
                                 tmp.taxLevel = tax.level,
                                 tmp.fixFormula = paste0(c("PrePostExp"), collapse = " + "), 
                                 tmp.randFormula = NULL,
                                 tmp.group = NULL
)


# Save R_object -----------------------------------------------------------


save_env(path.rObjects, ID = analysis.ID, extra_info = "microbiomeProcessing")


# End of doc ---------------------------------------------------------------------


