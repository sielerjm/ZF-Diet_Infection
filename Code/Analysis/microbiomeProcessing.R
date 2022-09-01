
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

# sample_data(ps.all) <- sample.data.frame(ps.all) %>%
#   rename(Body.Condition.Score = Body.Condition.Trad)


# Load Sample Data
df.all <- sample.data.frame(ps.all)
dt.all <- sample.data.table(ps.all) 



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
  filter(Timepoint == "3mpf")

# Datatable
dt.T0 <- dt.all %>%
  filter(Timepoint == "3mpf")

# Phyloseq Object
ps.T0 <- ps.all %>%
  subset_samples(Timepoint == "3mpf")

## Initial and Final ----------------------------------------------------

# Dataframe
df.conT0T1 <- df.all %>%
  filter(Treatment == "Control")

# Datatable
dt.conT0T1 <- dt.all %>%
  filter(Treatment == "Control")

# Phyloseq Object
ps.conT0T1 <- ps.all %>%
  subset_samples(Treatment == "Control")

## Final ----------------------------------------------------------------
# Dataframe
df.T1 <- df.all %>%
  filter(Timepoint == "6mpf")

# Datatable
dt.T1 <- dt.all %>%
  filter(Timepoint == "6mpf")

# Phyloseq Object
ps.T1 <- ps.all %>%
  subset_samples(Timepoint == "6mpf")

## Final and Control ----------------------------------------------------------------
# Dataframe
df.conT1 <- df.all %>%
  filter(Timepoint == "6mpf" & Exposure == "Unexposed")

# Datatable
dt.conT1 <- dt.all %>%
  filter(Timepoint == "6mpf" & Exposure == "Unexposed")

# Phyloseq Object
ps.conT1 <- ps.all %>%
  subset_samples(Timepoint == "6mpf" & Exposure == "Unexposed")

## Final and Exposed ------------------------------------------------
# - See what the effect of exposure had on final time point 

# Dataframe
df.expFin <- df.all %>%
  filter(Timepoint == "6mpf" & Exposure == "Exposed")

# Datatable
dt.expFin <- dt.all %>%
  filter(Timepoint == "6mpf" & Exposure == "Exposed")

# Phyloseq Object
ps.expFin <- ps.all %>%
  subset_samples(Timepoint == "6mpf" & Exposure == "Exposed")


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

## Additional
distList.conT0T1.Gemma <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "Gemma"), methods = "taxonomic", cores = num.cores)
distList.conT0T1.Watts <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "Watts"), methods = "taxonomic", cores = num.cores)
distList.conT0T1.ZIRC <- gen.dist.matrices(subset_samples(ps.conT0T1, Diet == "ZIRC"), methods = "taxonomic", cores = num.cores)
distList.conT1.ZIRC <- gen.dist.matrices(subset_samples(ps.conT1, Diet == "ZIRC"), methods = "taxonomic", cores = num.cores)

# Differential Abundance ------------------------------------------------


## Diet --------------------------------------------------------------------


## Exposure ----------------------------------------------------------------

pseq = ps.all
data = dt.all

# Relevel


levels(data$PrePostExp) <- factor(levels(data$PrePostExp), levels = c("Pre-exposure", "Unexposed", "Exposed"))
data$PrePostExp <- relevel(factor(data$PrePostExp), ref = "Unexposed")

levels(sample_data(pseq)$PrePostExp) <- factor(levels(sample_data(pseq)$PrePostExp),levels = c("Pre-exposure", "Unexposed", "Exposed"))
sample_data(pseq)$PrePostExp <- relevel(factor(sample_data(pseq)$PrePostExp), ref = "Unexposed")

# ADD INTERACTION
sample_data(pseq) <- microbiome::meta(pseq) %>%
  mutate(Diet.Exp = paste0(Diet,".", PrePostExp), .after = Age) 

data <- data %>%
  mutate(Diet.Exp = paste0(Diet,".", PrePostExp), .after = Age)

# Sanity check counts are correct

microbiome::meta(pseq) %>%
  group_by(Diet.Exp) %>%
  count()


# Genus level data
Genus_data = microbiome::aggregate_taxa(pseq, "Genus")


tax.level = Genus_data #phylum_data
tax.label = "Genus"
tax.data = taxa.data.table(tax.level)


out = ancombc(phyloseq = tax.level, 
              formula = "Diet*PrePostExp", 
              p_adj_method = "BH", 
              # lib_cut = 0, 
              # prv_cut = 0,
              group = "Diet", # Set variable if "struct_zero = T". (e.g. you had two groups treated or not with antibiotics. "Antibiotics")
              # struc_zero = F,  # Set true if you expect certain treatments to select against certain taxa (e.g. antibiotics, germ-free)
              neg_lb = ifelse(tax.label == "Genus", F, T),  # F is more conservative, Set to F with ASV/Genus level, Set T with higher tax (e.g. Family, Phylum)
              tol = 1e-05, 
              max_iter = 100,
              conserve = T, 
              alpha = 0.05, 
              global = T
              #assay_name = "counts"
              )

res.EXP = out$res
res_global.EXP = out$res_global



# Global

tab_w = res_global.EXP[, "W", drop = FALSE]
tab_p = res_global.EXP[, "p_val", drop = FALSE]
tab_q = res_global.EXP[, "q_val", drop = FALSE]
tab_diff = res_global.EXP[, "diff_abn", drop = FALSE]


# Log transform

samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 

# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(abundances(tax.level) + 1) 

# Adjust the log observed abundances
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)

# Prep for visualization

sig_taxa = res_global %>%
  tibble::rownames_to_column("Taxon") %>%
  dplyr::filter(diff_abn == TRUE) %>%
  .$Taxon

df_sig.EXP = as.data.frame(t(log_obs_abn_adj[sig_taxa, ])) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(data %>%
                     dplyr::select(Sample, Diet, PrePostExp, Diet.Exp),
                   by = "Sample") %>%
  dplyr::filter(!is.na(Diet)) %>%
  tidyr::pivot_longer(cols = -one_of("Sample", "Diet", "PrePostExp", "Diet.Exp"), 
                      names_to = "Taxon", values_to = "value")


