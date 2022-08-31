
# postDada2Processing -----------------------------------------------------

# This script cleans up the Phyloseq object for downstream analysis
# - Cleaning meta data

# Output Start time
start.time <- Sys.time()
cat(paste0("Running postDada2Processing.R: ", start_time))

print("test")


## Load Phyloseq Object ----------------------------------------------------
# Timestamp
cat(paste0("Running postDada2Processing.R > Load phyloseq object: ", round(Sys.time() - start.time, 5), " time has elapsed" ))

ps0 <- readRDS(paste0(path.data, "/Raw/Dada2_Output/phyloseq_cleaned.rds"))



## Rarefaction -------------------------------------------------------------
# Timestamp
cat(paste0("Running postDada2Processing.R > Rarefaction: ", round(Sys.time() - start.time, 5), " time has elapsed" ))


# summary(sample_sums(ps0))  # See sample sums summary

rarefaction.minimum <- 5000 
min.smpl.size <- min(sample_sums(ps0)[sample_sums(ps0) >= rarefaction.minimum])

ps.rar <- rarefy_even_depth(
  physeq = ps0,
  sample.size = min.smpl.size,
  trimOTUs = TRUE,
  rngseed = 42
) %>%
  subset_taxa(
    !is.na(Kingdom) & 
      Kingdom != "Eukaryota" &
      Order != "Chloroplast" &
      Family != "Mitochondria"
  )
ps.rar <- rename.NA.taxa(ps.rar)

## Center Log Ratio (CLR) --------------------------------------------------
# Timestamp
cat(paste0("Running postDada2Processing.R > CLR: ", round(Sys.time() - start.time, 5), " time has elapsed" ))


clr.mat <- gen.clr.matrix(asv.mat = otu.matrix(ps0), min_reads = min.smpl.size)
ps.clr <- prune_taxa(colnames(clr.mat), ps0)
otu_table(ps.clr) <- otu_table(clr.mat, taxa_are_rows = F)

ps.clr <- ps.clr %>% subset_taxa(
  !is.na(Kingdom) & 
    Kingdom != "Eukaryota" &
    Order != "Chloroplast" &
    Family != "Mitochondria")
  
ps.clr <- rename.NA.taxa(ps.clr)



## Refactor levels ---------------------------------------------------------
# Timestamp
cat(paste0("Running postDada2Processing.R > Refactoring Levels: ", round(Sys.time() - start.time, 5), " time has elapsed" ))

ps.rar <- ps_refactor_lvls(ps.rar)
ps.clr <- ps_refactor_lvls(ps.clr)


## Export Phyloseq Objects -------------------------------------------------


saveRDS(object = ps.rar,
        file = paste0(path.input, "/phyloseq_rarefied","_", Sys.Date(),".rds")
        )

saveRDS(object = ps.clr,
        file = paste0(path.input, "/phyloseq_clr","_", Sys.Date(),".rds")
)


# Output end time of script
cat(paste0("Finished running postDada2Processing.R: ", round(Sys.time() - start.time, 5), " time has elapsed" ))
