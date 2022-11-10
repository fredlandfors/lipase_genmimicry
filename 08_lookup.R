# Lookup stats of variants on CVD risk and clinical lipids ----
# install.packages("ieugwasr")
# install.packages("TwoSampleMR")

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get available outcomes
ao <- TwoSampleMR::available_outcomes()

outcomes <- c(
  "ebi-a-GCST005194", # CAD van der Haarst 2017
  "ieu-b-111", # Clinical TG UKBB
  "ieu-b-109", # Clinical HDL UKBB
  "ieu-b-110" # Clinical LDL UKBB
)

rsids <- c(
  "rs115849089", # LPL eQTL
  "rs116843064", # ANGPTL4 E40K
  "rs11207977", # ANGPTL3 eQTL
  "rs2278426", # ANGPTL8 R59W
  "rs77960347", # LIPG N396S 
  "rs1801177", # LPL D36N 
  "rs1800588" # LIPC promoter variant
)

lookup <- ieugwasr::associations(variants = rsids, id = outcomes)
lookup2 <- merge(
  lookup,
  ao,
  by = "id",
  suffixes = c(".lookup", ".ao")
)

# Save to file ----
write.table(
  lookup2, file = paste0(wd, "/data/lookup/lookup.txt"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

# Import data ----
lookup3 <- read.delim(
  file = paste0(wd, "/data/lookup/lookup.txt"), sep = "\t", 
)


