# Load libs and set wd ----
library(data.table)
wd <- "/Volumes/LACIE_SETUP/nightingale_ukbb_gwas/"

# Endothelial lipase (EL, LIPG) ----
ssi.lipg <- data.table(fread(file = paste0(wd, "/head.txt.gz")),
                       fread(file = "/Users/fredriklandfors/projekt/lplact_mimgwas/results/lipg_rs77960347_ukbb_mimgwas.txt"))
setnames(ssi.lipg, c("rsid", "chr.ssi", "pos.ssi", "ea.ssi", "nea.ssi", "af.ssi",
                     "id.ssi", "coef.ssi", "se.ssi", "pval.ssi", "r2.ssi"))

## vs. ANGPTL4 region +-200 kb: 19:8229011-8629011 ----
# import variants in ANGPTL4 region that were GW-level associated with any of the 249 metabolic parameters
angptl4_subset <- data.table::fread(file = "/Volumes/LACIE_SETUP/nightingale_ukbb_gwas/gene_subset/angptl4_subset.txt")
data.table::setnames(angptl4_subset, new = c("metid", "rsid", "chr", "pos", "ea", "nea", "beta", "se", "lp", "eaf"))
# get names of variants
angptl4_snps <- unique(paste0(angptl4_subset$rsid, "_", angptl4_subset$chr, "_",
                              angptl4_subset$pos, "_", angptl4_subset$ea, "_",
                              angptl4_subset$nea))
# Subset metabolite-associated snps from ssi df 
ssi.lipg.angptl4  <- ssi.lipg[id.ssi %in% angptl4_snps, ]
# Calculate distance to transcription start site (tss)
ssi.lipg.angptl4$tss.dist <- 8429011 - ssi.lipg.angptl4$pos.ssi
ssi.lipg.angptl4$abs.tss.dist <- abs(ssi.lipg.angptl4$tss.dist)
# Exonic region only: 11350295 - 11352619
ssi.lipg.angptl4.exonic  <- ssi.lipg.angptl4[pos.ssi > 8429011 & pos.ssi < 8439257, ]
# Write data to file
data.table::fwrite(ssi.lipg.angptl4, file = "./data/ssi_lipg_angptl4.txt")

# Hepatic lipase (HL, LIPC) ----
ssi.lipc <- data.table(fread(file = paste0(wd, "/head.txt.gz")),
                       fread(file = "/Users/fredriklandfors/projekt/lplact_mimgwas/results/lipc_rs1800588_ukbb_mimgwas.txt"))
setnames(ssi.lipc, c("rsid", "chr.ssi", "pos.ssi", "ea.ssi", "nea.ssi", "af.ssi",
                     "id.ssi", "coef.ssi", "se.ssi", "pval.ssi", "r2.ssi"))

## vs. ANGPTL4 region +-200 kb: 19:8229011-8629011 ----
# import variants in ANGPTL4 region that were GW-level associated with any of the 249 metabolic parameters
angptl4_subset <- data.table::fread(file = "/Volumes/LACIE_SETUP/nightingale_ukbb_gwas/gene_subset/angptl4_subset.txt")
data.table::setnames(angptl4_subset, new = c("metid", "rsid", "chr", "pos", "ea", "nea", "beta", "se", "lp", "eaf"))
# get names of variants
angptl4_snps <- unique(paste0(angptl4_subset$rsid, "_", angptl4_subset$chr, "_",
                              angptl4_subset$pos, "_", angptl4_subset$ea, "_",
                              angptl4_subset$nea))
# Subset metabolite-associated snps from ssi df 
ssi.lipc.angptl4  <- ssi.lipc[id.ssi %in% angptl4_snps, ]
# Calculate distance to transcription start site (tss)
ssi.lipc.angptl4$tss.dist <- 8429011 - ssi.lipc.angptl4$pos.ssi
ssi.lipc.angptl4$abs.tss.dist <- abs(ssi.lipc.angptl4$tss.dist)
# Exonic region only: 11350295 - 11352619
ssi.lipc.angptl4.exonic  <- ssi.lipc.angptl4[pos.ssi > 8429011 & pos.ssi < 8439257, ]
# Write data to file
data.table::fwrite(ssi.lipc.angptl4, file = "./data/ssi_lipc_angptl4.txt")

# End ----