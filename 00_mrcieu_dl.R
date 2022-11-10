# Get libs
# install.packages(readxl)
# install.packages(ieugwasr)
library(ieugwasr)

# 1_Nightingale_20.R ----
# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get gwas ids
gwas <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Nightingale_2020")

## Get effect of SNP on serum triglycerides ---- 
tg <- rbind(
  ieugwasr::associations("rs115849089", "met-d-Total_TG"), # LPL
  ieugwasr::associations("rs116843064", "met-d-Total_TG"),  # ANGPTL4
  ieugwasr::associations("rs11207977", "met-d-Total_TG"),  # ANGPTL3
  ieugwasr::associations("rs2278426", "met-d-Total_TG")  # ANGPTL8
)

## Import data ----
import_gwas <- function(rsid,
                        gwas_ids,
                        gene_name,
                        scalevar,
                        scalenum) {
  # get data from UKBB (requires internet connection)
  out <- ieugwasr::associations(rsid, gwas_ids)
  
  # define cat vars 
  out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
  
  # get std errors
  out$se_min <- out$beta - out$se
  out$se_max <- out$beta + out$se
  
  # scale with outcome var
  out[, paste0("beta.", scalevar)] <- out$beta / scalenum
  out[, paste0("se.", scalevar)] <- out$se / scalenum
  out[, paste0("se_min.", scalevar)] <- out$beta / scalenum - abs(out$se / scalenum)
  out[, paste0("se_max.", scalevar)] <- out$beta / scalenum + abs(out$se / scalenum)
  
  # rename
  names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
  
  # return
  return(out)
}

## Get rs115849089 [eQTL] (LPL) summary stats and return to df ----
lpl_nmr <- import_gwas(
  rsid = "rs115849089", 
  gwas_ids = gwas$gwas_id,
  gene_name = "LPL",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs115849089")$beta
)

## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df  ----
angptl4_nmr <- import_gwas(
  rsid = "rs116843064", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL4",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs116843064")$beta
)

## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df  ----
angptl3_nmr <- import_gwas(
  rsid = "rs11207977", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL3",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs11207977")$beta
)

## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df  ----
angptl8_nmr <- import_gwas(
  rsid = "rs2278426", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL8",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs2278426")$beta
)

## Merge frames ----
cbind_df <- cbind(lpl_nmr, angptl3_nmr, angptl4_nmr, angptl8_nmr)
cbind_df$class.LPL <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                             cbind_df$class.LPL)
cbind_df$class.ANGPTL3 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL4 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL8 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)

## Write to file ----
data.table::fwrite(cbind_df, file = paste0(wd, "/data/1_Nightingale20.txt"), sep = "\t")

# 2_Kettunen16.R ----
# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get gwas ids
gwas <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Kettunen_2016")

## Get effect of SNP on serum triglycerides ----
tg <- rbind(
  ieugwasr::associations("rs115849089", "met-c-934"), # LPL
  ieugwasr::associations("rs116843064", "met-c-934"),  # ANGPTL4
  ieugwasr::associations("rs11207977", "met-c-934"),  # ANGPTL3
  ieugwasr::associations("rs2278426", "met-c-934")  # ANGPTL8
)

## Import data ----
import_gwas <- function(rsid,
                        gwas_ids,
                        gene_name,
                        scalevar,
                        scalenum) {
  # get data from UKBB (requires internet connection)
  out <- ieugwasr::associations(rsid, gwas_ids)
  
  # define cat vars 
  out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
  
  # get std errors
  out$se_min <- out$beta - out$se
  out$se_max <- out$beta + out$se
  
  # scale with outcome var
  out[, paste0("beta.", scalevar)] <- out$beta / scalenum
  out[, paste0("se.", scalevar)] <- out$se / scalenum
  out[, paste0("se_min.", scalevar)] <- out$beta / scalenum - abs(out$se / scalenum)
  out[, paste0("se_max.", scalevar)] <- out$beta / scalenum + abs(out$se / scalenum)
  
  # rename
  names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
  
  # return
  return(out)
}

## Get rs115849089 [eQTL] (LPL) summary stats and return to df ----
lpl_nmr <- import_gwas(
  rsid = "rs115849089", 
  gwas_ids = gwas$gwas_id,
  gene_name = "LPL",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs115849089")$beta
)

## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df  ----
angptl4_nmr <- import_gwas(
  rsid = "rs116843064", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL4",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs116843064")$beta
)

## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df  ----
angptl3_nmr <- import_gwas(
  rsid = "rs11207977", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL3",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs11207977")$beta
)

## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df  ----
angptl8_nmr <- import_gwas(
  rsid = "rs2278426", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL8",
  scalevar = "tg",
  scalenum = -subset(tg, rsid == "rs2278426")$beta
)

# Merge frames
cbind_df <- cbind(lpl_nmr, angptl3_nmr, angptl4_nmr, angptl8_nmr)
cbind_df$class.LPL <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                             cbind_df$class.LPL)
cbind_df$class.ANGPTL3 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL4 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL8 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)

## Write to file ----
data.table::fwrite(cbind_df, file = paste0(wd, "/data/2_Kettunen20.txt"), sep = "\t")

# 3_twosample_analysis.R ----
# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get gwas ids
gwas_k <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Kettunen_2016")

## Get effect of SNP on serum triglycerides ----
tg_k <- rbind(
  ieugwasr::associations("rs115849089", "met-c-934"), # LPL
  ieugwasr::associations("rs116843064", "met-c-934"),  # ANGPTL4
  ieugwasr::associations("rs11207977", "met-c-934"),  # ANGPTL3
  ieugwasr::associations("rs2278426", "met-c-934"),  # ANGPTL8
  ieugwasr::associations("rs77960347", "met-c-934")  # LIPG 
)

## Import data ----
import_gwas_k <- function(rsid,
                          gwas,
                          gwas_ids,
                          gene_name,
                          scalevar,
                          scalenum) {
  # get data from UKBB (requires internet connection)
  out <- ieugwasr::associations(rsid, gwas_ids)
  
  # define cat vars 
  out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
  
  # get std errors
  out$se_min <- out$beta - out$se
  out$se_max <- out$beta + out$se
  
  # scale with outcome var
  out[, paste0("beta.", scalevar)] <- out$beta / scalenum
  out[, paste0("se.", scalevar)] <- out$se / scalenum
  out[, paste0("se_min.", scalevar)] <- out$beta / scalenum - abs(out$se / scalenum)
  out[, paste0("se_max.", scalevar)] <- out$beta / scalenum + abs(out$se / scalenum)
  
  # rename
  names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
  
  # return
  return(out)
}

## Get rs115849089 [eQTL] (LPL) summary stats and return to df ----
lpl_nmr_k <- import_gwas_k(
  rsid = "rs115849089", 
  gwas = gwas_k,
  gwas_ids = gwas_k$gwas_id,
  gene_name = "LPL",
  scalevar = "tg",
  scalenum = -subset(tg_k, rsid == "rs115849089")$beta
)

## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df  ----
angptl4_nmr_k <- import_gwas_k(
  rsid = "rs116843064", 
  gwas = gwas_k,
  gwas_ids = gwas_k$gwas_id,
  gene_name = "ANGPTL4",
  scalevar = "tg",
  scalenum = -subset(tg_k, rsid == "rs116843064")$beta
)

## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df  ----
angptl3_nmr_k <- import_gwas_k(
  rsid = "rs11207977", 
  gwas = gwas_k,
  gwas_ids = gwas_k$gwas_id,
  gene_name = "ANGPTL3",
  scalevar = "tg",
  scalenum = -subset(tg_k, rsid == "rs11207977")$beta
)

## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df  ----
angptl8_nmr_k <- import_gwas_k(
  rsid = "rs2278426", 
  gwas = gwas_k,
  gwas_ids = gwas_k$gwas_id,
  gene_name = "ANGPTL8",
  scalevar = "tg",
  scalenum = -subset(tg_k, rsid == "rs2278426")$beta
)

## Get rs77960347 [eQTL] (LIPG) summary stats and return to df  ----
LIPG_nmr_k <- import_gwas_k(
  rsid = "rs77960347", 
  gwas = gwas_k,
  gwas_ids = gwas_k$gwas_id,
  gene_name = "LIPG",
  scalevar = "tg",
  scalenum = -subset(tg_k, rsid == "rs77960347")$beta
)

# Merge frames
cbind_df_k <- cbind(lpl_nmr_k, angptl3_nmr_k, angptl4_nmr_k, angptl8_nmr_k, LIPG_nmr_k)
cbind_df_k$class.LPL <- ifelse(is.na(cbind_df_k$class.LPL), "Other metabolites", 
                               cbind_df_k$class.LPL)


## Import data: Nightingale ----
gwas_n <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Nightingale_2020")

## Get effect of SNP on serum triglycerides ----
tg_n <- rbind(
  ieugwasr::associations("rs115849089", "met-d-Total_TG"), # LPL
  ieugwasr::associations("rs116843064", "met-d-Total_TG"),  # ANGPTL4
  ieugwasr::associations("rs11207977", "met-d-Total_TG"),  # ANGPTL3
  ieugwasr::associations("rs2278426", "met-d-Total_TG"),  # ANGPTL8
  ieugwasr::associations("rs77960347", "met-d-Total_TG")  # LIPG 
)

## Import data ----
import_gwas_n <- function(rsid,
                          gwas,
                          gwas_ids,
                          gene_name,
                          scalevar,
                          scalenum) {
  # get data from UKBB (requires internet connection)
  out <- ieugwasr::associations(rsid, gwas_ids)
  
  # define cat vars 
  out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
  
  # get std errors
  out$se_min <- out$beta - out$se
  out$se_max <- out$beta + out$se
  
  # scale with outcome var
  out[, paste0("beta.", scalevar)] <- out$beta / scalenum
  out[, paste0("se.", scalevar)] <- out$se / scalenum
  out[, paste0("se_min.", scalevar)] <- out$beta / scalenum - abs(out$se / scalenum)
  out[, paste0("se_max.", scalevar)] <- out$beta / scalenum + abs(out$se / scalenum)
  
  # rename
  names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
  
  # return
  return(out)
}

## Get rs115849089 [eQTL] (LPL) summary stats and return to df ----
lpl_nmr_n <- import_gwas_n(
  rsid = "rs115849089",
  gwas = gwas_n,
  gwas_ids = gwas_n$gwas_id,
  gene_name = "LPL",
  scalevar = "tg",
  scalenum = -subset(tg_n, rsid == "rs115849089")$beta
)

## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df  ----
angptl4_nmr_n <- import_gwas_n(
  rsid = "rs116843064", 
  gwas = gwas_n,
  gwas_ids = gwas_n$gwas_id,
  gene_name = "ANGPTL4",
  scalevar = "tg",
  scalenum = -subset(tg_n, rsid == "rs116843064")$beta
)

## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df  ----
angptl3_nmr_n <- import_gwas_n(
  rsid = "rs11207977", 
  gwas = gwas_n,
  gwas_ids = gwas_n$gwas_id,
  gene_name = "ANGPTL3",
  scalevar = "tg",
  scalenum = -subset(tg_n, rsid == "rs11207977")$beta
)

## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df  ----
angptl8_nmr_n <- import_gwas_n(
  rsid = "rs2278426", 
  gwas = gwas_n,
  gwas_ids = gwas_n$gwas_id,
  gene_name = "ANGPTL8",
  scalevar = "tg",
  scalenum = -subset(tg_n, rsid == "rs2278426")$beta
)

## Get rs77960347 [N396S] (LIPG) summary stats and return to df  ----
LIPG_nmr_n <- import_gwas_n(
  rsid = "rs77960347", 
  gwas = gwas_n,
  gwas_ids = gwas_n$gwas_id,
  gene_name = "LIPG",
  scalevar = "tg",
  scalenum = -subset(tg_n, rsid == "rs77960347")$beta
)

# Merge frames
cbind_df_n <- cbind(lpl_nmr_n, angptl3_nmr_n, angptl4_nmr_n, angptl8_nmr_n, LIPG_nmr_n)
cbind_df_n$class.LPL <- ifelse(is.na(cbind_df_n$class.LPL),
                               "Other metabolites", 
                               cbind_df_n$class.LPL)

## Save to image ----
save.image(file = paste0(wd, "/data/3_twosample_analysis.rData"))

# 4_lipg_Nightingale20.R ----
# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get gwas ids
gwas <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Nightingale_2020")

## Get effect of SNP on total plasma cholesterol ----
kol <- rbind(
  ieugwasr::associations("rs115849089", "met-d-Total_C"), # LPL
  ieugwasr::associations("rs116843064", "met-d-Total_C"),  # ANGPTL4
  ieugwasr::associations("rs11207977", "met-d-Total_C"),  # ANGPTL3
  ieugwasr::associations("rs2278426", "met-d-Total_C"),  # ANGPTL8
  ieugwasr::associations("rs77960347", "met-d-Total_C")  # LIPG
)


## Import data ----
import_gwas <- function(rsid,
                        gwas_ids,
                        gene_name,
                        scalevar,
                        scalenum) {
  # get data from UKBB (requires internet connection)
  out <- ieugwasr::associations(rsid, gwas_ids)
  
  # define cat vars 
  out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
  
  # get std errors
  out$se_min <- out$beta - out$se
  out$se_max <- out$beta + out$se
  
  # scale with outcome var
  out[, paste0("beta.", scalevar)] <- out$beta / scalenum
  out[, paste0("se.", scalevar)] <- out$se / scalenum
  out[, paste0("se_min.", scalevar)] <- out$beta / scalenum - abs(out$se / scalenum)
  out[, paste0("se_max.", scalevar)] <- out$beta / scalenum + abs(out$se / scalenum)
  
  # rename
  names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
  
  # return
  return(out)
}

## Get rs115849089 [eQTL] (LPL) summary stats and return to df ----
lpl_nmr <- import_gwas(
  rsid = "rs115849089", 
  gwas_ids = gwas$gwas_id,
  gene_name = "LPL",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs115849089")$beta
)

## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df  ----
angptl4_nmr <- import_gwas(
  rsid = "rs116843064", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL4",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs116843064")$beta
)

## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df  ----
angptl3_nmr <- import_gwas(
  rsid = "rs11207977", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL3",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs11207977")$beta
)

## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df  ----
angptl8_nmr <- import_gwas(
  rsid = "rs2278426", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL8",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs2278426")$beta
)

## Get rs77960347 [N396S] (LIPG) summary stats and return to df  ----
LIPG_nmr <- import_gwas(
  rsid = "rs77960347", 
  gwas_ids = gwas$gwas_id,
  gene_name = "LIPG",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs77960347")$beta
)

## Merge frames ----
cbind_df <- cbind(lpl_nmr, angptl3_nmr, angptl4_nmr, angptl8_nmr, LIPG_nmr)
cbind_df$class.LPL <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                             cbind_df$class.LPL)
cbind_df$class.ANGPTL3 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL4 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL8 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.LIPG <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                              cbind_df$class.LPL)

## Save to file ----
data.table::fwrite(cbind_df, file = paste0(wd, "/data/4_lipg_Nightingale20.txt"), sep = "\t")

# 5_lipg_Kettunen16.R ----
# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get gwas ids
gwas <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Kettunen_2016")

## Get effect of SNP on HDL cholesterol ----
kol <- rbind(
  ieugwasr::associations("rs115849089", "met-c-933"), # LPL
  ieugwasr::associations("rs116843064", "met-c-933"),  # ANGPTL4
  ieugwasr::associations("rs11207977", "met-c-933"),  # ANGPTL3
  ieugwasr::associations("rs2278426", "met-c-933"),  # ANGPTL8
  ieugwasr::associations("rs77960347", "met-c-933")  # LIPG
)

## Import data ----
import_gwas <- function(rsid,
                        gwas_ids,
                        gene_name,
                        scalevar,
                        scalenum) {
  # get data from UKBB (requires internet connection)
  out <- ieugwasr::associations(rsid, gwas_ids)
  
  # define cat vars 
  out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
  
  # get std errors
  out$se_min <- out$beta - out$se
  out$se_max <- out$beta + out$se
  
  # scale with outcome var
  out[, paste0("beta.", scalevar)] <- out$beta / scalenum
  out[, paste0("se.", scalevar)] <- out$se / scalenum
  out[, paste0("se_min.", scalevar)] <- out$beta / scalenum - abs(out$se / scalenum)
  out[, paste0("se_max.", scalevar)] <- out$beta / scalenum + abs(out$se / scalenum)
  
  # rename
  names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
  
  # return
  return(out)
}

## Get rs115849089 [eQTL] (LPL) summary stats and return to df ----
lpl_nmr <- import_gwas(
  rsid = "rs115849089", 
  gwas_ids = gwas$gwas_id,
  gene_name = "LPL",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs115849089")$beta
)

## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df  ----
angptl4_nmr <- import_gwas(
  rsid = "rs116843064", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL4",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs116843064")$beta
)

## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df  ----
angptl3_nmr <- import_gwas(
  rsid = "rs11207977", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL3",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs11207977")$beta
)

## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df  ----
angptl8_nmr <- import_gwas(
  rsid = "rs2278426", 
  gwas_ids = gwas$gwas_id,
  gene_name = "ANGPTL8",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs2278426")$beta
)

## Get rs77960347 [] (LIPG) summary stats and return to df  ----
LIPG_nmr <- import_gwas(
  rsid = "rs77960347", 
  gwas_ids = gwas$gwas_id,
  gene_name = "LIPG",
  scalevar = "kol",
  scalenum = -subset(kol, rsid == "rs77960347")$beta
)

## Merge frames ----
cbind_df <- cbind(lpl_nmr, angptl3_nmr, angptl4_nmr, angptl8_nmr, LIPG_nmr)
cbind_df$class.LPL <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                             cbind_df$class.LPL)
cbind_df$class.ANGPTL3 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL4 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.ANGPTL8 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                                 cbind_df$class.LPL)
cbind_df$class.LIPG <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites", 
                              cbind_df$class.LPL)

## Write to file ----
data.table::fwrite(cbind_df, file = paste0(wd, "/data/5_lipg_Kettunen20.txt"), sep = "\t")

# End ----
