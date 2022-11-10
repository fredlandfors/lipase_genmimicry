# Supplemental tables ----
library(data.table)
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

## Table S4 ----
cbind_df_1 <- data.table::fread(file = paste0(wd, "/data/1_Nightingale20.txt"))
cbind_df_2 <- data.table::fread(file = paste0(wd, "/data/4_lipg_Nightingale20.txt"))
cbind_df_3 <- data.table::fread(file = paste0(wd, "/data/12_lipc_Nightingale20.txt"), sep = "\t")

table_s4 <- as.data.frame(cbind(cbind_df_1, cbind_df_2, cbind_df_3))
names(table_s4) <- make.names(names(table_s4))
table_s4_fin <- table_s4[
  c("gwas_id.LPL", "trait.LPL",
    
    "rsid.LPL", "beta.LPL", "se.LPL", "p.LPL", "beta.tg.LPL", "se.tg.LPL",
    "beta.kol.LPL", "se.kol.LPL",
    
    "rsid.ANGPTL3", "beta.ANGPTL3", "se.ANGPTL3", "p.ANGPTL3", "beta.tg.ANGPTL3", "se.tg.ANGPTL3",
    "beta.kol.ANGPTL3", "se.kol.ANGPTL3",
    
    "rsid.ANGPTL4", "beta.ANGPTL4", "se.ANGPTL4", "p.ANGPTL4", "beta.tg.ANGPTL4", "se.tg.ANGPTL4",
    "beta.kol.ANGPTL4", "se.kol.ANGPTL4",
    
    "rsid.ANGPTL8", "beta.ANGPTL8", "se.ANGPTL8", "p.ANGPTL8", "beta.tg.ANGPTL8", "se.tg.ANGPTL8",
    "beta.kol.ANGPTL8", "se.kol.ANGPTL8",
    
    "rsid.LIPG", "beta.LIPG", "se.LIPG", "p.LIPG", "beta.kol.LIPG", "se.kol.LIPG",
    
    "rsid.LIPC", "beta.LIPC", "se.LIPC", "p.LIPC", "beta.tg.LIPC", "se.tg.LIPC"
    )
]
names(table_s4_fin)[1:2] <- c("mrc_ieu_id", "trait")
data.table::fwrite(table_s4_fin, file = "./tables/suppl_table_4.txt", quote = F,
                   sep = "\t", dec = ",")

## Table S5 ----
cbind_df_4 <- data.table::fread(file = paste0(wd, "/data/2_Kettunen20.txt"))
cbind_df_5 <- data.table::fread(file = paste0(wd, "/data/5_lipg_Kettunen20.txt"))
cbind_df_6 <- data.table::fread(file = paste0(wd, "/data/16_lipc_Kettunen16.txt"), sep = "\t")

table_s5 <- as.data.frame(cbind(cbind_df_4, cbind_df_5, cbind_df_6))
names(table_s5) <- make.names(names(table_s5))
table_s5_fin <- table_s5[
  c("gwas_id.LPL", "trait.LPL",
    
    "rsid.LPL", "beta.LPL", "se.LPL", "p.LPL", "beta.tg.LPL", "se.tg.LPL",
    "beta.kol.LPL", "se.kol.LPL",
    
    "rsid.ANGPTL3", "beta.ANGPTL3", "se.ANGPTL3", "p.ANGPTL3", "beta.tg.ANGPTL3", "se.tg.ANGPTL3",
    "beta.kol.ANGPTL3", "se.kol.ANGPTL3",
    
    "rsid.ANGPTL4", "beta.ANGPTL4", "se.ANGPTL4", "p.ANGPTL4", "beta.tg.ANGPTL4", "se.tg.ANGPTL4",
    "beta.kol.ANGPTL4", "se.kol.ANGPTL4",
    
    "rsid.ANGPTL8", "beta.ANGPTL8", "se.ANGPTL8", "p.ANGPTL8", "beta.tg.ANGPTL8", "se.tg.ANGPTL8",
    "beta.kol.ANGPTL8", "se.kol.ANGPTL8",
    
    "rsid.LIPG", "beta.LIPG", "se.LIPG", "p.LIPG", "beta.kol.LIPG", "se.kol.LIPG",
    
    "rsid.LIPC", "beta.LIPC", "se.LIPC", "p.LIPC", "beta.tg.LIPC", "se.tg.LIPC"
  )
]
names(table_s5_fin)[1:2] <- c("mrc_ieu_id", "trait")
data.table::fwrite(table_s5_fin, file = "./tables/suppl_table_5.txt", quote = F,
                   sep = "\t", dec = ",")

# End ----
# View(cbind(cbind_df_1[, "gwas_id.LPL"], cbind_df_2[, "gwas_id.LPL"], cbind_df_3[, "gwas_id.LPL"]))
# View(cbind(cbind_df_4[, "gwas_id.LPL"], cbind_df_5[, "gwas_id.LPL"], cbind_df_6[, "gwas_id.LPL"]))