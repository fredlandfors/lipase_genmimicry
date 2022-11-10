# Merge with UKBB clinical lipids on ~450,000 individuals ----
library(data.table)
library(ggplot2)

# Get clinical lipids data ----
wd <- "/Volumes/LACIE_SETUP/"

# GW significance threshold
gws_thresh <- -log10(0.05/12321875)

# import variants in ANGPTL4 region that were GW-level associated with any of the 249 metabolic parameters ----
night_angptl4 <- data.table::fread(file = "/Volumes/LACIE_SETUP/nightingale_ukbb_gwas/angptl4_subset/angptl4_subset.txt")
data.table::setnames(night_angptl4, new = c("metid", "rsid", "chr", "pos", "ea", "nea", "beta", "se", "lp", "eaf"))
night_angptl4 <- night_angptl4[rsid != ".", ]
angptl4_snps <- unique(night_angptl4$rsid)

# Import genome-wide mimicry analysis data ----
ssi <- data.table(fread(file = paste0(wd, "nightingale_ukbb_gwas/head.txt.gz")),
                  fread(file = "/Users/fredriklandfors/projekt/lplact_mimgwas/results/ukbb_result.txt"))
setnames(ssi, c("rsid", "chr.ssi", "pos.ssi", "ea.ssi", "nea.ssi", "af.ssi",
                "id.ssi", "coef.ssi", "se.ssi", "pval.ssi", "r2.ssi"))

# Subset GW-significant variants for any of the parameters
ssi_angptl4 <- ssi[rsid %in% angptl4_snps & rsid != ".",]

# 19:8575400-C and 19:8575400-G were both coded as rs2967574 in the data set, where C was associated with  metabolites and G
# was not. Therefore, im removing the G variant from the data set.
ssi_angptl4 <- ssi_angptl4[!(rsid == "rs2967574" & ea.ssi == "G")] 

# Get ld matrix ----
# ldm <- ieugwasr::ld_matrix(ssi_angptl4$rsid, with_alleles = FALSE, pop = "EUR")
#data.table::fwrite(as.data.frame(ldm), file = "/Users/fredriklandfors/projekt/lipase_genmimicry/data/lipc/angptl4_ldm.txt", sep = "\t")
ldm <- data.table::fread(file = "/Users/fredriklandfors/projekt/lipase_genmimicry/data/lipc/angptl4_ldm.txt")
ldm <- as.matrix(ldm)
ldm_2 <- ldm^2
ldm_2 <- data.frame(ldm_2)
ldm_3 <- ldm_2["rs116843064"]
ldm_4 <- data.frame(
  rsid = colnames(ldm_2),
  ld.rs116843064 = ldm_3$rs116843064
)

# Make plot df ----
plot_df <- merge(
  ssi_angptl4,
  ldm_4,
  by = "rsid",
  all.x = TRUE
)
plot_df$r2.ssi = plot_df$r2.ssi * 100
plot_df$rs116843064 <- ifelse(plot_df$rsid == "rs116843064", "rs116843064-A", "")
plot_df$E40K <- ifelse(plot_df$rsid == "rs116843064", TRUE, FALSE)

## Scatter plot: R2 vs. pos ----
p1 <- ggplot(data = plot_df, aes(x = pos.ssi, y = r2.ssi)) + 
  ggtitle("LPL activity R\U00B2 of variants in the ANGPTL4 region (19:8,229,011-8,629,011)") + 
  ylab("Explained [rs115849089-A] LPL activity (R\U00B2, %)") +
  xlab("Chromosomal position (base pairs)") +
  geom_point(aes(color = ld.rs116843064, shape = E40K)) +
  geom_text(aes(label = rs116843064), vjust = -0.1, hjust = -0.05,
            family = "Helvetica", colour = "black", fontface = "bold", size = 2.5) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(
    breaks = seq(8300000, 8600000, 100000),
    limits = c(8229011, 8629011),
    labels = c("8,300,000", "8,400,000", "8,500,000", "8,600,000")
  ) +
  scale_color_gradientn(
    colors = c("#00008b", "#87ceeb", "#006400", "#ffa500", "#FF0000"),
    name = "LD (r\U00B2)",
    n.breaks = 6
    ) + 
  theme_classic() + 
  theme(
    plot.title = element_text(family = "Helvetica", hjust = 0.5, size = 10,
                              colour = "black", face = "bold"),
    axis.text = element_text(family = "Helvetica", colour = "black", face = "plain", size = 8),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(family = "Helvetica", colour = "black", face = "bold", size = 8)
    )

# Write to file
ggsave(plot = p1,
       path = "/Users/fredriklandfors/projekt/lipase_genmimicry/plots/sensitivity/angptl4",
       filename = "angptl4_ldm.pdf",
       width = 5*1.5, height = 5, units = "in",
       dpi = 1100)

## Gene tracks plot ----
library(ggbio)
library(Homo.sapiens)
x1 <- "chr19:8229011-8629011"
x2 <- as(x1, "GRanges")
class(Homo.sapiens)
p.txdb <- autoplot(Homo.sapiens, which = x2, 
                   # stat = "reduce",
                   rownames.label = TRUE,
                   label.size = 1.5,
                   label.color = "black")

# Write to file
pdf(
  file = "/Users/fredriklandfors/projekt/lipase_genmimicry/plots/sensitivity/angptl4/organismdb_track.pdf",
  width = 5*1.5 * 0.7945, 
  height = 2
) 
p.txdb + 
  theme_classic(
  ) + 
  scale_x_continuous(
    breaks = seq(8300000, 8600000, 100000),
    limits = c(8229011, 8629011),
    #labels = c("8.3 Mb", "8.4 Mb", "8.5 Mb", "8.6 Mb"),
    labels = c("", "", "", "")
  ) 
dev.off()
 
# end ----
