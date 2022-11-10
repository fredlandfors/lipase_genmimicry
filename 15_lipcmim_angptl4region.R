# LIPC activity vs. angptl4 region plot
library(data.table)
library(ggplot2)

# Get LIPC mimicry data ----
ssi.lipc.angptl4 <- data.table::fread(file = "./data/ssi_lipc_angptl4.txt")
ssi.lipc.angptl4 <- ssi.lipc.angptl4[rsid != ".", ]

# Get ld matrix ----
# ldm <- ieugwasr::ld_matrix(ssi_angptl4$rsid, with_alleles = FALSE, pop = "EUR")
#data.table::fwrite(as.data.frame(ldm), file = "/Users/fredriklandfors/projekt/lipase_genmimicry/data/lipc/angptl4_ldm.txt", sep = "\t")
ldm <- data.table::fread(file = "/Users/fredriklandfors/projekt/lipase_genmimicry/data/angptl4/angptl4_ldm.txt")
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
  ssi.lipc.angptl4,
  ldm_4,
  by = "rsid",
  all.x = TRUE
)
plot_df$r2.ssi = plot_df$r2.ssi * 100
plot_df$rs116843064 <- ifelse(plot_df$rsid == "rs116843064", "rs116843064-A", "")
plot_df$E40K <- ifelse(plot_df$rsid == "rs116843064", TRUE, FALSE)

## Scatter plot: R2 vs. pos ----
p1 <- ggplot(data = plot_df, aes(x = pos.ssi, y = r2.ssi)) + 
  ggtitle("HL (LIPC) activity R\U00B2 of variants in the ANGPTL4 region (19:8,229,011-8,629,011)") + 
  ylab("Explained LIPC activity [rs1800588-T] (R\U00B2 %)") +
  xlab("Chromosomal position (base pairs)") +
  geom_point(aes(color = ld.rs116843064, shape = E40K)) +
  geom_text(aes(label = rs116843064), vjust = 2, hjust = 0,
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
       path = "/Users/fredriklandfors/projekt/lipase_genmimicry/plots/sensitivity/lipc_angptl4/",
       filename = "lipc_angptl4.pdf",
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
  file = "/Users/fredriklandfors/projekt/lipase_genmimicry/plots/sensitivity/lipc_angptl4/lipca4_organismdb_track.pdf",
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

## end ----
