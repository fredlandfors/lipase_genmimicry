# Get libs ----
# Clear environment
rm(list = ls())

# install.packages(readxl)
# install.packages(ieugwasr)
library(ieugwasr)
library(ggplot2)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Multiple comparisons correction ----
m_comp = 9

# Get gwas ids
gwas <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Kettunen_2016")

## Get effect of SNP on total plasma cholesterol ----
tg <- rbind(
  ieugwasr::associations("rs115849089", "met-c-934"), # LPL
  ieugwasr::associations("rs116843064", "met-c-934"), # ANGPTL4
  ieugwasr::associations("rs1800588", "met-c-934") # LIPC
)

# Import data ----
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

## Get LPL rs116843064 [eQTL] ----
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

## Get rs1800588 [promoter variant] (HL/LIPC) summary stats and return to df  ---
LIPC_nmr <- import_gwas(
    rsid = "rs1800588",
    gwas_ids = gwas$gwas_id,
    gene_name = "LIPC",
    scalevar = "tg",
    scalenum = -subset(tg, rsid == "rs1800588")$beta
  )

## Merge frames ----
cbind_df <- cbind(angptl4_nmr,  LIPC_nmr, lpl_nmr)
cbind_df$class.ANGPTL4 <- ifelse(is.na(cbind_df$class.ANGPTL4), "Other metabolites",
                                 cbind_df$class.ANGPTL4)
cbind_df$class.LIPC <- ifelse(is.na(cbind_df$class.LIPC), "Other metabolites",
                                cbind_df$class.LIPC)
cbind_df$class.LPL <- ifelse(is.na(cbind_df$class.LIPC), "Other metabolites",
                              cbind_df$class.LIPC)

## Save to file ----
data.table::fwrite(cbind_df, file = paste0(wd, "/data/16_lipc_Kettunen16.txt"), sep = "\t")

## Load data ----
cbind_df <- data.table::fread(file = paste0(wd, "/data/16_lipc_Kettunen16.txt"), sep = "\t")

# Remove scaling param
cbind_df <- subset(cbind_df, gwas_id.LPL != "met-c-934")

# Plots and regressions ----
## Graphical parameters ----
### scatter_plot ----
scatter.width = 2.4
scatter.height  = 2.4

scatter.hvline.color = "grey25"
scatter.hline.size = 0.3
scatter.vline.size = 0.3
scatter.abline.size = 0.3
scatter.smooth.size = 0.3
scatter.linerange.size = 0.3
scatter.axis.line = element_line(size = 0.4)
scatter.axis.ticks = element_line(size = 0.4)
scatter.plot.title = element_text(family = "Helvetica", size = 6, colour = "black", face = "bold", hjust = 0.5,
                                  margin = margin(1,1,1,1))
scatter.plot.subtitle = element_text(family = "Helvetica", size = 5, colour = "black", face = "italic", hjust = 0.5,
                                     margin = margin(1,1,2,1))
scatter.legend.text = element_text(family = "Helvetica", size = 4, colour = "black", face = "plain")
scatter.legend.key.size = unit(0, 'points')
scatter.legend.margin = margin(0, 0, 0, 0)
scatter.axis.title = element_text(family = "Helvetica", size = 5, face = "plain")
scatter.axis.text.x = element_text(family = "Helvetica", size = 5, colour = "black", face = "plain")
scatter.axis.text.y = element_text(family = "Helvetica", size = 5, colour = "black", face = "plain")
scatter.xlims = c(-3.5, 3.5)
scatter.ylims = c(-3.5, 3.5)
scatter.point.size = 1
scatter.stroke = 0.25

## Plots ----
scatter_plot <- function(df,
                         title = "",
                         subtitle = "",
                         x_var, y_var,
                         scale_var = "",
                         x_gene, y_gene,
                         x_lab, y_lab,
                         x_lim = c(-1, 1),
                         y_lim = c(-1, 1),
                         x_breaks = seq(-1, 1, 1),
                         y_breaks = seq(-1, 1, 1),
                         point_size = 2) {
  out <- ggplot(data = df, mapping = aes_string(x = x_var, y = y_var)) + 
    ggtitle(label = paste0(title),
            subtitle = paste0(subtitle)) +
    xlab(paste0(x_lab, "\n[1-SD effect on parameter per 1-SD TG change]")) + 
    ylab(paste0(y_lab, "\n[1-SD effect on parameter per 1-SD TG change]")) + 
    scale_x_continuous(limits = x_lim, breaks = x_breaks) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks) +
    geom_hline(yintercept = 0, size = scatter.hline.size, color = scatter.hvline.color, linetype = "solid") + 
    geom_vline(xintercept = 0, size = scatter.vline.size, color = scatter.hvline.color, linetype = "solid") +
    geom_abline(intercept = 0, slope = -1, color = "grey50", linetype = "dashed",
                size = scatter.abline.size) +
    geom_linerange(aes_string(
      xmin = paste0("se_min.tg.", scale_var, x_gene), 
      xmax = paste0("se_max.tg.", scale_var, x_gene)
    ), color = "grey75", size = scatter.linerange.size
    ) +
    geom_linerange(aes_string(
      ymin = paste0("se_min.tg.", scale_var, y_gene), 
      ymax = paste0("se_max.tg.", scale_var, y_gene) 
    ), color = "grey75", size = scatter.linerange.size
    ) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black",
                linetype = "dashed", size = scatter.smooth.size, fullrange = TRUE,
                level = 1 - 0.05 / m_comp) +
    geom_point(aes_string(fill = paste0("class.", x_gene)), 
               size = point_size, stroke = scatter.stroke,
               alpha = 1, shape = 21, color = "black") +
    theme_classic() +
    theme(
      plot.title = scatter.plot.title,
      plot.subtitle = scatter.plot.subtitle,
      legend.title = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1, 0.6),
      legend.key.size = scatter.legend.key.size,
      legend.text = scatter.legend.text,
      legend.margin = scatter.legend.margin,
      axis.title = scatter.axis.title,
      axis.text.x = scatter.axis.text.x,
      axis.text.y = scatter.axis.text.y,
      axis.line = scatter.axis.line,
      axis.ticks = scatter.axis.ticks
    ) +
    scale_fill_manual(breaks = c("ApoA-I", "ApoB", "VLDL", "IDL", "LDL",
                                 "HDL", "Fatty acids", "Other lipids", "Other metabolites"),
                      values = c("white", "black", "#00A08A", "#F98400", "#F2AD00",
                                 "#5BBCD6", "#FD6467", "#046C9A", "grey50")
    ) + 
    coord_fixed(ratio = 1)
  return(out)
}

### LIPC vs. ANGPTL4 ----
p1 <- scatter_plot(
  cbind_df,
  title = "HUMAN PLASMA: HL (LIPC) vs. ANGPTL4",
  subtitle = "Systemic effects on metabolic parameters",
  y_var = "beta.tg.LIPC",
  x_var = "beta.tg.ANGPTL4",
  scale_var = "",
  y_gene = "LIPC",
  x_gene = "ANGPTL4",
  y_lab = "LIPC (rs1800588-T [promoter variant])",
  x_lab = "ANGPTL4 (rs116843064-A [E40K])",
  y_lim = c(-9, 9),
  x_lim = c(-9, 9),
  y_breaks = seq(-9, 9, 1),
  x_breaks = seq(-9, 9, 1),
  point_size = scatter.point.size
)
### write to file ----
ggsave(plot = p1,
       path = paste0(wd, "/plots/"),
       filename = "16_lipc_kettunen16.pdf",
       width = scatter.width, height = scatter.height,
       units = "in", dpi = 1100)

# Regression analyses ----
f1 <- lm(formula = beta.tg.LIPC ~ beta.tg.ANGPTL4, data = cbind_df)
s1 <- summary(f1)

## Get coef confidence intervals ----
c1 <- confint(f1, level = 1 - 0.05 / m_comp)
