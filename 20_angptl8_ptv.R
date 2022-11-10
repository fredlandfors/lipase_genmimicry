# ANGPTL8 PTV analysis ----
# Clear environment
rm(list = ls())

# Get libs
library(readxl)
library(ggplot2)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Import data ----
## LPL data ----
lpl_eqtl <- data.table::fread(file = paste0(wd, "/data/1_Nightingale20.txt"))
lpl_eqtl <- as.data.frame(lpl_eqtl)
lpl_eqtl$short.id <- sapply(
  lpl_eqtl$gwas_id.LPL,
  function(x) {
    strsplit(x, "met-d-")[[1]][2]
  }
)

## ANGPTL8 PTV data ----
angptl8_ptv <- readxl::read_xls(
  path = "./data/angptl8_ptv/angptl8_ptv_metabolites.xls"
)
angptl8_ptv <- as.data.frame(angptl8_ptv)
names(angptl8_ptv) <- c("hg38_pos.angptl8_ptv", "short.id", "full.id.angptl8_ptv", 
                        "beta.angptl8_ptv", "se.angptl8_ptv", "p.angptl8_ptv")
angptl8_ptv$se_min <- angptl8_ptv$beta.angptl8_ptv - angptl8_ptv$se.angptl8_ptv
angptl8_ptv$se_max <- angptl8_ptv$beta.angptl8_ptv + angptl8_ptv$se.angptl8_ptv
scalenum.angptl8_ptv <- 0.4621327
angptl8_ptv$beta.tg.angptl8_ptv <- angptl8_ptv$beta.angptl8_ptv / scalenum.angptl8_ptv
angptl8_ptv$se.tg.angptl8_ptv <- angptl8_ptv$se.angptl8_ptv / scalenum.angptl8_ptv
angptl8_ptv$se_min.tg.angptl8_ptv <- angptl8_ptv$beta.angptl8_ptv / scalenum.angptl8_ptv - angptl8_ptv$se.angptl8_ptv / scalenum.angptl8_ptv
angptl8_ptv$se_max.tg.angptl8_ptv <- angptl8_ptv$beta.angptl8_ptv / scalenum.angptl8_ptv + angptl8_ptv$se.angptl8_ptv / scalenum.angptl8_ptv


## Merge dfs ----
cbind_df <- merge(
  x = lpl_eqtl,
  y = angptl8_ptv,
  by = "short.id",
  all.x = T,
  all.y = F
)
cbind_df$class.angptl8_ptv <- cbind_df$class.LPL

# Plot ----
## Multiple comparisons correction ----
m_comp = 9

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
scatter.xlims = c(-2, 2)
scatter.ylims = c(-2, 2)
scatter.point.size = 1
scatter.stroke = 0.25

## func ----
scatter_plot <- function(df,
                         title = "",
                         subtitle = "",
                         x_var, y_var,
                         scale_var = "TG",
                         x_gene, y_gene,
                         x_lab, y_lab,
                         x_lim = c(-1, 1),
                         y_lim = c(-1, 1),
                         point_size = 2) {
  out <- ggplot(data = df, mapping = aes_string(x = x_var, y = y_var)) + 
    ggtitle(label = paste0(title),
            subtitle = paste0(subtitle)) +
    xlab(paste0(x_lab, ")\n[1-SD effect on parameter per 1-SD TG change]")) + 
    ylab(paste0(y_lab, ")\n[1-SD effect on parameter per 1-SD TG change]")) + 
    scale_x_continuous(limits = x_lim, breaks = seq(-2, 2, 1)) +
    scale_y_continuous(limits = y_lim, breaks = seq(-2, 2, 1)) +
    geom_hline(yintercept = 0, size = scatter.hline.size, color = scatter.hvline.color, linetype = "solid") + 
    geom_vline(xintercept = 0, size = scatter.vline.size, color = scatter.hvline.color, linetype = "solid") +
    geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed",
                size = scatter.abline.size) +
    geom_linerange(aes_string(
      xmin = paste0("se_min.", scale_var, x_gene), 
      xmax = paste0("se_max.", scale_var, x_gene)
    ), color = "grey75", size = scatter.linerange.size
    ) +
    geom_linerange(aes_string(
      ymin = paste0("se_min.", scale_var, y_gene), 
      ymax = paste0("se_max.", scale_var, y_gene) 
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
      legend.position = c(1, 0.05),
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

## draw ----
### angptl4 ----
angptl8_plot <- scatter_plot(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL8 PTV",
  subtitle = "Systemic effects on metabolic parameters per 1-SD TG change",
  y_var = "beta.tg.LPL",
  x_var = "beta.tg.angptl8_ptv",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "angptl8_ptv",
  x_lab = "ANGPTL8 PTV (rs145464906-T [Q121X])",
  y_lab = "LPL (rs115849089-A [eQTL])",
  y_lim = scatter.xlims,
  x_lim = scatter.ylims,
  point_size = scatter.point.size
)
### Save ----
ggsave(plot = angptl8_plot,
       path = paste0(wd, "/plots/"),
       filename = "20_angptl8_ptv.pdf",
       width = scatter.width, height = scatter.height, units = "in",
       dpi = 1100)

# Summary stats ----
l1 <- lm(data = cbind_df, formula = beta.tg.LPL ~ beta.tg.angptl8_ptv)
s1 <- summary(l1)
c1 <- confint(l1)
