# Clear environment
rm(list = ls())

# Get libs and metadata ----
library(ieugwasr)
library(ggplot2)
library(plot3D)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get data ----
cbind_df <- data.table::fread(file = paste0(wd, "/data/5_lipg_Kettunen20.txt"))
cbind_df <- as.data.frame(cbind_df)

# Remove scaling param
cbind_df <- subset(cbind_df, gwas_id.LPL != "met-c-933")

# Multiple comparisons correction ----
m_comp = 9

# LIPG rs77960347 [N396S] decreases LIPG activity, harmonize betas to alternative allele effect
cbind_df$beta.kol.LIPG = -cbind_df$beta.kol.LIPG
cbind_df$se_min.kol.LIPG = -cbind_df$se_min.kol.LIPG
cbind_df$se_max.kol.LIPG = -cbind_df$se_max.kol.LIPG

# Graphical parameters ----
## scatter_plot ----
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

# Scatter plot ----
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
    xlab(paste0(x_lab, "\n[1-SD effect on parameter per 1-SD TC change]")) + 
    ylab(paste0(y_lab, "\n[1-SD effect on parameter per 1-SD TC change]")) + 
    scale_x_continuous(limits = x_lim, breaks = x_breaks) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks) +
    geom_hline(yintercept = 0, size = scatter.hline.size, color = scatter.hvline.color, linetype = "solid") + 
    geom_vline(xintercept = 0, size = scatter.vline.size, color = scatter.hvline.color, linetype = "solid") +
    geom_abline(intercept = 0, slope = -1, color = "grey50", linetype = "dashed",
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

## LIPG vs. ANGPTL8 ----
EL_angptl8_plot <- scatter_plot(
  cbind_df,
  title = "HUMAN PLASMA: EL (LIPG) vs. ANGPTL8 R59W",
  subtitle = "Systemic effects on metabolic parameters",
  y_var = "beta.kol.LIPG",
  x_var = "beta.kol.ANGPTL8",
  scale_var = "kol.",
  y_gene = "LIPG",
  x_gene = "ANGPTL8",
  y_lab = "LIPG (rs77960347-G [N396S])",
  x_lab = "ANGPTL8 (rs2278426-T [R59W])",
  y_lim = c(-3, 3),
  x_lim = c(-3, 3),
  y_breaks = seq(-3, 3, 1),
  x_breaks = seq(-3, 3, 1),
  point_size = scatter.point.size
)

## Save ----
ggsave(plot = EL_angptl8_plot,
       path = paste0(wd, "/plots/05_lipg_Kettunen20"),
       filename = "lipg_a8_Kettunen16.pdf",
       width = scatter.width, height = scatter.height, units = "in",
       dpi = 1100)

# Residual plots ----
residual_plot <- function(df,
                          title = "",
                          subtitle = "",
                          x_var, y_var,
                          class,
                          x_lab,
                          y_lab,
                          x_lim = c(-1, 1),
                          y_lim = c(-1, 1),
                          x_breaks = seq(-1, 1, 1),
                          y_breaks = seq(-1, 1, 1),
                          point_size = 2) {
  out <- ggplot(data = df, mapping = aes_string(x = x_var, y = y_var)) + 
    ggtitle(label = paste0(title),
            subtitle = paste0(subtitle)) +
    xlab(paste0(x_lab, "\n[1-SD effect on parameter per 1-SD TC change]")) + 
    ylab(paste0(y_lab, "\n[1-SD effect on parameter per 1-SD TC change]")) + 
    scale_x_continuous(limits = x_lim, breaks = x_breaks) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks) +
    geom_hline(yintercept = 0, size = scatter.hline.size, color = scatter.hvline.color, linetype = "solid") + 
    geom_vline(xintercept = 0, size = scatter.vline.size, color = scatter.hvline.color, linetype = "solid") +
    geom_abline(intercept = 0, slope = -1, color = "grey50", linetype = "dashed",
                size = scatter.abline.size) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black",
                linetype = "dashed", size = scatter.smooth.size, fullrange = TRUE,
                level = 1 - 0.05 / m_comp) +
    geom_point(aes_string(fill = class), 
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

## Get residuals ----
# y = LPL
# x = ANGPTL3
# s = ANGPTL8
# z = EL (LIPG)
cbind_df$res_zy <- lm(formula = beta.kol.LIPG ~ beta.kol.LPL, data = cbind_df)$residuals
cbind_df$res_sy <- lm(formula = beta.kol.ANGPTL8 ~ beta.kol.LPL, data = cbind_df)$residuals
cbind_df$res_xy <- lm(formula = beta.kol.ANGPTL3 ~ beta.kol.LPL, data = cbind_df)$residuals

## LIPG | LPL vs. (ANGPTL3 | LPL) ----
res_2 <- residual_plot(
  cbind_df,
  title = "HUMAN PLASMA: EL (LIPG) vs. ANGPTL3 | LPL",
  subtitle = "Systemic effects on metabolic parameters",
  y_var = "res_zy",
  x_var = "res_xy",
  class = "class.LPL",
  y_lab = "LIPG | LPL",
  x_lab = "ANGPTL3 | LPL",
  y_lim = c(-3, 3),
  x_lim = c(-3, 3),
  y_breaks = seq(-3, 3, 1),
  x_breaks = seq(-3, 3, 1),
  point_size = scatter.point.size
)

## Save ----
ggsave(plot = res_2,
       path = paste0(wd, "/plots/05_lipg_Kettunen20"),
       filename = "zy_xy_Kettunen16.pdf",
       width = scatter.width, height = scatter.height, units = "in",
       dpi = 1100)

# Summary stats ----
f1 <- lm(formula = beta.kol.LIPG ~ beta.kol.ANGPTL8, data = cbind_df)
f2 <- lm(formula = res_zy ~ res_xy, data = cbind_df)
f3 <- lm(formula = beta.kol.LIPG ~ beta.kol.ANGPTL3 + beta.kol.LPL, data = cbind_df)

s1 <- summary(f1)
s2 <- summary(f2)
s3 <- summary(f3)

## Get coef confidence intervals ----
c1 <- confint(f1, level = 1 - 0.05 / m_comp)
c2 <- confint(f2, level = 1 - 0.05 / m_comp)
c3 <- confint(f3, level = 1 - 0.05 / m_comp)
