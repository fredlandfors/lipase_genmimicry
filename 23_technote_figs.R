# Clear environment
rm(list = ls())

# Get libs 
library(ggplot2)
library(ieugwasr)
library(readxl)
library(data.table)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Multiple comparisons correction ----
m_comp = 9

## LPL vs ANGPTL3, ANGPTL4, ANGPTL8; LIPG vs ANGPTL8 R59W, ANGPTL3 | LPL, ANGPTL4; LIPC vs ANGPTL4 ----
cbind_df_n <- data.table::fread(file = paste0(wd, "/data/4_lipg_Nightingale20.txt"))
cbind_df_n <- as.data.frame(cbind_df_n)
cbind_df_k <- data.table::fread(file = paste0(wd, "/data/5_lipg_Kettunen20.txt"))
cbind_df_k <- as.data.frame(cbind_df_k)
lipc_n <- as.data.frame(data.table::fread(file = paste0(wd, "/data/12_lipc_Nightingale20.txt"), sep = "\t"))
lipc_n <- lipc_n[, grep(".LIPC", names(lipc_n))]
lipc_k <- as.data.frame(data.table::fread(file = paste0(wd, "/data/16_lipc_Kettunen16.txt"), sep = "\t"))
lipc_k <- lipc_k[, grep(".LIPC", names(lipc_k))]
cbind_df_n <- cbind(cbind_df_n, lipc_n)
cbind_df_k <- cbind(cbind_df_k, lipc_k)
cbind_df_n$res_lpla3 <- lm(formula = beta.ANGPTL3 ~ beta.LPL, data = cbind_df_n)$residuals
cbind_df_k$res_lpla3 <- lm(formula = beta.ANGPTL3 ~ beta.LPL, data = cbind_df_k)$residuals

names(cbind_df_k) <- sapply(names(cbind_df_k), function(x){
  paste0(x, ".k")
})
names(cbind_df_n) <- sapply(names(cbind_df_n), function(x){
  paste0(x, ".n")
})

# Get meta data
gwas_id.n <- readxl::read_xlsx(path = "./data/gwas_id.xlsx", sheet = "Nightingale_2020")
gwas_id.k <- readxl::read_xlsx(path = "./data/gwas_id.xlsx", sheet = "Kettunen_2016")

# Bind meta data to dfs
cbind_df_n$short.id <- sapply(cbind_df_n$gwas_id.LPL.n, function(x) {
  strsplit(x, "-")[[1]][3]
})
gwas_id.n$short.id <- sapply(gwas_id.n$gwas_id, function(x) {
  strsplit(x, "-")[[1]][3]
})
cbind_df_k <- merge(
  gwas_id.k,
  cbind_df_k,
  by.x = "gwas_id",
  by.y = "gwas_id.LPL.k"
)

# Merge dfs
cbind_df_m <- merge(
  cbind_df_n,
  cbind_df_k,
  by.x = "short.id",
  by.y = "short.id2"
)

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

### av_plot ----
av.width = 2.4
av.height  = 2.4

av.hvline.color = "grey25"
av.hline.size = 0.3
av.vline.size = 0.3
av.abline.size = 0.3
av.smooth.size = 0.3
av.linerange.size = 0.3
av.axis.line = element_line(size = 0.4)
av.axis.ticks = element_line(size = 0.4)
av.plot.title = element_text(family = "Helvetica", size = 6, colour = "black", face = "bold", hjust = 0.5,
                             margin = margin(1,1,1,1))
av.plot.subtitle = element_text(family = "Helvetica", size = 5, colour = "black", face = "italic", hjust = 0.5,
                                margin = margin(1,1,2,1))
av.legend.text = element_text(family = "Helvetica", size = 5, colour = "black", face = "plain")
av.legend.key.size = unit(5, 'points')
av.legend.margin = margin(0, 0, 0, 0)
av.axis.title = element_text(family = "Helvetica", size = 5, face = "plain")
av.axis.text.x = element_text(family = "Helvetica", size = 5, colour = "black", face = "plain")
av.axis.text.y = element_text(family = "Helvetica", size = 5, colour = "black", face = "plain")
av.xlims = c(-2, 2)
av.ylims = c(-2, 2)
av.point.size = 1
av.stroke = 0.25

## functions ----
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
    xlab(paste0(x_lab, "\n[Effect size (scaled beta)]")) + 
    ylab(paste0(y_lab, "\n[Effect size (scaled beta)]")) + 
    scale_x_continuous(limits = x_lim, breaks = x_breaks) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks) +
    geom_hline(yintercept = 0, size = scatter.hline.size, color = scatter.hvline.color, linetype = "solid") + 
    geom_vline(xintercept = 0, size = scatter.vline.size, color = scatter.hvline.color, linetype = "solid") +
    geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed",
                size = scatter.abline.size) +
    # geom_linerange(aes_string(
    #   xmin = paste0("se_min.", scale_var, x_gene), 
    #   xmax = paste0("se_max.", scale_var, x_gene)
    # ), color = "grey75", size = scatter.linerange.size
    # ) +
    # geom_linerange(aes_string(
    #   ymin = paste0("se_min.", scale_var, y_gene), 
    #   ymax = paste0("se_max.", scale_var, y_gene) 
    # ), color = "grey75", size = scatter.linerange.size
    # ) +
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
      legend.position = c(1, 0.025),
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
av_plot_1 <- function(df,
                      title = "",
                      subtitle = "",
                      x_var, y_var,
                      scale_var = "TG",
                      x_gene, y_gene, s_gene,
                      x_mut, y_mut, s_mut,
                      x_rsid = "rsid.x", 
                      y_rsid = "rsid.y",
                      s_rsid = "rsid.s",
                      x_lim = c(-1, 1),
                      y_lim = c(-1, 1),
                      point_size = 2) {
  out <- ggplot(data = df, mapping = aes_string(x = x_var, y = y_var)) + 
    ggtitle(label = paste0(title),
            subtitle = paste0(subtitle)) +
    xlab(paste0(x_gene,
                "\n[Effect size (scaled beta)]")) + 
    ylab(paste0(y_gene,
                #" (", y_rsid,"-", df[1, paste0("ea.", y_gene)],
                #" ", y_mut, ")",
                "  |  ", 
                s_gene,
                #" (", s_rsid,"-", df[1, paste0("ea.", s_gene)], 
                #" ", s_mut, ")",
                "\n[Effect size (scaled beta)]")) + 
    scale_x_continuous(limits = x_lim, breaks = seq(-3, 3, 1)) +
    scale_y_continuous(limits = y_lim, breaks = seq(-3, 3, 1)) +
    geom_hline(yintercept = 0, size = av.hline.size, color = "black", linetype = "solid") + 
    geom_vline(xintercept = 0, size = av.vline.size, color = "black", linetype = "solid") +
    geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed",
                size = av.abline.size) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black",
                linetype = "dashed", size = av.smooth.size, fullrange = TRUE,
                level = 1 - 0.05 / m_comp) +
    geom_point(aes_string(fill = paste0("class.", x_gene)), 
               size = point_size, stroke = av.stroke,
               alpha = 1, shape = 21, color = "black") +
    theme_classic() +
    theme(
      plot.title = av.plot.title,
      plot.subtitle = av.plot.subtitle,
      legend.title = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1, 0.65),
      legend.key.size = av.legend.key.size,
      legend.text = av.legend.text,
      legend.margin = av.legend.margin,
      axis.title = av.axis.title,
      axis.text.x = av.axis.text.x,
      axis.text.y = av.axis.text.y,
      axis.line = av.axis.line,
      axis.ticks = av.axis.ticks
    ) +
    scale_fill_manual(breaks = c("ApoA-I", "ApoB", "VLDL", "IDL", "LDL",
                                 "HDL", "Fatty acids", "Other lipids", "Other metabolites"),
                      values = c("white", "black", "#00A08A", "#F98400", "#F2AD00",
                                 "#5BBCD6", "#FD6467", "#046C9A", "grey50")
    )
  return(out)
}

## Plots ----
cbind_df_n$scaled.beta.LPL.n <- scale(cbind_df_n$beta.LPL.n)
cbind_df_n$scaled.beta.LIPC.n <- scale(cbind_df_n$beta.LIPC.n)
cbind_df_n$scaled.beta.LIPG.n <- scale(cbind_df_n$beta.LIPG.n)
cbind_df_n$scaled.beta.ANGPTL3.n <- scale(cbind_df_n$beta.ANGPTL3.n)
cbind_df_n$scaled.beta.res.n <- resid(
  lm(scaled.beta.ANGPTL3.n ~ scaled.beta.LPL.n, data = cbind_df_n)
)

### LPL vs. LIPC ----
p1 <- scatter_plot(
  cbind_df_n,
  title = "HUMAN PLASMA: LPL vs. HL (LIPC)",
  subtitle = "Systemic effects on metabolic parameters",
  y_var = "scaled.beta.LPL.n",
  x_var = "scaled.beta.LIPC.n",
  scale_var = "",
  y_gene = "LPL.n",
  x_gene = "LIPC.n",
  y_lab = "LPL (rs115849089-A [eQTL])",
  x_lab = "LIPC (rs1800588-T [promoter variant])",
  y_lim = c(-3.5, 3.5),
  x_lim = c(-3.5, 3.5),
  y_breaks = seq(-3, 3, 1),
  x_breaks = seq(-3, 3, 1),
  point_size = scatter.point.size
) + theme(legend.position = "none")
p1

### LPL vs. LIPG ----
p2 <- scatter_plot(
  cbind_df_n,
  title = "HUMAN PLASMA: LPL vs. EL (LIPG)",
  subtitle = "Systemic effects on metabolic parameters",
  y_var = "scaled.beta.LPL.n",
  x_var = "scaled.beta.LIPG.n",
  scale_var = "",
  y_gene = "LPL.n",
  x_gene = "LIPG.n",
  y_lab = "LPL (rs115849089-A [eQTL])",
  x_lab = "LIPG (rs77960347-G [N396S])",
  y_lim = c(-3.5, 3.5),
  x_lim = c(-3.5, 3.5),
  y_breaks = seq(-3, 3, 1),
  x_breaks = seq(-3, 3, 1),
  point_size = scatter.point.size
) + theme(legend.position = "none")
p2

### LPL | ANGPTL3 vs. LIPC ----
p3 <- av_plot_1(
  cbind_df_n,
  title = "HUMAN PLASMA: LPL | ANGPTL3 vs. LIPC",
  subtitle = "Conditional systemic effects on metabolic parameters",
  y_var = "scaled.beta.res.n",
  x_var = "scaled.beta.LIPC.n",
  scale_var = "",
  y_gene = "LPL.n",
  x_gene = "LIPC.n",
  s_gene = "ANGPTL3.n",
  y_rsid = "",
  x_rsid = "",
  s_rsid = "",
  y_mut = "[eQTL]",
  x_mut = "[promoter variant]",
  s_mut = "[eQTL]",
  y_lim = c(-3.5, 3.5),
  x_lim = c(-3.5, 3.5),
  point_size = av.point.size
) + theme(legend.position = "none")
p3

### write to file ----
ggsave(plot = p1,
       path = paste0(wd, "/plots/"),
       filename = "23_technote_fig_1.pdf",
       width = scatter.width, height = scatter.height,
       units = "in", dpi = 1100)
ggsave(plot = p2,
       path = paste0(wd, "/plots/"),
       filename = "23_technote_fig_2.pdf",
       width = scatter.width, height = scatter.height,
       units = "in", dpi = 1100)
ggsave(plot = p3,
       path = paste0(wd, "/plots/"),
       filename = "23_technote_fig_3.pdf",
       width = scatter.width, height = scatter.height,
       units = "in", dpi = 1100)

# Regression analyses ----
f1 <- lm(formula = beta.LPL.n ~ beta.LIPC.n, data = cbind_df_n)
f2 <- lm(formula = beta.LPL.n ~ beta.LIPG.n, data = cbind_df_n)
f3 <- lm(formula = res_lpla3.n ~ beta.LIPC.n, data = cbind_df_n)

s1 <- summary(f1)
s2 <- summary(f2)
s3 <- summary(f3)

## Get coef confidence intervals ----
c1 <- confint(f1, level = 1 - 0.05 / m_comp)
c2 <- confint(f2, level = 1 - 0.05 / m_comp)
c3 <- confint(f3, level = 1 - 0.05 / m_comp)

# end ----