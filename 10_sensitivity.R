# Sensitivity analysis 2 ----
# Clear environment
rm(list = ls())

# install.packages("ieugwasr")
# install.packages("readxl")
library(ggplot2)
library(plot3D)

# Set working directory ----
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get gwas ids ----
gwas <- readxl::read_xlsx("./data/gwas_id.xlsx", sheet = "Nightingale_2020")

# Multiple comparisons correction ----
m_comp = 9

# # Get effect of SNP on serum triglycerides
# tg <- ieugwasr::associations(
#   c("rs1801177", "rs116843064", "rs11207977", "rs2278426", "rs77960347"),
#   "met-d-Total_TG"
# )
# # Get effect of SNP on serum cholesterol
# kol <- ieugwasr::associations(
#   c("rs1801177", "rs116843064", "rs11207977", "rs2278426", "rs77960347"),
#   "met-d-Total_C"
# )
# 
# ## Import data
# import_gwas <- function(rsid,
#                         gwas_ids,
#                         gene_name,
#                         scalevar1,
#                         scalevar2,
#                         scalenum1,
#                         scalenum2) {
#   # get data from UKBB (requires internet connection)
#   out <- ieugwasr::associations(rsid, gwas_ids)
# 
#   # define cat vars
#   out <- merge(gwas[c("gwas_id", "class")], out, by.x = "gwas_id", by.y = "id")
# 
#   # get std errors
#   out$se_min <- out$beta - out$se
#   out$se_max <- out$beta + out$se
# 
#   # scale with outcome var 1
#   out[, paste0("beta.", scalevar1)] <- out$beta / scalenum1
#   out[, paste0("se.", scalevar1)] <- out$se / scalenum1
#   out[, paste0("se_min.", scalevar1)] <- out$beta / scalenum1 - abs(out$se / scalenum1)
#   out[, paste0("se_max.", scalevar1)] <- out$beta / scalenum1 + abs(out$se / scalenum1)
# 
#   # scale with outcome var 2
#   out[, paste0("beta.", scalevar2)] <- out$beta / scalenum2
#   out[, paste0("se.", scalevar2)] <- out$se / scalenum2
#   out[, paste0("se_min.", scalevar2)] <- out$beta / scalenum2 - abs(out$se / scalenum2)
#   out[, paste0("se_max.", scalevar2)] <- out$beta / scalenum2 + abs(out$se / scalenum2)
# 
#   # rename
#   names(out) <- sapply(names(out), function(x) {paste0(x, ".", gene_name)})
# 
#   # return
#   return(out)
# }
# 
# ## Get rs1801177 [D36N] (LPL) summary stats and return to df
# lpl_nmr <- import_gwas(
#   rsid = "rs1801177",
#   gwas_ids = gwas$gwas_id,
#   gene_name = "LPL",
#   scalevar1 = "tg",
#   scalevar2 = "kol",
#   scalenum1 = -subset(tg, rsid == "rs1801177")$beta,
#   scalenum2 = -subset(kol, rsid == "rs1801177")$beta
# )
# 
# ## Get rs116843064 [E40K] (ANGPTL4) summary stats and return to df
# angptl4_nmr <- import_gwas(
#   rsid = "rs116843064",
#   gwas_ids = gwas$gwas_id,
#   gene_name = "ANGPTL4",
#   scalevar1 = "tg",
#   scalevar2 = "kol",
#   scalenum1 = -subset(tg, rsid == "rs116843064")$beta,
#   scalenum2 = -subset(kol, rsid == "rs116843064")$beta
# )
# 
# ## Get rs11207977 [eQTL] (ANGPTL3) summary stats and return to df
# angptl3_nmr <- import_gwas(
#   rsid = "rs11207977",
#   gwas_ids = gwas$gwas_id,
#   gene_name = "ANGPTL3",
#   scalevar1 = "tg",
#   scalevar2 = "kol",
#   scalenum1 = -subset(tg, rsid == "rs11207977")$beta,
#   scalenum2 = -subset(kol, rsid == "rs11207977")$beta
# )
# 
# ## Get rs2278426 [R59W] (ANGPTL8) summary stats and return to df
# angptl8_nmr <- import_gwas(
#   rsid = "rs2278426",
#   gwas_ids = gwas$gwas_id,
#   gene_name = "ANGPTL8",
#   scalevar1 = "tg",
#   scalevar2 = "kol",
#   scalenum1 = -subset(tg, rsid == "rs2278426")$beta,
#   scalenum2 = -subset(kol, rsid == "rs2278426")$beta
# )
# 
# ## Get rs77960347 [N396S] (LIPG) summary stats and return to df
# LIPG_nmr <- import_gwas(
#   rsid = "rs77960347",
#   gwas_ids = gwas$gwas_id,
#   gene_name = "LIPG",
#   scalevar1 = "tg",
#   scalevar2 = "kol",
#   scalenum1 = -subset(tg, rsid == "rs77960347")$beta,
#   scalenum2 = -subset(kol, rsid == "rs77960347")$beta
# )
# 
# ## Merge frames
# cbind_df <- cbind(lpl_nmr, angptl3_nmr, angptl4_nmr, angptl8_nmr, LIPG_nmr)
# cbind_df$class.LPL <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites",
#                              cbind_df$class.LPL)
# cbind_df$class.ANGPTL3 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites",
#                                  cbind_df$class.LPL)
# cbind_df$class.ANGPTL4 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites",
#                                  cbind_df$class.LPL)
# cbind_df$class.ANGPTL8 <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites",
#                                  cbind_df$class.LPL)
# cbind_df$class.LIPG <- ifelse(is.na(cbind_df$class.LPL), "Other metabolites",
#                               cbind_df$class.LPL)
# 
# ## Write data to file
# write.table(
#   cbind_df, file = paste0(wd, "/data/sensitivity/sensitivity2.txt"), sep = "\t",
#   quote = FALSE, row.names = FALSE
# )

# Import data ----
cbind_df <- data.table::fread(
  file = paste0(wd, "/data/sensitivity/sensitivity2.txt"), sep = "\t"
)
cbind_df <- as.data.frame(cbind_df)

# Remove scaling param
cbind_df <- subset(cbind_df, gwas_id.LPL != "met-d-Total_TG")

# LPL rs1801177-A act in the opposite direction compared to rs115849089-A
# harmonize effect estimates according to alternative allele effect
cbind_df$beta.tg.LPL = -cbind_df$beta.tg.LPL
cbind_df$se_min.tg.LPL = -cbind_df$se_min.tg.LPL
cbind_df$se_max.tg.LPL = -cbind_df$se_max.tg.LPL
cbind_df$beta.kol.LPL = -cbind_df$beta.kol.LPL
cbind_df$se_min.kol.LPL = -cbind_df$se_min.kol.LPL
cbind_df$se_max.kol.LPL = -cbind_df$se_max.kol.LPL
# LIPG rs77960347 [N396S] decreases LIPG activity
# harmonize effect estimates according to alternative allele effect
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

## av_plots ----
av.width = 2.4
av.height = 2.4

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

## threed_plot ----
threed.width = 5*1.186296
threed.height = 5
threed.cex = 0.7
threed.lims = c(-2, 2)

# 2D plots ----
scatter_plot <- function(df,
                         title = "",
                         subtitle = "",
                         x_var, y_var,
                         scale_var = "TG",
                         x_gene, y_gene,
                         x_mut, y_mut,
                         x_rsid = "rsid.x", 
                         y_rsid = "rsid.y",
                         x_lim = c(-1, 1),
                         y_lim = c(-1, 1),
                         point_size = 2) {
  out <- ggplot(data = df, mapping = aes_string(x = x_var, y = y_var)) + 
    ggtitle(label = paste0(title),
            subtitle = paste0(subtitle)) +
    xlab(paste0(x_gene," (", x_rsid,"-", df[1, paste0("ea.", x_gene)], 
                " ", x_mut, ")\n[1-SD effect on parameter per 1-SD TG change]")) + 
    ylab(paste0(y_gene," (", y_rsid,"-", df[1, paste0("ea.", y_gene)],
                " ", y_mut, ")\n[1-SD effect on parameter per 1-SD TG change]")) + 
    scale_x_continuous(limits = x_lim, breaks = seq(-3, 3, 1)) +
    scale_y_continuous(limits = y_lim, breaks = seq(-3, 3, 1)) +
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
      legend.key.size = av.legend.key.size,
      legend.text = av.legend.text,
      legend.margin = av.legend.margin,
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

## angptl4 ----
angptl4_plot <- scatter_plot(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL4",
  subtitle = "Systemic effects on metabolic parameters per 1-SD TG change",
  y_var = "beta.tg.LPL",
  x_var = "beta.tg.ANGPTL4",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "ANGPTL4",
  y_rsid = "rs1801177",
  x_rsid = "rs116843064",
  y_mut = "[D36N]",
  x_mut = "[E40K]",
  y_lim = scatter.xlims,
  x_lim = scatter.ylims,
  point_size = scatter.point.size
)

## angptl3 ----
angptl3_plot <- scatter_plot(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL3",
  subtitle = "Systemic effects on metabolic parameters per 1-SD TG change",
  y_var = "beta.tg.LPL",
  x_var = "beta.tg.ANGPTL3",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "ANGPTL3",
  y_rsid = "rs1801177",
  x_rsid = "rs11207977",
  y_mut = "[D36N]",
  x_mut = "[eQTL]",
  y_lim = scatter.xlims,
  x_lim = scatter.ylims,
  point_size = scatter.point.size
)

## angptl8 ----
angptl8_plot <- scatter_plot(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL8 R59W",
  subtitle = "Systemic effects on metabolic parameters per 1-SD TG change",
  y_var = "beta.tg.LPL",
  x_var = "beta.tg.ANGPTL8",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "ANGPTL8",
  y_rsid = "rs1801177",
  x_rsid = "rs2278426",
  y_mut = "[D36N]",
  x_mut = "[R59W]",
  y_lim = scatter.xlims,
  x_lim = scatter.ylims,
  point_size = scatter.point.size
) 

## write to file ----
ggsave(plot = angptl4_plot,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "angptl4_Nightingale20.pdf",
       width = scatter.width, height = scatter.height, units = "in", dpi = 1100)
ggsave(plot = angptl3_plot,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "angptl3_Nightingale20.pdf",
       width = scatter.width, height = scatter.height, units = "in", dpi = 1100)
ggsave(plot = angptl8_plot,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "angptl8_Nightingale20.pdf",
       width = scatter.width, height = scatter.height, units = "in", dpi = 1100)


# AV plots ----
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
                #" (", x_rsid,"-", df[1, paste0("ea.", x_gene)], 
                #" ", x_mut, ")",
                "  |  ", 
                s_gene,
                #" (", s_rsid,"-", df[1, paste0("ea.", s_gene)], 
                #" ", s_mut, ")",
                "\n[1-SD effect on metabolite per 1-SD TG change]")) + 
    ylab(paste0(y_gene,
                #" (", y_rsid,"-", df[1, paste0("ea.", y_gene)],
                #" ", y_mut, ")",
                "  |  ", 
                s_gene,
                #" (", s_rsid,"-", df[1, paste0("ea.", s_gene)], 
                #" ", s_mut, ")",
                "\n[1-SD effect on metabolite per 1-SD TG change]")) + 
    scale_x_continuous(limits = x_lim, breaks = seq(x_lim[1], x_lim[2], 0.5)) +
    scale_y_continuous(limits = y_lim, breaks = seq(y_lim[1], y_lim[2], 0.5)) +
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
      legend.position = c(1, 0.05),
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
av_plot_2 <- function(df,
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
                #" (", x_rsid,"-", df[1, paste0("ea.", x_gene)], 
                #" ", x_mut, ")",
                "  |  ", 
                s_gene,
                #" (", s_rsid,"-", df[1, paste0("ea.", s_gene)], 
                #" ", s_mut, ")",
                "\n[1-SD effect on metabolite per 1-SD TG change]")) + 
    ylab(paste0(y_gene,
                #" (", y_rsid,"-", df[1, paste0("ea.", y_gene)],
                #" ", y_mut, ")",
                "  |  ", 
                s_gene,
                #" (", s_rsid,"-", df[1, paste0("ea.", s_gene)], 
                #" ", s_mut, ")",
                "\n[1-SD effect on metabolite per 1-SD TG change]")) + 
    scale_x_continuous(limits = x_lim, breaks = seq(x_lim[1], x_lim[2], 0.5)) +
    scale_y_continuous(limits = y_lim, breaks = seq(y_lim[1], y_lim[2], 0.5)) +
    geom_hline(yintercept = 0, size = av.hline.size, color = "black", linetype = "solid") + 
    geom_vline(xintercept = 0, size = av.vline.size, color = "black", linetype = "solid") +
    geom_abline(intercept = 0, slope = -1, color = "grey50", linetype = "dashed",
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
      legend.position = c(1, 0.6),
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

## regressions ----
# y = LPL
# x = ANGPTL3
# s = ANGPTL8
res_ys <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL8, data = cbind_df)$residuals
res_xs <- lm(formula = beta.tg.ANGPTL3 ~ beta.tg.ANGPTL8, data = cbind_df)$residuals
res_yx <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3, data = cbind_df)$residuals
res_sx <- lm(formula = beta.tg.ANGPTL8 ~ beta.tg.ANGPTL3, data = cbind_df)$residuals
cbind_df <- cbind(cbind_df, res_yx, res_sx, res_ys, res_xs)

## angptl3 + angptl8 ----
av_plot_a3_a8 <- av_plot_2(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL3 | ANGPTL8",
  subtitle = "Conditional systemic effects on metabolic parameters",
  y_var = "res_ys",
  x_var = "res_xs",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "ANGPTL3",
  s_gene = "ANGPTL8",
  y_rsid = "rs1801177",
  x_rsid = "rs11207977",
  s_rsid = "rs2278426",
  y_mut = "[D36N]",
  x_mut = "[eQTL]",
  s_mut = "[R59W]",
  y_lim = av.xlims,
  x_lim = av.ylims,
  point_size = av.point.size
)

av_plot_a8_a3 <- av_plot_1(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL8 | ANGPTL3",
  subtitle = "Conditional systemic effects on metabolic parameters",
  y_var = "res_yx",
  x_var = "res_sx",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "ANGPTL8",
  s_gene = "ANGPTL3",
  y_rsid = "rs1801177",
  x_rsid = "rs2278426",
  s_rsid = "rs11207977",
  y_mut = "[D36N]",
  x_mut = "[R59W]",
  s_mut = "[eQTL]",
  y_lim = av.xlims,
  x_lim = av.ylims,
  point_size = av.point.size
)
# 
## write to file ----
ggsave(plot = av_plot_a3_a8,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "avplot1_Nightingale20.pdf",
       width = av.width, height = av.height, units = "in", dpi = 1100)
ggsave(plot = av_plot_a8_a3,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "avplot2_Nightingale20.pdf",
       width = av.width, height = av.height, units = "in", dpi = 1100)

# 3D plots ----
threed_plot <- function(x, y, z, col_grp,
                        theta, phi,
                        bty) {
  # Compute the linear regression 
  fit <- lm(z ~ x + y)
  # create a grid from the x and y values (min to max) and predict values for every point
  # this will become the regression plane
  grid.lines = 30
  x.pred <- seq(min(x), max(x), length.out = grid.lines)
  y.pred <- seq(min(y), max(y), length.out = grid.lines)
  xy <- expand.grid(x = x.pred, y = y.pred)
  z.pred <- matrix(predict(fit, newdata = xy), 
                   nrow = grid.lines, 
                   ncol = grid.lines)
  # create the fitted points for droplines to the surface
  fitpoints <- predict(fit)
  
  # scatter plot with regression plane
  plot3D::scatter3D(
    x, y, z, 
    xlim = threed.lims, 
    ylim = threed.lims,
    zlim = threed.lims,
    pch = 19, cex = threed.cex, cex.lab = threed.cex, cex.axis = threed.cex, tck = -.02,
    nticks = 7,
    ticktype = "detailed",
    
    colvar = col_grp,
    col = rev(c("white", "black", "#00A08A", "#F98400", "#F2AD00",
                "#5BBCD6", "#FD6467", "#046C9A", "grey50")),
    colkey = FALSE,
    # colkey = list(
    #   at = seq(1,9,1)*8/9 + 0.5,
    #   side = 4,
    #   addlines = TRUE,
    #   cex.axis = threed.cex,
    #   dist = -0.08,
    #   length = 0.5,
    #   width = 0.5,
    #   labels = rev(c("ApoA-I", "ApoB", "VLDL", "IDL", "LDL",
    #                  "HDL", "Fatty acids", "Other lipids", "Other metabolites"))
    # ),
    
    theta = theta, 
    phi = phi,
    bty = bty,
    
    xlab = "\nANGPTL3\n[1-SD effect per 1-SD TG change]", 
    ylab = "\n\n\nANGPTL8\n[1-SD effect per 1-SD TG change]", 
    zlab = "\n\nLPL\n[1-SD effect per 1-SD TG change]",  
    surf = list(x = x.pred, y = y.pred, z = z.pred,  
                facets = TRUE, 
                fit = fitpoints,
                col = ramp.col(
                  col = c("black","black"), 
                  n = 300,
                  alpha = 0.1
                ),
                border = "black"), 
    main = "")
}


## write to file ----
cbind_df$col_grp <- rep(NA, nrow(cbind_df))
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "ApoA-I", 9, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "ApoB", 8, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "VLDL", 7, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "IDL", 6, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "LDL", 5, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "HDL", 4, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "Fatty acids", 3, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "Other lipids", 2, cbind_df$col_grp)
cbind_df$col_grp <- ifelse(cbind_df$class.LPL == "Other metabolites", 1, cbind_df$col_grp)

### +155 degrees rotation ----
pdf(
  file = paste0(wd, "/plots/sensitivity/main_angptl3-angptl8_Nightingale20_rot155dgr.pdf"),
  width =  threed.width,
  height = threed.height
)
threed_plot(x = cbind_df$beta.tg.ANGPTL3,
            y = cbind_df$beta.tg.ANGPTL8,
            z = cbind_df$beta.tg.LPL,
            col_grp = cbind_df$col_grp,
            theta = 165,
            phi = 10,
            bty = "g")
dev.off()

# Summary stats ----
## Regression analyses ----
s1 <- summary(lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3, data = cbind_df))
s2 <- summary(lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL4, data = cbind_df))
s3 <- summary(lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL8, data = cbind_df))
s4 <- summary(lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3 + beta.tg.ANGPTL8, data = cbind_df))
s5 <- summary(lm(formula = res_ys ~ res_xs, data = cbind_df))
s6 <- summary(lm(formula = res_yx ~ res_sx, data = cbind_df))
l1 <- list(s1, s2, s3, s4, s5, s6)

## Estimate very low p values ----
# angptl3
Rmpfr::mpfr(
  2*pt(
    q = abs(s1$coefficients[2, "t value"]), 
    df = s1$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)
# angptl4
Rmpfr::mpfr(
  2*pt(
    q = abs(s2$coefficients[2, "t value"]), 
    df = s2$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)
# angptl8
Rmpfr::mpfr(
  2*pt(
    q = abs(s3$coefficients[2, "t value"]), 
    df = s3$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)
# angptl3+angptl8
Rmpfr::mpfr(
  2*pt(
    q = abs(s4$coefficients[2, "t value"]), 
    df = s4$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)
Rmpfr::mpfr(
  2*pt(
    q = abs(s4$coefficients[3, "t value"]), 
    df = s4$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)
# res_ys ~ res_xs
Rmpfr::mpfr(
  2*pt(
    q = abs(s5$coefficients[2, "t value"]), 
    df = s5$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)
# res_yx ~ res_sx
Rmpfr::mpfr(
  2*pt(
    q = abs(s6$coefficients[2, "t value"]), 
    df = s6$df[2],
    lower.tail = FALSE
  ), 
  precBits = 100
)

# Bar plot ----
## Make bar df ----
bar_df <- data.frame(
  order = factor(
    c("a", "d", "b", "c"),
    labels = c("ANGPTL3", "ANGPTL8", "ANGPTL3 + \n ANGPTL8", "ANGPTL4")
  ),
  prot = c("ANGPTL3", "ANGPTL4", "ANGPTL8", "ANGPTL3.8"),
  cohort = c("Derivation", "Derivation", "Derivation", "Derivation"),
  r2 = c(
    s1$r.squared,
    s2$r.squared,
    s3$r.squared,
    s4$r.squared
  )
)
bar_df <- rbind(bar_df, bar_df)
bar_df$r2[5:8] <- 1 - bar_df$r2[1:4]
bar_df$var_explained <- c(rep("y", 4), rep("n", 4))

## Graphical parameters ----
bar.width = 2.4
bar.height = 2.4

bar.binwidth = 0.7
bar.errorwidth = 0.4
bar.errorsize = 0.3

bar.plot.title = element_text(family = "Helvetica", size = 7, colour = "black", face = "bold", hjust = 0.5,
                              margin = margin(0,0,2,0))
bar.plot.subtitle = element_text(family = "Helvetica", size = 6, colour = "black", face = "italic", hjust = 0.5,
                                 margin = margin(0,0,0,0))
bar.axis.text.x = element_text(family = "Helvetica", size = 6, colour = "black", face = "plain",
                               margin = margin(0,0,0,0))
bar.axis.text.y = element_text(family = "Helvetica", size = 6, colour = "black", face = "plain")
bar.axis.title.x = element_blank()
bar.axis.title.y = element_text(family = "Helvetica", size = 6, face = "plain")
bar.strip.text = element_text(family = "Helvetica", size = 5, colour = "black", face = "bold")

## Subset ----
bar_df <- subset(bar_df, var_explained == "y")

## Draw plot ----
stackedbar <- ggplot(bar_df, aes(x = cohort, y = r2, fill = var_explained)) + 
  ggtitle(label = "HUMAN PLASMA: Explained variance",
          subtitle = "LPL enhancement vs. ANGPTL3-4-8 disinhibition of LPL") +
  geom_hline(yintercept = 0.25, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 0.50, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 1.00, linetype = "dotted", color = "grey75") + 
  geom_bar(stat = "identity", 
           #position = "fill",
           width = bar.binwidth, 
           color = "black", alpha = 1) + 
  scale_fill_manual(breaks = c("y", "n"), values = c("grey25", "white")) +
  xlab("") +
  ylab("Explained [genetic] LPL activity (R\U00B2, %)") +
  scale_y_continuous(labels = c(0, 25, 50, 75, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = bar.axis.title.x,
    axis.title.y = bar.axis.title.y,
    plot.subtitle = bar.plot.subtitle,
    axis.text.x = bar.axis.text.x,
    axis.text.y = bar.axis.text.x,
    plot.title = bar.plot.title,
    strip.text = bar.strip.text
  ) +
  facet_wrap(~order, strip.position = "bottom", nrow = 1) +
  guides(x = guide_axis(angle = 45), y = guide_axis(angle = 90))

## Write to file ----
ggsave(plot = stackedbar,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "stackedbar.pdf",
       width = bar.width, height = bar.height, units = "in", dpi = 1100)

# LIPG ~ ANGPTL3 + ANGPTL8 | LPL residuals analysis ----
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
    xlab(paste0(x_lab, "\n[1-SD effect on parameter 1-SD TG change]")) + 
    ylab(paste0(y_lab, "\n[1-SD effect on parameter 1-SD TG change]")) + 
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
cbind_df$res_zy.kol <- lm(formula = beta.kol.LIPG ~ beta.kol.LPL, data = cbind_df)$residuals
cbind_df$res_xy.kol <- lm(formula = beta.kol.ANGPTL3 ~ beta.kol.LPL, data = cbind_df)$residuals

## LIPG | LPL vs. (ANGPTL3 | LPL) ----
res_2 <- residual_plot(
  cbind_df,
  title = "HUMAN PLASMA: EL (LIPG) vs. ANGPTL3 | LPL",
  subtitle = "Systemic effects on metabolic parameters",
  y_var = "res_zy.kol",
  x_var = "res_xy.kol",
  class = "class.LPL",
  y_lab = "LIPG | LPL",
  x_lab = "ANGPTL3 | LPL",
  y_lim = c(-4, 4),
  x_lim = c(-4, 4),
  y_breaks = seq(-4, 4, 1),
  x_breaks = seq(-4, 4, 1),
  point_size = scatter.point.size
)

## Save ----
ggsave(plot = res_2,
       path = paste0(wd, "/plots/sensitivity/"),
       filename = "zy_xy_Nightingale20.pdf",
       width = scatter.width, height = scatter.height, units = "in",
       dpi = 1100)

# ANGPTL8 PTV ----
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

# LPL data ----
cbind_df_2 <- cbind_df
cbind_df_2$short.id <- sapply(
  cbind_df_2$gwas_id.LPL,
  function(x) {
    strsplit(x, "met-d-")[[1]][2]
  }
)

## Merge dfs ----
cbind_df_m <- merge(
  x = cbind_df_2,
  y = angptl8_ptv,
  by = "short.id",
  all.x = T,
  all.y = F
)
cbind_df_m$class.angptl8_ptv <- cbind_df_m$class.LPL

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
      legend.position = c(1, 0.65),
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
  cbind_df_m,
  title = "HUMAN PLASMA: LPL D36N vs. ANGPTL8 PTV",
  subtitle = "Systemic effects on metabolic parameters per 1-SD TG change",
  y_var = "beta.tg.LPL",
  x_var = "beta.tg.angptl8_ptv",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "angptl8_ptv",
  x_lab = "ANGPTL8 PTV (rs145464906-T [Q121X])",
  y_lab = "LPL (rs1801177-A [D36N])",
  y_lim = scatter.xlims,
  x_lim = scatter.ylims,
  point_size = scatter.point.size
)
### Save ----
ggsave(plot = angptl8_plot,
       path = paste0(wd, "/plots/sensitivity"),
       filename = "20_d36n_angptl8_ptv.pdf",
       width = scatter.width, height = scatter.height, units = "in",
       dpi = 1100)


# Regression analyses ----
f1 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3, data = cbind_df)
f2 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL4, data = cbind_df)
f3 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL8, data = cbind_df)
f4 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3 + beta.tg.ANGPTL8, data = cbind_df)
f5 <- lm(formula = beta.kol.LIPG ~ beta.kol.ANGPTL8, data = cbind_df)
f6 <- lm(formula = res_zy.kol ~ res_xy.kol, data = cbind_df)
f7 <- lm(formula = beta.tg.LPL ~ beta.tg.angptl8_ptv, data = cbind_df_m)

s1 <- summary(f1)
s2 <- summary(f2)
s3 <- summary(f3)
s4 <- summary(f4)
s5 <- summary(f5)
s6 <- summary(f6)
s7 <- summary(f7)

## Get coef confidence intervals ----
c1 <- confint(f1, level = 1 - 0.05 / m_comp)
c2 <- confint(f2, level = 1 - 0.05 / m_comp)
c3 <- confint(f3, level = 1 - 0.05 / m_comp)
c4 <- confint(f4, level = 1 - 0.05 / m_comp)
c5 <- confint(f5, level = 1 - 0.05 / m_comp)
c6 <- confint(f6, level = 1 - 0.05 / m_comp)
c7 <- confint(f7, level = 1 - 0.05 / m_comp)

# end ----