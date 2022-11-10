# Clear environment
rm(list = ls())

# Get libs and metadata ----
library(ieugwasr)
library(ggplot2)
library(plot3D)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get data ----
cbind_df <- data.table::fread(file = paste0(wd, "/data/1_Nightingale20.txt"))
cbind_df <- as.data.frame(cbind_df)

# Remove scaling param
cbind_df <- subset(cbind_df, gwas_id.LPL != "met-d-Total_TG")

# Multiple comparisons correction ----
m_comp = 9

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
  y_rsid = "rs115849089",
  x_rsid = "rs116843064",
  y_mut = "[eQTL]",
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
  y_rsid = "rs115849089",
  x_rsid = "rs11207977",
  y_mut = "[eQTL]",
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
  y_rsid = "rs115849089",
  x_rsid = "rs2278426",
  y_mut = "[eQTL]",
  x_mut = "[R59W]",
  y_lim = scatter.xlims,
  x_lim = scatter.ylims,
  point_size = scatter.point.size
) 

## write to file ----
ggsave(plot = angptl4_plot,
       path = paste0(wd, "/plots/01_Nightingale20"),
       filename = "angptl4_Nightingale20.pdf",
       width = scatter.width, height = scatter.height, units = "in", dpi = 1100)
ggsave(plot = angptl3_plot,
       path = paste0(wd, "/plots/01_Nightingale20"),
       filename = "angptl3_Nightingale20.pdf",
       width = scatter.width, height = scatter.height, units = "in", dpi = 1100)
ggsave(plot = angptl8_plot,
       path = paste0(wd, "/plots/01_Nightingale20"),
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
                "\n[1-SD effect on parameter per 1-SD TG change]")) + 
    ylab(paste0(y_gene,
                #" (", y_rsid,"-", df[1, paste0("ea.", y_gene)],
                #" ", y_mut, ")",
                "  |  ", 
                s_gene,
                #" (", s_rsid,"-", df[1, paste0("ea.", s_gene)], 
                #" ", s_mut, ")",
                "\n[1-SD effect on parameter per 1-SD TG change]")) + 
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
av_plot_a3_a8 <- av_plot_1(
  cbind_df,
  title = "HUMAN PLASMA: LPL vs. ANGPTL3 | ANGPTL8",
  subtitle = "Conditional systemic effects on metabolic parameters",
  y_var = "res_ys",
  x_var = "res_xs",
  scale_var = "tg.",
  y_gene = "LPL",
  x_gene = "ANGPTL3",
  s_gene = "ANGPTL8 R59W",
  y_rsid = "rs115849089",
  x_rsid = "rs11207977",
  s_rsid = "rs2278426",
  y_mut = "[eQTL]",
  x_mut = "[eQTL]",
  s_mut = "[R59W]",
  y_lim = av.xlims,
  x_lim = av.ylims,
  point_size = av.point.size
)

## write to file ----
ggsave(plot = av_plot_a3_a8,
       path = paste0(wd, "/plots/01_Nightingale20"),
       filename = "avplot1_Nightingale20.pdf",
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
    ylab = "\n\n\nANGPTL8 R59W\n[1-SD effect per 1-SD TG change]", 
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

### -45 degrees rotation ----
pdf(
  file = paste0(wd, "/plots/3dscatters/angptl3-angptl8_Nightingale20_rot45dgr.pdf"),
  width =  threed.width, 
  height = threed.height
) 
threed_plot(x = cbind_df$beta.tg.ANGPTL3, 
            y = cbind_df$beta.tg.ANGPTL8,
            z = cbind_df$beta.tg.LPL,
            col_grp = cbind_df$col_grp,
            theta = -45,
            phi = 10,
            bty = "g")
dev.off()

### -13 degrees rotation ----
pdf(
  file = paste0(wd, "/plots/3dscatters/angptl3-angptl8_Nightingale20_rot13dgr.pdf"),
  width =  threed.width, 
  height = threed.height
) 
threed_plot(x = cbind_df$beta.tg.ANGPTL3, 
            y = cbind_df$beta.tg.ANGPTL8,
            z = cbind_df$beta.tg.LPL,
            col_grp = cbind_df$col_grp,
            theta = -13 + 180 - 12,
            phi = 10,
            bty = "g")
dev.off()

### +32 degrees rotation ----
pdf(
  file = paste0(wd, "/plots/3dscatters/angptl3-angptl8_Nightingale20_rot32dgr.pdf"),
  width =  threed.width, 
  height = threed.height
) 
threed_plot(x = cbind_df$beta.tg.ANGPTL3, 
            y = cbind_df$beta.tg.ANGPTL8,
            z = cbind_df$beta.tg.LPL,
            col_grp = cbind_df$col_grp,
            theta = -13 + 45,
            phi = 10,
            bty = "g")
dev.off()

### +155 degrees rotation ----
pdf(
  file = paste0(wd, "/plots/Nightingale20/main_angptl3-angptl8_Nightingale20_rot155dgr.pdf"),
  width =  threed.width, 
  height = threed.height
) 
threed_plot(x = cbind_df$beta.tg.ANGPTL3, 
            y = cbind_df$beta.tg.ANGPTL8,
            z = cbind_df$beta.tg.LPL,
            col_grp = cbind_df$col_grp,
            theta = 155,
            phi = 10,
            bty = "g")
dev.off()

# Summary stats ----
## Regression analyses ----
f1 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3, data = cbind_df)
f2 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL4, data = cbind_df)
f3 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL8, data = cbind_df)
f4 <- lm(formula = beta.tg.LPL ~ beta.tg.ANGPTL3 + beta.tg.ANGPTL8, data = cbind_df)
f5 <- lm(formula = res_ys ~ res_xs, data = cbind_df)
f6 <- lm(formula = res_yx ~ res_sx, data = cbind_df)

s1 <- summary(f1)
s2 <- summary(f2)
s3 <- summary(f3)
s4 <- summary(f4)
s5 <- summary(f5)
s6 <- summary(f6)

## Get coef confidence intervals ----
c1 <- confint(f1, level = 1 - 0.05 / m_comp)
c2 <- confint(f2, level = 1 - 0.05 / m_comp)
c3 <- confint(f3, level = 1 - 0.05 / m_comp)
c4 <- confint(f4, level = 1 - 0.05 / m_comp)
c5 <- confint(f5, level = 1 - 0.05 / m_comp)
c6 <- confint(f6, level = 1 - 0.05 / m_comp)