# Nightingale-Kettunen barplot
# Get libs and metadata ----
library(ieugwasr)
library(ggplot2)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Get data ----
load(file = paste0(wd, "/data/3_twosample_analysis.rData"))

# Remove scaling param
cbind_df_n <- subset(cbind_df_n, gwas_id.LPL != "met-d-Total_TG")
cbind_df_k <- subset(cbind_df_k, gwas_id.LPL != "met-c-934")

## Construct plot df ----
bar_df <- data.frame(
  order = factor(
    c("a", "d", "b", "c", "a", "d", "b", "c"),
    labels = c("ANGPTL3", "ANGPTL8", "ANGPTL3 + \n ANGPTL8", "ANGPTL4")
  ),
  prot = c("ANGPTL3", "ANGPTL4", "ANGPTL8", "ANGPTL3.8",
           "ANGPTL3", "ANGPTL4", "ANGPTL8", "ANGPTL3.8"),
  cohort = c("Derivation", "Derivation", "Derivation", "Derivation",
             "Validation", "Validation", "Validation", "Validation")
  ,
  r2 = c(
    summary(lm(data = cbind_df_n, formula = beta.tg.LPL ~ beta.tg.ANGPTL3))$r.squared,
    summary(lm(data = cbind_df_n, formula = beta.tg.LPL ~ beta.tg.ANGPTL4))$r.squared,
    summary(lm(data = cbind_df_n, formula = beta.tg.LPL ~ beta.tg.ANGPTL8))$r.squared,
    summary(lm(data = cbind_df_n, formula = beta.tg.LPL ~ beta.tg.ANGPTL3 + beta.tg.ANGPTL8))$r.squared,
    summary(lm(data = cbind_df_k, formula = beta.tg.LPL ~ beta.tg.ANGPTL3))$r.squared,
    summary(lm(data = cbind_df_k, formula = beta.tg.LPL ~ beta.tg.ANGPTL4))$r.squared,
    summary(lm(data = cbind_df_k, formula = beta.tg.LPL ~ beta.tg.ANGPTL8))$r.squared,
    summary(lm(data = cbind_df_k, formula = beta.tg.LPL ~ beta.tg.ANGPTL3 + beta.tg.ANGPTL8))$r.squared
  )
)
bar_df <- rbind(bar_df, bar_df)
bar_df$r2[9:16] <- 1 - bar_df$r2[1:8]
bar_df$var_explained <- c(rep("y", 8), rep("n", 8))

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
       path = paste0(wd, "/plots/03_twosample_analysis"),
       filename = "stackedbar.pdf",
       width = bar.width, height = bar.height, units = "in", dpi = 1100)

# end ----