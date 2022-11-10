# 22_crossover.R ----
# Clear environment
rm(list = ls())

# Get libs 
library(ggplot2)
library(ieugwasr)
library(readxl)
library(data.table)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

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

# # Get residuals for conditional analyses
# cbind_df_m$res_lpla3.n <- lm(formula = beta.ANGPTL3.n ~ beta.LPL.n, data = cbind_df_m)$residuals
# cbind_df_m$res_lpla3.k <- lm(formula = beta.ANGPTL3.k ~ beta.LPL.k, data = cbind_df_m)$residuals
# cbind_df_n$res_lpla3.n <- lm(formula = beta.ANGPTL3.n ~ beta.LPL.n, data = cbind_df_n)$residuals
# cbind_df_k$res_lpla3.k <- lm(formula = beta.ANGPTL3.k ~ beta.LPL.k, data = cbind_df_k)$residuals

## Make barplot dfs ----
### UKB (y) vs. Kett (x) ----
bar_df_1 <- data.frame(
  reg = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL vs.\nANGPTL3", "LPL vs.\nANGPTL4", "LPL vs.\nANGPTL8 R59W", "LPL vs.\nANGPTL3 +\nANGPTL8 R59W", 
               "EL vs.\nANGPTL8 R59W", "EL vs.\nANGPTL3 | LPL", "EL vs.\nANGPTL4", "HL vs.\nANGPTL4")
  ),
  enz = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL", "LPL", "LPL", "LPL", "EL", "EL", "EL", "HL")
  ),
  cohort = factor(
    rep("c", 8),
    labels = "Deriv. (y) vs. valid. (x)"
  ),
  r2 = c(
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL3.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL4.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL8.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL3.k + beta.ANGPTL8.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ beta.ANGPTL8.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ res_lpla3.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ beta.ANGPTL4.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPC.n ~ beta.ANGPTL4.k, data = cbind_df_m))$r.squared
  )
)

### Kett (y) vs. UKB (x) ----
bar_df_2 <- data.frame(
  reg = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL vs.\nANGPTL3", "LPL vs.\nANGPTL4", "LPL vs.\nANGPTL8 R59W", "LPL vs.\nANGPTL3 +\nANGPTL8 R59W", 
               "EL vs.\nANGPTL8 R59W", "EL vs.\nANGPTL3 | LPL", "EL vs.\nANGPTL4", "HL vs.\nANGPTL4")
  ),
  enz = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL", "LPL", "LPL", "LPL", "EL", "EL", "EL", "HL")
  ),
  cohort = factor(
    rep("c", 8),
    labels = "Valid. (y) vs. deriv. (x)"
  ),
  r2 = c(
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL3.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL4.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL8.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL3.n + beta.ANGPTL8.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.k ~ beta.ANGPTL8.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.k ~ res_lpla3.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.k ~ beta.ANGPTL4.n, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPC.k ~ beta.ANGPTL4.n, data = cbind_df_m))$r.squared
  )
)

### Derivation ----
bar_df_3 <- data.frame(
  reg = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL vs.\nANGPTL3", "LPL vs.\nANGPTL4", "LPL vs.\nANGPTL8 R59W", "LPL vs.\nANGPTL3 +\nANGPTL8 R59W", 
               "EL vs.\nANGPTL8 R59W", "EL vs.\nANGPTL3 | LPL", "EL vs.\nANGPTL4", "HL vs.\nANGPTL4")
  ),
  enz = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL", "LPL", "LPL", "LPL", "EL", "EL", "EL", "HL")
  ),
  cohort = factor(
    rep("a", 8),
    labels = "Deriv."
  ),
  r2 = c(
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL3.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL4.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL8.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LPL.n ~ beta.ANGPTL3.n + beta.ANGPTL8.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ beta.ANGPTL8.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ res_lpla3.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ beta.ANGPTL4.n, data = cbind_df_n))$r.squared,
    summary(lm(formula = beta.LIPC.n ~ beta.ANGPTL4.n, data = cbind_df_n))$r.squared
  )
)

### Validation ----
bar_df_4 <- data.frame(
  reg = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL vs.\nANGPTL3", "LPL vs.\nANGPTL4", "LPL vs.\nANGPTL8 R59W", "LPL vs.\nANGPTL3 +\nANGPTL8 R59W", 
               "EL vs.\nANGPTL8 R59W", "EL vs.\nANGPTL3 | LPL", "EL vs.\nANGPTL4", "HL vs.\nANGPTL4")
  ),
  enz = factor(
    c("a", "b", "c", "d", "e", "f", "g", "h"),
    labels = c("LPL", "LPL", "LPL", "LPL", "EL", "EL", "EL", "HL")
  ),
  cohort = factor(
    rep("b", 8),
    labels = "Valid."
  ),
  r2 = c(
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL3.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL4.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL8.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LPL.k ~ beta.ANGPTL3.k + beta.ANGPTL8.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LIPG.k ~ beta.ANGPTL8.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LIPG.k ~ res_lpla3.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LIPG.k ~ beta.ANGPTL4.k, data = cbind_df_k))$r.squared,
    summary(lm(formula = beta.LIPC.k ~ beta.ANGPTL4.k, data = cbind_df_k))$r.squared
  )
)

### Bind together ----
bar_df <- rbind(bar_df_3, bar_df_4, bar_df_1, bar_df_2)


## Draw plot ----
### Graphical parameters ----
bar.width = 2.4
bar.height = 2.4
bar.binwidth = 0.7
bar.errorwidth = 0.4
bar.errorsize = 0.3
bar.plot.title = element_text(family = "Helvetica", size = 8, colour = "black", face = "bold", hjust = 0.5,
                              margin = margin(0,0,2,0))
bar.plot.subtitle = element_text(family = "Helvetica", size = 7, colour = "black", face = "italic", hjust = 0.5,
                                 margin = margin(0,0,0,0))
bar.axis.text.x = element_text(family = "Helvetica", size = 7, colour = "black", face = "plain",
                               margin = margin(0,0,0,0))
bar.axis.text.y = element_text(family = "Helvetica", size = 7, colour = "black", face = "plain")
bar.axis.title.x = element_blank()
bar.axis.title.y = element_text(family = "Helvetica", size = 7, face = "plain")
bar.strip.text = element_text(family = "Helvetica", size = 6, colour = "black", face = "bold")

### ggplot ----
stackedbar <- ggplot(bar_df, aes(x = cohort, y = r2)) + 
  ggtitle(label = "HUMAN PLASMA: Explained variance",
          subtitle = "LPL, EL, & HL enhancement vs. ANGPTL3-4-8 suppression, using 112 overlapping NMR parameters") +
  geom_hline(yintercept = 0.25, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 0.50, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 1.00, linetype = "dotted", color = "grey75") + 
  geom_bar(stat = "identity", 
           #position = "fill",
           width = bar.binwidth, 
           color = "black", alpha = 1) + 
  #scale_fill_manual(breaks = c("y", "n"), values = c("grey25", "white")) +
  xlab("") +
  ylab("Explained [genetic] LPL activity (R\U00B2, %)") +
  #scale_y_continuous(labels = c(0, 25, 50, 75, 100)) +
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
  facet_wrap(~reg, strip.position = "bottom", nrow = 1) +
  guides(x = guide_axis(angle = 45), y = guide_axis(angle = 90))

### Save plot to file ----
ggsave(plot = stackedbar,
       path = paste0(wd, "/plots/"),
       filename = "22_crossover.pdf",
       width = 7, height = 4, units = "in", dpi = 1100)

# 22_Crossover2.R ----
## Make barplot df ----
bar_df_5 <- data.frame(
  reg = factor(
    c("a", "b", "c", "d", "e", "f"),
    labels = c("LPL (deriv.) vs.\nLPL (valid.)", 
               "EL (deriv.) vs.\nEL (valid.)",
               "HL (deriv.) vs.\nHL (valid.)",
               "ANGPTL3 (deriv.) vs.\nANGPTL3 (valid.)",
               "ANGPTL4 (deriv.) vs.\nANGPTL4 (valid.)",
               "ANGPTL8 R59W (deriv.) vs.\nANGPTL8 R59W (valid.)")
  ),
  cohort = factor(
    rep("c", 6),
    labels = "Derivation ~ Validation"
  ),
  r2 = c(
    summary(lm(formula = beta.LPL.n ~ beta.LPL.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPG.n ~ beta.LIPG.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.LIPC.n ~ beta.LIPC.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.ANGPTL3.n ~ beta.ANGPTL3.k + beta.ANGPTL8.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.ANGPTL4.n ~ beta.ANGPTL4.k, data = cbind_df_m))$r.squared,
    summary(lm(formula = beta.ANGPTL8.n ~ beta.ANGPTL8.k, data = cbind_df_m))$r.squared
  )
)

## Draw plot ----
### Graphical parameters ----
bar.width = 2.4
bar.height = 2.4
bar.binwidth = 0.7
bar.errorwidth = 0.4
bar.errorsize = 0.3
bar.plot.title = element_text(family = "Helvetica", size = 8, colour = "black", face = "bold", hjust = 0.5,
                              margin = margin(0,0,2,0))
bar.plot.subtitle = element_text(family = "Helvetica", size = 7, colour = "black", face = "italic", hjust = 0.5,
                                 margin = margin(0,0,0,0))
bar.axis.text.x = element_text(family = "Helvetica", size = 6, colour = "black", face = "bold",
                               margin = margin(0,0,0,0))
bar.axis.text.y = element_text(family = "Helvetica", size = 7, colour = "black", face = "plain")
bar.axis.title.x = element_blank()
bar.axis.title.y = element_text(family = "Helvetica", size = 7, face = "plain")
bar.strip.text = element_text(family = "Helvetica", size = 6, colour = "black", face = "bold")

### ggplot ----
stackedbar2 <- ggplot(bar_df_5, aes(x = reg, y = r2)) + 
  ggtitle(label = "HUMAN PLASMA: Explained variance",
          subtitle = "LPL, EL, & HL enhancement, and ANGPTL3-4-8 suppression") +
  geom_hline(yintercept = 0.25, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 0.50, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "grey75") + 
  geom_hline(yintercept = 1.00, linetype = "dotted", color = "grey75") + 
  geom_bar(stat = "identity", 
           #position = "fill",
           width = 0.5, 
           color = "black", alpha = 1) + 
  #scale_fill_manual(breaks = c("y", "n"), values = c("grey25", "white")) +
  xlab("") +
  ylab("Explained [genetic] LPL activity (R\U00B2, %)") +
  #scale_y_continuous(labels = c(0, 25, 50, 75, 100)) +
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
  #facet_wrap(~reg, strip.position = "bottom", nrow = 1) +
  guides(x = guide_axis(angle = 45), y = guide_axis(angle = 90))

### Save plot to file ----
ggsave(plot = stackedbar2,
       path = paste0(wd, "/plots/"),
       filename = "22_crossover2.pdf",
       width = 3, height = 4, units = "in", dpi = 1100)

# end ----