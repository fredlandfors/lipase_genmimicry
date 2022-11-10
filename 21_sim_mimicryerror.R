# Simulate genetic mimicry bias induced by small sample estimation error ----
# Clear environment
rm(list = ls())

# Get libs
library(data.table)
library(readxl)
library(MASS)
library(ggplot2)

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

## Get data ----
### Import real NMR data set ----
nmr <- readxl::read_xlsx(path = "/Users/fredriklandfors/projekt_data/2019-03-11_VIPVIZA-ITC_data/14766.xlsx",
                         sheet = "dataMatrix")
nmr <- as.data.frame(nmr)

## Import gene effects on each variable
gen <- data.table::fread(file = paste0(wd, "/data/1_Nightingale20.txt"))
gen <- as.data.frame(gen)
gen$metname <- sapply(
  gen$gwas_id.LPL,
  function(x){
    strsplit(x, split = "-")[[1]][3]
  }
)

## align data sets
gen2 <- subset(gen, metname %in% names(nmr)[-1])
nmr2 <- nmr[ ,gen2$metname]
nmr2[nmr2 == 0] <- NA
nmr2 <- na.omit(nmr2)

## Simulate LPL eQTL and ANGPTL4 E40K genotypes with known correlations from the UKBB derivation cohort ----
### Select vars ----
sim_vars <- c(names(nmr2)[-1])

### std select vars ----
nmr3 <- nmr2[sim_vars]
nmr4 <- sapply(nmr3, function(x) {scale(x)})

### Get correlation with genotypes ----
gt_effs <- subset(gen, metname %in% sim_vars)
gt_effs <- gt_effs[match(sim_vars, gt_effs$metname),]

### Calc covariance matrix ---
nmr_cov <- cov(nmr4)
nmr_cor <- cov2cor(nmr_cov)
c.lpl <- gt_effs$beta.LPL
c.a4 <- gt_effs$beta.ANGPTL4
r.lpl <- c(t(c.lpl), 1, 0)
r.a4 <- c(t(c.a4), 0, 1)
nmr_cor2 <- cbind(nmr_cor, c.lpl, c.a4)
nmr_cor2 <- rbind(nmr_cor2, r.lpl, r.a4)
nmr_cov2 <- nmr_cor2
row.names(nmr_cov2)[length(sim_vars) + 1] <- "gt_lpl"
colnames(nmr_cov2)[length(sim_vars) + 1] <- "gt_lpl"
row.names(nmr_cov2)[length(sim_vars) + 2] <- "gt_a4"
colnames(nmr_cov2)[length(sim_vars) + 2] <- "gt_a4"
nmr_cov3 <- lqmm::make.positive.definite(nmr_cov2)

### Simulate data using mvrnorm ----
nmr_sim <- MASS::mvrnorm(n = 120000, mu = rep(0, length(sim_vars) + 2), Sigma = nmr_cov3)
nmr_sim <- data.frame(nmr_sim)

### Add random noise to NMR vars
nmr_sim2 <- sapply(
  nmr_sim[!names(nmr_sim) %in% c("gt_lpl", "gt_a4")],
  function(x) {
    x + rnorm(length(x), 0, 4)
})
nmr_sim2 <- as.data.frame(nmr_sim2)
names(nmr_sim2) <- names(nmr_sim[!names(nmr_sim) %in% c("gt_lpl", "gt_a4")])
nmr_sim <- cbind(nmr_sim2, nmr_sim[names(nmr_sim) %in% c("gt_lpl", "gt_a4")])

summary(lm(formula = M_VLDL_TG ~ gt_lpl,data = nmr_sim))
summary(lm(formula = M_VLDL_TG ~ gt_a4,data = nmr_sim))

### Choose pre-specified sample sizes and run regression for each step ----
steps <- c(100, seq(500, 2500, 500), seq(5000, 30000, 2500), seq(40000, 120000, 10000))
n_steps <- length(steps)
n_nmr <- 1:nrow(nmr_sim)
gt_1 <- "gt_lpl"
gt_2 <- "gt_a4"
boots <-  1000
r2 <- matrix(0, nrow = n_steps, ncol = boots)

#### Define lm extra stats func ----
se_pval = function(object, tss) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  Qr <- z$qr
  r = z$residuals
  rss <- sum(r^2)
  r2 <- 1 - rss / tss 
  resvar <- rss/rdf
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  est <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  pval = 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  se_pval = rbind(se, pval, r2); names(se_pval) = colnames(pval)
  return(se_pval)
}

#### Run for loop ----
for (boot in seq(1, boots)) {
  for (i in seq(1, n_steps)) {
    # Subset random sample
    df_indx <- sample(n_nmr, size = steps[i])
    subset_1 <- nmr_sim[df_indx,]
    
    # Run first-stage beta estimation regressions
    ## 1
    xdatl1 <- lapply(subset_1[sim_vars], function(x){x})
    outcome1 <- subset_1[, gt_1]
    tss1 <- var(outcome1) * (length(outcome1) - 1)
    regs1 = lapply(
      names(xdatl1),
      function(id) {
        X = cbind(intercept = 1, x = xdatl1[[id]])
        fit = lm.fit(X, outcome1)
        # c(id = id, coef(fit)['x'], se_pval(fit)[, 'x']) ## Could add ID here..
        c(coef(fit)['x'], se_pval(fit, tss1)[, 'x'])
      })
    regs1 = signif(do.call("rbind", regs1), digits = 4) ## numeric matrix
    ## 2
    xdatl2 <- lapply(subset_1[sim_vars], function(x){x})
    outcome2 <- subset_1[, gt_2]
    tss2 <- var(outcome2) * (length(outcome2) - 1)
    regs2 = lapply(
      names(xdatl2),
      function(id) {
        X = cbind(intercept = 1, x = xdatl2[[id]])
        fit = lm.fit(X, outcome2)
        # c(id = id, coef(fit)['x'], se_pval(fit)[, 'x']) ## Could add ID here..
        c(coef(fit)['x'], se_pval(fit, tss2)[, 'x'])
      })
    regs2 = signif(do.call("rbind", regs2), digits = 4) ## numeric matrix
    
    # Run second-stage genetic mimicry regressions
    lmdf <- data.frame(
      y = regs1[,1],
      x = regs2[,1]
    )
    r2[i, boot] <- summary(lm(formula = y~x, data = lmdf))$r.squared
  }
  message(paste0("Boot no. ", boot, " completed on: ", Sys.time(), "."))
}

### Make mtrx into df ----
r2_2 <- data.frame(
  r2 = c(r2),
  step = rep(steps, boots)
)

### Plot ----
r2_2$step_fct <- factor(r2_2$step)
r2_2$step_fct <- reorder(r2_2$step_fct, r2_2$step)
p1 <- ggplot(data = r2_2, aes(x = step_fct, y = r2*100)) + 
  ggtitle(label = "VIOLIN DENSITY PLOT: Effect of sample size on the estimated R\U00B2 (true R\U00B2 = 0.99)"#,
          #subtitle = "Distributions of derived R\U00B2, using two variants with a true of R\U00B2 = 99 %, and 1000 random samplings per sample size step"
          ) +
  ylab("R\U00B2 (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  xlab("Sample size") +
  geom_violin(
    draw_quantiles = 0.5,
    scale = "width",
    width = 0.5
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(family = "Helvetica", size = 8, colour = "black", face = "bold", hjust = 0.5,
                              margin = margin(1,1,1,1)),
    plot.subtitle = element_text(family = "Helvetica", size = 8, colour = "black", face = "italic", hjust = 0.5,
                                 margin = margin(1,1,2,1)),
    legend.title = element_blank(),
    axis.title = element_text(family = "Helvetica", size = 8, face = "bold"),
    axis.text.x = element_text(family = "Helvetica", size = 7, colour = "black", face = "plain"),
    axis.text.y = element_text(family = "Helvetica", size = 7, colour = "black", face = "plain"),
    axis.line = element_line(size = 0.4),
    axis.ticks = element_line(size = 0.4)
  ) +
  guides(x = guide_axis(angle = 45), y = guide_axis(angle = 90))
ggsave(plot = p1,
       path = paste0(wd, "/plots/sim"),
       filename = "21_sim_mimicryerror.pdf",
       width = 5, height = 3, units = "in", dpi = 1100)

# end ----
