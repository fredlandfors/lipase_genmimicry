# Coloc analysis ----
library(data.table)
library(coloc)
library(ieugwasr)
# install.packages("devtools")
# devtools::install_github("boxiangliu/locuscomparer")

# also requires bash software wget, bgzip, bcftools and awk

# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

## Download liver expression data and annotation reference ----
# system(paste0("cd ", wd, "/data/coloc"))
# system("wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/Liver.allpairs.txt.gz")
# system("wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/all_snp_gene_associations/Whole_Blood.allpairs.txt.gz")
# system("wget https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz")

# ANGPTL3  ----

## subset 200 kb in ANGPTL3 region ----

# rs11207977 is located at 1:62,977,307 in GRCh37/hg19
# Cut a 0.2 Mb window in liver expression data
62977307 - 200000
62977307 + 200000

# system(paste0("cd ", wd, "/data/coloc"))
# system("gzcat Liver.allpairs.txt.gz | awk 'NR>1 {print $2}' | awk -F'[_/]' 'BEGIN{OFS="\t"; print "#CHROM","POS","REF","ALT", "VER"}{$1=$1}1' > addvar_Liver.allpairs.txt; bgzip -i addvar_Liver.allpairs.txt;")
# system("paste <(gzcat Liver.allpairs.txt.gz) <(gzcat addvar_Liver.allpairs.txt.gz) | awk '$10 == 1 && $11 >= 62777307 && $11 <= 63177307' > angptl3/fin_Liver.allpairs.txt; bgzip -i angptl3/fin_Liver.allpairs.txt;")
# 
# Subset annotation frame using awk
# system("gzcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz | awk '$1 == 1 && $2 >= 62777307 && $2 <= 63177307' > angptl3/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt; bgzip -i angptl3/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt;")    

## Import to R ----
chr1_liver_anno <- data.table::fread(
  file = paste0(wd, "/data/coloc/angptl3/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz")
)
data.table::setnames(
  chr1_liver_anno,
  c("chr",	"variant_pos",	"variant_id",	"ref",	"alt",	"num_alt_per_site",
    "rs_id_dbSNP147_GRCh37p13")
)
chr1_liver_expression <- data.table::fread(
  file = paste0(wd, "/data/coloc/angptl3/fin_Liver.allpairs.txt.gz")
)
data.table::setnames(
  chr1_liver_expression,
  c("gene_id",	"variant_id",	"tss_distance",	"ma_samples",	"ma_count",	"maf",
    "pval_nominal",	"slope",	"slope_se", "chr", "pos", "ref", "alt", "hg")
)
merge1 <- data.table::merge.data.table(
  x = chr1_liver_anno,
  y = chr1_liver_expression,
  by = "variant_id",
  suffixes = c(".anno", ".expr")
)

# Subset ANGPTL3 (ENSG00000132855)
merge2 <- merge1[gene_id %like% "ENSG00000132855"]

## PLOT 1: GLGC GWAS summary stats ----
### DL and write to file ----
# tg <- ieugwasr::associations(
#   variants = merge2$rs_id_dbSNP147_GRCh37p13,
#   id = "ieu-a-302"
# )
# data.table::fwrite(tg, file = paste0(wd, "/data/coloc/angptl3/glgc/tg.txt"))

# import saved data
tg <- data.table::fread(file = paste0(wd, "/data/coloc/angptl3/glgc/tg.txt"))
tg$maf <- ifelse(tg$eaf > 0.5, 1 - tg$eaf, tg$eaf)

### remove duplicate SNPs ----
duplicates <- names(table(tg$target_snp)[table(tg$target_snp) == 2])
tg <- tg[!rsid %in% duplicates]

### Merge triglyceride data with gene expression data ----
merge3 <- data.table::merge.data.table(
  x = tg,
  y = merge2,
  by.x = "target_snp",
  by.y = "rs_id_dbSNP147_GRCh37p13",
  suffixes = c(".tg", ".liver_expr")
)
merge3 <- merge3[!is.na(eaf)]

# check_ea <- merge3[, c("ea", "nea", "alt.expr", "ref.expr", "beta", "slope", "eaf", "maf")]
# check_af <-  merge3[, c("ea", "nea", "alt.anno", "ref.anno", "alt.expr", "ref.expr", "maf", "raf")]

### Run coloc analysis ----
gwas <- list(
  beta = merge3$beta,
  varbeta = merge3$se^2,
  snp = merge3$target_snp,
  position = merge3$position,
  type = "quant",
  N = 177861,
  MAF = merge3$maf.tg,
  pvalues = merge3$p
)

eqtl <- list(
  beta = merge3$slope,
  varbeta = merge3$slope_se^2,
  snp = merge3$target_snp,
  position = merge3$variant_pos,
  type = "quant",
  N = 188,
  MAF = merge3$maf.liver_expr,
  pvalues = merge3$pval_nominal
)

check_dataset(gwas)
check_dataset(eqtl)

result <- coloc.abf(
  dataset1 = gwas,
  dataset2 = eqtl
)
data.table::fwrite(
  data.frame(result$summary),
  file = paste0(wd, "/data/coloc/angptl3/glgc/coloc_summary.txt"),
  row.names = TRUE, sep = "\t"
)
data.table::fwrite(
  result$results, 
  file = paste0(wd, "/data/coloc/angptl3/glgc/coloc_results.txt"),
  row.names = FALSE, sep = "\t"
)
data.table::fwrite(
  data.frame(result$priors),
  file = paste0(wd, "/data/coloc/angptl3/glgc/coloc_priors.txt"),
  row.names = TRUE, sep = "\t"
)

### Draw locuscomparer plot and write to file ----
library(locuscomparer)
gwas_p <- data.frame(
  rsid = merge3$target_snp,
  pval = merge3$p
)
eqtl_p <- data.frame(
  rsid = merge3$target_snp,
  pval = merge3$pval_nominal
)
gwas_fn = paste0(wd, "/data/coloc/angptl3/glgc/colocplot_gwas.tsv")
eqtl_fn = paste0(wd, "/data/coloc/angptl3/glgc/colocplot_eqtl.tsv")
write.table(gwas_p, file = gwas_fn, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eqtl_p, file = eqtl_fn, sep = "\t", quote = FALSE, row.names = FALSE)
p1 <- locuscompare(
  in_fn1 = gwas_fn,
  in_fn2 = eqtl_fn,
  title1 = 'Tot-TG: GLGC (Willer) 2013 GWAS',
  title2 = 'eQTL: GTEx v7 liver ANGPTL3 expression',
  snp = "rs11207977",
  population = "EUR", 
  combine = TRUE,
  legend = TRUE,
  legend_position = c("bottomright", "topright", "topleft"), 
  lz_ylab_linebreak = FALSE, 
  genome = c("hg19")
)
ggplot2::ggsave(
  plot = p1,
  path = paste0(wd, "/plots/coloc/angptl3/glgc"),
  filename = "coloc.pdf",
  width = 14, height = 7, units = "in", dpi = 1100
)

# LPL ----

# Clear environment
rm(list = ls())

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

## subset 200 kb in ANGPTL3 region ----
# requires wget, bgzip and awk
# rs115849089 is located at 8:19,912,370 in GRCh37/hg19

# Cut a 0.2 Mb window in whole blood expression data
19912370 - 200000
19912370 + 200000

# system(paste0("cd ", wd, "/data/coloc"))
# system("gzcat Whole_Blood.allpairs.txt.gz | awk 'NR>1 {print $2}' | awk -F'[_/]' 'BEGIN{OFS="\t"; print "#CHROM","POS","REF","ALT", "VER"}{$1=$1}1' > addvar_Whole_Blood.allpairs.txt; bgzip -i addvar_Whole_Blood.allpairs.txt;")
# system("paste <(gzcat Whole_Blood.allpairs.txt.gz) <(gzcat addvar_Whole_Blood.allpairs.txt.gz) | awk '$10 == 8 && $11 >= 19712370 && $11 <= 20112370' > lpl/fin_Whole_Blood.allpairs.txt; bgzip -i lpl/fin_Whole_Blood.allpairs.txt;")
# 
# Subset annotation frame using awk
# system("gzcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz | awk '$1 == 8 && $2 >= 19712370 && $2 <= 20112370' > lpl/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt; bgzip -i lpl/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt;")

## Import to R ----
chr8_blood_anno <- data.table::fread(
  file = paste0(wd, "/data/coloc/lpl/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz")
)
data.table::setnames(
  chr8_blood_anno,
  c("chr",	"variant_pos",	"variant_id",	"ref",	"alt",	"num_alt_per_site",
    "rs_id_dbSNP147_GRCh37p13")
)
ch8_blood_expression <- data.table::fread(
  file = paste0(wd, "/data/coloc/lpl/fin_Whole_Blood.allpairs.txt.gz")
)
data.table::setnames(
  ch8_blood_expression,
  c("gene_id",	"variant_id",	"tss_distance",	"ma_samples",	"ma_count",	"maf",
    "pval_nominal",	"slope",	"slope_se", "chr", "pos", "ref", "alt", "hg")
)
merge1 <- data.table::merge.data.table(
  x = chr8_blood_anno,
  y = ch8_blood_expression,
  by = "variant_id",
  suffixes = c(".anno", ".expr")
)

# Subset LPL (ENSG00000175445)
merge2 <- merge1[gene_id %like% "ENSG00000175445"]

## PLOT 2: Kettunen GWAS summary stats ----
### DL and write to file ----
# tg <- ieugwasr::associations(
#   variants = merge2$rs_id_dbSNP147_GRCh37p13,
#   id = "met-c-934"
# )
# data.table::fwrite(tg, file = paste0(wd, "/data/coloc/lpl/kettunen/tg.txt"))

## import saved data
tg <- data.table::fread(file = paste0(wd, "/data/coloc/lpl/kettunen/tg.txt"))
tg$maf <- ifelse(tg$eaf > 0.5, 1 - tg$eaf, tg$eaf)
tg$p <- ifelse(tg$p == 0, 2.225074e-308, tg$p)

### remove duplicate SNPs ----
duplicates <- names(table(tg$target_snp)[table(tg$target_snp) == 2])
tg <- tg[!rsid %in% duplicates]

### Merge triglyceride data with gene expression data ----
merge3 <- data.table::merge.data.table(
  x = tg,
  y = merge2,
  by.x = "target_snp",
  by.y = "rs_id_dbSNP147_GRCh37p13",
  suffixes = c(".tg", ".liver_expr")
)
merge3 <- merge3[!is.na(eaf)]

# check_af <-  merge3[, c("ea", "nea", "alt.anno", "ref.anno", "alt.expr", "ref.expr", "eaf", "maf.tg")]

### Run coloc analysis ----
gwas <- list(
  beta = merge3$beta,
  varbeta = merge3$se^2,
  snp = merge3$target_snp,
  position = merge3$position,
  type = "quant",
  N = 21545,
  MAF = merge3$maf.tg,
  pvalues = merge3$p
)

eqtl <- list(
  beta = merge3$slope,
  varbeta = merge3$slope_se^2,
  snp = merge3$target_snp,
  position = merge3$variant_pos,
  type = "quant",
  N = 188,
  MAF = merge3$maf.liver_expr,
  pvalues = merge3$pval_nominal
)

check_dataset(gwas)
check_dataset(eqtl)

result <- coloc.abf(
  dataset1 = gwas,
  dataset2 = eqtl
)
data.table::fwrite(
  data.frame(result$summary),
  file = paste0(wd, "/data/coloc/lpl/kettunen/coloc_summary.txt"),
  row.names = TRUE, sep = "\t"
)
data.table::fwrite(
  result$results, 
  file = paste0(wd, "/data/coloc/lpl/kettunen/coloc_results.txt"),
  row.names = FALSE, sep = "\t"
)
data.table::fwrite(
  data.frame(result$priors),
  file = paste0(wd, "/data/coloc/lpl/kettunen/coloc_priors.txt"),
  row.names = TRUE, sep = "\t"
)

### Draw locuscomparer plot and write to file ----
library(locuscomparer)
gwas_p <- data.frame(
  rsid = merge3$target_snp,
  pval = merge3$p
)
eqtl_p <- data.frame(
  rsid = merge3$target_snp,
  pval = merge3$pval_nominal
)
gwas_fn = paste0(wd, "/data/coloc/lpl/kettunen/colocplot_gwas.tsv")
eqtl_fn = paste0(wd, "/data/coloc/lpl/kettunen/colocplot_eqtl.tsv")
write.table(gwas_p, file = gwas_fn, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eqtl_p, file = eqtl_fn, sep = "\t", quote = FALSE, row.names = FALSE)
p1 <- locuscompare(
  in_fn1 = gwas_fn,
  in_fn2 = eqtl_fn,
  snp = "rs115849089",
  title1 = 'Tot-TG: Kettunen 2016 GWAS',
  title2 = 'eQTL: GTEx v7 whole blood LPL expression',
  population = "EUR", 
  combine = TRUE,
  legend = TRUE,
  legend_position = c("bottomright", "topright", "topleft"), 
  lz_ylab_linebreak = FALSE, 
  genome = c("hg19")
)

ggplot2::ggsave(
  plot = p1,
  path = paste0(wd, "/plots/coloc/lpl/kettunen"),
  filename = "coloc.pdf",
  width = 14, height = 7, units = "in", dpi = 1100
)


