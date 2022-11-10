# LD analysis ----
library(LDlinkR)
library(data.table)
library(coloc)
# also requires bash software wget, bgzip and awk

# Set working directory
wd <- "/Users/fredriklandfors/projekt/lipase_genmimicry"

# Access token from https://ldlink.nci.nih.gov/?tab=apiaccess are required. ----
# access_token = # Get access token via instructions from https://ldlink.nci.nih.gov
pop = "EUR"

# Get LD R2  for  variant using LDlinkR -----
## ANGPTL8: rs2278426 -----
# angptl8 <- LDlinkR::LDproxy(
#   snp = "rs2278426",
#   pop = pop,
#   r2d = "r2",
#   token = access_token,
#   file = read.delim(paste0(wd, "/data/coloc/angptl8_ld/angptl8.txt"))
# )

# Read saved file
angptl8 <- data.table::fread(file = paste0(wd, "/data/coloc/angptl8_ld/angptl8.txt"))
angptl8$chr <- 19
angptl8$pos <- sapply(angptl8$Coord, function(x) {
    strsplit(x, ":")[[1]][2]
  })
angptl8$Alleles2 <- sapply(angptl8$Alleles, function(x) {
   y <- gsub("\\(", "", x)
   z <- gsub("\\)", "", y)
   return(z)
  })
angptl8$ref <- sapply(angptl8$Alleles2, function(x) {
  strsplit(x, "/")[[1]][1]
  })
angptl8$alt <- sapply(angptl8$Alleles2, function(x) {
    strsplit(x, "/")[[1]][2]
  })
angptl8$variant_id <- paste0(angptl8$chr, "_", angptl8$pos, "_",
                             angptl8$ref, "_", angptl8$alt, "_b37")

# Import GTEx v7 chromosome 19 EUR expression data ----

## Get data and preprocess in bash ----
# requires wget, bgzip and awk
# rs2278426 is located at 19:11350488 in GRCh37/hg19

# Cut a 1 Mb window liver expression data
# 11350488 - 1000000
# 11350488 + 1000000
# system("gzcat Liver.allpairs.txt.gz | awk 'NR>1 {print $2}' | awk -F'[_/]' 'BEGIN{OFS="\t"; print "#CHROM","POS","REF","ALT", "VER"}{$1=$1}1' > addvar_Liver.allpairs.txt; bgzip -i addvar_Liver.allpairs.txt;")
# system("paste <(gzcat Liver.allpairs.txt.gz) <(gzcat addvar_Liver.allpairs.txt.gz) | awk '$10 == 19 && $11 >= 10350488 && $11 <= 12350488' > angptl8_ld/fin_Liver.allpairs.txt; bgzip -i angptl8_ld/fin_Liver.allpairs.txt;")
# 
# Subset annotation frame using awk
# system("gzcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz | awk '$1 == 19 && $2 >= 10350488 && $2 <= 12350488' > angptl8_ld/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt; bgzip -i angptl8_ld/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt;")

## Import to R ----
chr19_liver_anno <- data.table::fread(
  file = paste0(wd, "/data/coloc/angptl8_ld/fin_GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz")
)
data.table::setnames(
  chr19_liver_anno,
  c("chr",	"variant_pos",	"variant_id",	"ref",	"alt",	"num_alt_per_site",
    "rs_id_dbSNP147_GRCh37p13")
)
chr19_liver_expression <- data.table::fread(
  file = paste0(wd, "/data/coloc/angptl8_ld/fin_Liver.allpairs.txt.gz")
)
data.table::setnames(
  chr19_liver_expression,
  c("gene_id",	"variant_id",	"tss_distance",	"ma_samples",	"ma_count",	"maf",
    "pval_nominal",	"slope",	"slope_se", "chr", "pos", "ref", "alt", "hg")
)
merge1 <- data.table::merge.data.table(
  x = chr19_liver_anno,
  y = chr19_liver_expression,
  by = "variant_id",
  suffixes = c(".anno", ".expr")
)

# Subset LDLR (ENSG00000130164)
merge2 <- merge1[gene_id %like% "ENSG00000130164"]

# Merge linkage data with gene expression data
merge3 <- data.table::merge.data.table(
  x = angptl8,
  y = merge2,
  by = "variant_id",
  suffixes = c(".ldlink", ".expr")
)

# Subset R2 > 0.1
merge4 <- merge3[R2 >= 0.1]
merge4$pval_adjusted <- p.adjust(merge4$pval_nominal, method = "bonferroni")

# Save to file
data.table::fwrite(merge4, file = paste0(wd, "/data/coloc/angptl8_ld/data_out.txt"))
  