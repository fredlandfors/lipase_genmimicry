# Select only SNPs near gene that were GW-associated with any of the 249 parameters ----
## File index ----
files <- readxl::read_xlsx("/Volumes/LACIE_SETUP/nightingale_ukbb_gwas/nmr_gwas_ids.xlsx")

# GW significance threshold
gws_thresh <- -log10(0.05/12321875) # = 8.391707

# ANGPTL4 +-200 kb: 19:8229011-8629011 ----
## Make bash script ----
scr_a4 <- data.frame(
  v1 = sapply(
    files$gwas_id,
    function(x) {
      paste0(
        "bcftools query -f '", x, "\\t%ID\\t%CHROM\\t%POS\\t%ALT\\t%REF[\\t%ES\\t%SE\\t%LP\\t%AF]\\n' ",
        x,
        ".vcf.gz -r 19:8229011-8629011 | awk -F \"\\t\" '{ if($9 >= 8.391707) { print } }' >> gene_subset/angptl4_ukbb_subset.txt"
      )
    }
  )
)
data.table::fwrite(scr_a4, file = "./13A_select_angptl4.sh", row.names = FALSE,
                   col.names = FALSE, quote = FALSE)

# ANGPTL3 +- 200 kb: 1:62863191-63263191 ----
## Make bash script ----
scr_a3 <- data.frame(
  v1 = sapply(
    files$gwas_id,
    function(x) {
      paste0(
        "bcftools query -f '", x, "\\t%ID\\t%CHROM\\t%POS\\t%ALT\\t%REF[\\t%ES\\t%SE\\t%LP\\t%AF]\\n' ",
        x,
        ".vcf.gz -r 1:62863191-63263191 | awk -F \"\\t\" '{ if($9 >= 8.391707) { print } }' >> gene_subset/angptl3_ukbb_subset.txt"
      )
    }
  )
)
data.table::fwrite(scr_a3, file = "./13B_select_angptl3.sh", row.names = FALSE,
                   col.names = FALSE, quote = FALSE)

# ANGPTL8 +- 200 kb: 19:11150295-11552619 ----
## Make bash script ----
scr_a8 <- data.frame(
  v1 = sapply(
    files$gwas_id,
    function(x) {
      paste0(
        "bcftools query -f '", x, "\\t%ID\\t%CHROM\\t%POS\\t%ALT\\t%REF[\\t%ES\\t%SE\\t%LP\\t%AF]\\n' ",
        x,
        ".vcf.gz -r 19:11150295-11552619 | awk -F \"\\t\" '{ if($9 >= 8.391707) { print } }' >> gene_subset/angptl8_ukbb_subset.txt"
      )
    }
  )
)
data.table::fwrite(scr_a8, file = "./13C_select_angptl8.sh", row.names = FALSE,
                   col.names = FALSE, quote = FALSE)

# End ----
