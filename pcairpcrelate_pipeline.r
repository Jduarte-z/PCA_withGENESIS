library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(gdsfmt)
library(BiocParallel)

# ============================================================
# CONFIG — the only place you should need to edit
# ============================================================

# ── Input: raw genotype data prefix — pipeline entry point ─────────────────
input_pfile   <- "inputPfile"
input_format  <- "pfile"   

# ── Data preparation (plink2) ──────────────────────────────────────────────
plink2_bin          <- "plink2"      
highld_regions_file <- "highLD_regions_grindeLab_hg38.tsv"  
maf_min             <- 0.01        
indep_window        <- 200           
indep_step          <- 50          
indep_r2            <- 0.2           

# Optional: exclude high-LD / long-range-LD regions before PCA.
#   "yes" -> write the regions tsv and drop those regions (extra prep step)
#   "no"  -> keep all regions; the tsv and the exclusion step are skipped
removeHighLDregions <- "yes"

# Intermediate plink2 output prefixes. Each branch uses its own names so the
# high-LD-removed and high-LD-kept variants never overwrite one another.
prep_ld <- "ld"                                                       
if (removeHighLDregions == "yes") {
  prep_noHighLD <- "outBfile_noHighLD_regions"                         
  prep_maf      <- "outBfile_noHighLD_regions_commonVars"                
  prep_pruned   <- "outBfile_noHighLD_regions_commonVars_indep-pairwise" 
} else {
  prep_maf      <- "outBfile_withHighLD_regions_commonVars"               
  prep_pruned   <- "outBfile_withHighLD_regions_commonVars_indep-pairwise" 
}
# KING uses the raw data as bed. If the input is already bed, use it directly;
# if it's a pfile, it is converted to this prefix in prep step 5.
king_bed      <- if (input_format == "bfile") input_pfile else "inputPfile_bed"

# ── PLINK prefixes consumed by the R analysis (derived from prep outputs) ──
pca_plink_prefix  <- prep_pruned     
king_plink_prefix <- king_bed        

# ── Outputs: GDS files ─────────────────────────────────────────────────────
pca_gds  <- "onlyTyped4PCA.gds"
king_gds <- "onlyTyped4King.gds"

# ── Outputs: RDS files ─────────────────────────────────────────────────────
king_mat_rds     <- "KINGmat_rds"                      
pcair_r1_rds     <- "mypcair_1round_results_rds_1stRound" 
pcrel_r1_rds     <- "mycprelate_1stround"                 
pcrel_r1_mat_rds <- "mypcrel_1round_mat_rds"              
pcair_r2_rds     <- "mypcair_r2_results_rds"              
pcrel_r2_rds     <- "mypcrel_r2_rds"                      
#pending to make the pcrel_r2_mat_rds that is in kinship scale 2, for the h2 project

# ── Analysis parameters ────────────────────────────────────────────────────
n_cores           <- 2      
n_pcs             <- 10     
snp_block_size    <- 10000  
sample_block_size <- 7000   

# ── Output filename patterns (sprintf style) ───────────────────────────────
# PC-pair scatter plots: %02d, %02d = the two PC numbers on the axes
pcpair_png_pattern_r1 <- "plot_PC%02d_PC%02d.png"
pcpair_png_pattern_r2 <- "plot_PC%02d_PC%02d_2ndround.png"
# SNP–PC correlation plots/tables: %02d = PC number
snpcorr_png_pattern_r1 <- "snpcorr_1stRound_PC%02d.png"
snpcorr_tsv_pattern_r1 <- "snpcorr_1stRound_PC%02d.tsv"
snpcorr_png_pattern_r2 <- "snpcorr_2ndRound_PC%02d.png"
snpcorr_tsv_pattern_r2 <- "snpcorr_2ndRound_PC%02d.tsv"

# ── Plot dimensions (pixels / dpi) ─────────────────────────────────────────
pcpair_plot_width  <- 1500;  pcpair_plot_height  <- 1500; pcpair_plot_res  <- 300
snpcorr_plot_width <- 1800;  snpcorr_plot_height <- 700;  snpcorr_plot_res <- 150

# ============================================================
# HELPERS — generic, no hardcoded values below this block
# ============================================================

# Convert a PLINK bed/bim/fam triple (given by prefix) to a GDS file,
# then open it, print a summary, and close it again.
bed2gds <- function(prefix, out_gds) {
  snpgdsBED2GDS(
    bed.fn    = paste0(prefix, ".bed"),
    bim.fn    = paste0(prefix, ".bim"),
    fam.fn    = paste0(prefix, ".fam"),
    out.gdsfn = out_gds
  )
  gds <- snpgdsOpen(out_gds)
  snpgdsSummary(out_gds)
  snpgdsClose(gds)
}

# Does a PLINK bed/bim/fam triple already exist for this prefix?
bed_exists <- function(prefix) all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))

# Run one plink2 command; stop the pipeline if it exits non-zero.
run_plink <- function(args) {
  args <- as.character(args)
  message("plink2 ", paste(args, collapse = " "))
  status <- system2(plink2_bin, args)
  if (status != 0L)
    stop("plink2 failed (exit ", status, "): ", paste(args, collapse = " "))
}
#endregion

# ============================================================
# DATA PREPARATION — plink2 (formerly processGeneticDataForPCA.sh)
# ============================================================
# Builds the two bed datasets the analysis consumes:
#   pca_plink_prefix  (optionally high-LD removed, maf-filtered, LD-pruned) -> PCA
#   king_plink_prefix (raw pfile converted to bed)                          -> KING
# High-LD removal is controlled by removeHighLDregions in CONFIG.
# Skipped automatically when both already exist, so re-runs are cheap.

# #region - R : plink2 data prep
if (bed_exists(pca_plink_prefix) && bed_exists(king_plink_prefix)) {
  message("Prep outputs already present — skipping plink2 data preparation.")
} else {
  # Read the raw data as a pfile or an already-existing bfile
  input_flag <- if (input_format == "bfile") "--bfile" else "--pfile"

  # 1) Optionally drop high-LD / long-range-LD regions from the raw data -> bed.
  #    When removeHighLDregions = "no", this whole step (and the tsv) is skipped
  #    and the maf filter reads straight from the raw input instead.
  if (removeHighLDregions == "yes") {
    # High-LD regions to exclude (Grinde lab, hg38).
    # Columns: chrom  start  end  label   (tab-separated, no header)
    highld_regions <- c(
      "1\t47761741\t51822307\tanderson1_price1_michigan1",
      "2\t129125957\t139525961\ttopmedLCT_michigan3_priveceliac1_price3_privepopres1_raskalct",
      "2\t182309767\t189427029\tmichigan4_anderson3_price4",
      "3\t47483506\t49987563\tanderson4_michigan5_price5",
      "3\t83368159\t86868160\tanderson5_michigan6_price6",
      "3\t161899518\t163699518\tpriveceliac4",
      "5\t98636396\t101136397\tmichigan9_price9",
      "5\t129636408\t132636409\tanderson7_michigan10_price10",
      "5\t136136412\t139136412\tmichigan11_price11",
      "6\t23691793\t38924246\tpriveceliac2_topmedMHC_raskahla_fellay2_anderson8_michigan12_price12_privepopres2",
      "6\t139637170\t142137170\tanderson10_michigan14_price14",
      "8\t6455071\t13598120\tpriveceliac3_privepopres3_topmedinversion_anderson12_fellay3_michigan16_price16_raskainv",
      "8\t110918595\t113918595\tanderson14_michigan18_price18",
      "11\t88127184\t91127184\tanderson16_michigan21_price21",
      "12\t110577812\t113099475\tprice23_michigan23",
      "14\t47061047\t47961047\tpriveceliac5",
      "17\t42394456\t46567318\ttopmedinversion",
      "20\t33948533\t36438183\tanderson18_michigan24_price24"
    )
    writeLines(highld_regions, highld_regions_file)

    run_plink(c(input_flag, input_pfile,
                "--exclude", "bed1", highld_regions_file,
                "--make-bed", "--out", prep_noHighLD))
    maf_input <- c("--bfile", prep_noHighLD)   
  } else {
     # maf filter reads the raw input directly in case high-LD regions are kept
    message("removeHighLDregions = 'no' — keeping high-LD regions; tsv and exclusion step skipped.")
    maf_input <- c(input_flag, input_pfile)   
  }

  # 2) Minor-allele-frequency filter
  run_plink(c(maf_input,
              "--maf", maf_min,
              "--make-bed", "--out", prep_maf))

  # 3) LD pruning: build the list of variants to prune
  run_plink(c("--bfile", prep_maf,
              "--indep-pairwise", indep_window, indep_step, indep_r2,
              "--out", prep_ld))

  # 4) Apply the prune list -> LD-pruned bed (pca_plink_prefix)
  run_plink(c("--bfile", prep_maf,
              "--exclude", paste0(prep_ld, ".prune.out"),
              "--make-bed", "--out", prep_pruned))

  # 5) Convert the raw pfile to bed for KING / SNP-correlations (king_plink_prefix).
  #    Skipped when the input is already bed — king_bed then points at it directly.
  if (input_format != "bfile") {
    run_plink(c("--pfile", input_pfile,
                "--make-bed", "--out", king_bed))
  }
}
# ============================================================
# DATA PREPARATION — convert plink files to gds files that R can read
# ============================================================
bed2gds(pca_plink_prefix,  pca_gds)
bed2gds(king_plink_prefix, king_gds)

# ============================================================
# Run KING — under their robust method, the kinship coefficient is computed to identify the set of unrelated individuals
# ============================================================
# KING is time-consuming: reuse a matrix from a previous run if one exists,
# otherwise compute it and cache it to king_mat_rds.
if (file.exists(king_mat_rds)) {
  message("Reusing existing KING matrix: ", king_mat_rds)
  KINGmat <- readRDS(king_mat_rds)
} else {
  message("Computing KING matrix -> ", king_mat_rds)
  gds_king <- snpgdsOpen(king_gds)

  king <- snpgdsIBDKING(gds_king, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,type=c("KING-robust"), family.id=NULL, verbose=TRUE)

  KINGmat <- king$kinship
  rownames(KINGmat) <- colnames(KINGmat) <- king$sample.id
  snpgdsClose(gds_king)

  saveRDS(KINGmat, king_mat_rds)
}

# KING matrix (square numeric matrix)
dim(KINGmat)                     # n_samples x n_samples
KINGmat[1:5, 1:5]               # top-left corner



# ============================================================
# Run the first Round of PC-AiR using the input from KING
# ============================================================
#load the gds file
geno_reader <- GdsGenotypeReader(filename = pca_gds)
geno_data <- GenotypeData(geno_reader)

#load the KING matrix
KINGmat <- readRDS(king_mat_rds)

#run pcair
mypcair_1round <- pcair(geno_data, kinobj = KINGmat, divobj = KINGmat, num.cores=n_cores, eigen.cnt = n_pcs)
saveRDS(mypcair_1round, pcair_r1_rds)
close(geno_data)
summary(mypcair_1round)

#plot the pcs
for (i in seq(1, n_pcs - 1, by = 2)) {
  png(sprintf(pcpair_png_pattern_r1, i, i + 1),
      width = pcpair_plot_width, height = pcpair_plot_height, res = pcpair_plot_res)
  plot(mypcair_1round, vx = i, vy = i + 1)
  dev.off()
}


# ============================================================
# Run the first Round of correlations of genotypes with PCs
# ============================================================
#consider that is the full dataset, no ld exclusion, no ld pruninng, no maf filtering, nothing of that


genofile <- snpgdsOpen(king_gds)

chr    <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos    <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

cr1 <- snpgdsPCACorr(
  pcaobj     = mypcair_1round$vectors[, 1:n_pcs],
  gdsobj     = genofile,
  eig.which  = 1:n_pcs,
  num.thread = n_cores
)

snpgdsClose(genofile)

# align annotation to exactly the SNPs snpgdsPCACorr used
snp_idx_1round     <- match(cr1$snp.id, snp_id)
chr_used_1round    <- chr[snp_idx_1round]
pos_used_1round   <- pos[snp_idx_1round]
snp_id_used_1round <- snp_id[snp_idx_1round]

# chromosome x-axis: midpoint SNP index per chromosome, labeled chr1, chr2...
chr_levels_1round <- sort(unique(chr_used_1round))
chr_mids_1round   <- tapply(seq_along(chr_used_1round), chr_used_1round, median)
chr_labels_1round <- paste0("chr", names(chr_mids_1round))

# alternating colors by chromosome
# one color per chromosome
chr_colors        <- palette()[((as.integer(chr_levels_1round) - 1) %% length(palette())) + 1]
names(chr_colors) <- chr_levels_1round
col_vec           <- chr_colors[as.character(chr_used_1round)]
for (i in 1:n_pcs) {
  abs_corr_1round <- abs(cr1$snpcorr[i, ])

  # PNG with chromosome x-axis
  png(sprintf(snpcorr_png_pattern_r1, i),
      width = snpcorr_plot_width, height = snpcorr_plot_height, res = snpcorr_plot_res)
  plot(abs_corr_1round,
       col  = col_vec,
       pch  = 20,
       cex  = 0.3,
       ylim = c(0, 1),
       xaxt = "n",
       xlab = "Chromosome",
       ylab = paste0("PC", i, " Correlation"),
       main = paste0("PC-AiR 1st round | PC", i, " — SNP correlation"))
  axis(1, at = chr_mids_1round, labels = chr_labels_1round, las = 2, cex.axis = 0.7)
  dev.off()
  # TSV sorted by absolute correlation descending
  df_1round <- data.frame(
    snp_id     = snp_id_used_1round,
    chromosome = chr_used_1round,
    position   = pos_used_1round,
    abs_corr_1round   = abs_corr_1round
  )
  df_1round <- df_1round[order(df_1round$abs_corr_1round, decreasing = TRUE), ]
  write.table(df_1round, sprintf(snpcorr_tsv_pattern_r1, i),
              sep = "\t", row.names = FALSE, quote = FALSE)
}



# ============================================================
# Run the first Round of PC-Relate using the input from the first round of PC-Air
# ============================================================
mypcair_1round <- readRDS(pcair_r1_rds)
geno_reader <- GdsGenotypeReader(filename = pca_gds)
geno_data   <- GenotypeData(geno_reader)
geno_iter   <- GenotypeBlockIterator(geno_data, snpBlock = snp_block_size)

mypcrel_1round <- pcrelate(geno_iter, pcs = mypcair_1round$vectors[, 1:n_pcs],
                            sample.block.size = sample_block_size,
                            training.set = mypcair_1round$unrels,
                            BPPARAM = BiocParallel::MulticoreParam(n_cores))
saveRDS(mypcrel_1round, pcrel_r1_rds)
close(geno_data)

# ── Convert for round 2 PC-AiR input: scaleKin=1 ──────────────────────────
# scaleKin=1 keeps the kinship-coefficient scale expected by pcair's kinobj
mypcrel_1round_mat <- pcrelateToMatrix(mypcrel_1round,
                                       scaleKin = 1,
                                       thresh=NULL,
                                       verbose = TRUE)
saveRDS(mypcrel_1round_mat, pcrel_r1_mat_rds)
# ============================================================
# Run the second round of PC-AiR using the kinships from PC-Relate round 1 that are more accurate than those from KING. 
#however, the divergence signal is still required to come from KING
# ============================================================
# kinobj  → PC-Relate r1 (accurate recent relatedness)
# divobj  → KING (retains divergence signal across ancestry groups)
KINGmat <- readRDS(king_mat_rds)
geno_reader  <- GdsGenotypeReader(filename = pca_gds)
geno_data    <- GenotypeData(geno_reader)

mypcair_r2 <- pcair(
  geno_data,
  kinobj    = mypcrel_1round_mat,
  divobj    = KINGmat,
  num.cores = n_cores,
  eigen.cnt = n_pcs
)
saveRDS(mypcair_r2, pcair_r2_rds)
close(geno_data)
summary(mypcair_r2)

###plot the new PCs
for (i in seq(1, n_pcs - 1, by = 2)) {
  png(sprintf(pcpair_png_pattern_r2, i, i + 1),
      width = pcpair_plot_width, height = pcpair_plot_height, res = pcpair_plot_res)
  plot(mypcair_r2, vx = i, vy = i + 1)
  dev.off()
}
# ============================================================
# Run the secound Round of correlations of genotypes with PCs
# ============================================================

genofile <- snpgdsOpen(king_gds)

chr    <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos    <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

cr2 <- snpgdsPCACorr(
  pcaobj     = mypcair_r2$vectors[, 1:n_pcs],
  gdsobj     = genofile,
  eig.which  = 1:n_pcs,
  num.thread = n_cores
)

snpgdsClose(genofile)

# align annotation to exactly the SNPs snpgdsPCACorr used
snp_idx_2r     <- match(cr2$snp.id, snp_id)
chr_used_2r    <- chr[snp_idx_2r]
pos_used_2r   <- pos[snp_idx_2r]
snp_id_used_2r <- snp_id[snp_idx_2r]

# chromosome x-axis: midpoint SNP index per chromosome, labeled chr1, chr2...
chr_levels_2r <- sort(unique(chr_used_2r))
chr_mids_2r   <- tapply(seq_along(chr_used_2r), chr_used_2r, median)
chr_labels_2r <- paste0("chr", names(chr_mids_2r))

# alternating colors by chromosome
# one color per chromosome
chr_colors        <- palette()[((as.integer(chr_levels_2r) - 1) %% length(palette())) + 1]
names(chr_colors) <- chr_levels_2r
col_vec           <- chr_colors[as.character(chr_used_2r)]
for (i in 1:n_pcs) {
  abs_corr_2r <- abs(cr2$snpcorr[i, ])
  # PNG with chromosome x-axis
  png(sprintf(snpcorr_png_pattern_r2, i),
      width = snpcorr_plot_width, height = snpcorr_plot_height, res = snpcorr_plot_res)
  plot(abs_corr_2r,
       col  = col_vec,
       pch  = 20,
       cex  = 0.3,
       ylim = c(0, 1),
       xaxt = "n",
       xlab = "Chromosome",
       ylab = paste0("PC", i, " Correlation"),
       main = paste0("PC-AiR 2nd round | PC", i, " — SNP correlation"))
  axis(1, at = chr_mids_2r, labels = chr_labels_2r, las = 2, cex.axis = 0.7)
  dev.off()

  # TSV sorted by absolute correlation descending
  df_2r <- data.frame(
    snp_id     = snp_id_used_2r,
    chromosome = chr_used_2r,
    position   = pos_used_2r,
    abs_corr_2r   = abs_corr_2r
  )
  df_2r <- df_2r[order(df_2r$abs_corr_2r, decreasing = TRUE), ]
  write.table(df_2r, sprintf(snpcorr_tsv_pattern_r2, i),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ============================================================
# Run the second round of PC-Relate, this is for the h2 project
# ============================================================
geno_reader <- GdsGenotypeReader(filename = pca_gds)
geno_data   <- GenotypeData(geno_reader)
geno_iter   <- GenotypeBlockIterator(geno_data, snpBlock = snp_block_size)
mypcrel_r2 <- pcrelate(
  geno_iter,
  pcs          = mypcair_r2$vectors[, 1:n_pcs],
  training.set = mypcair_r2$unrels,
  BPPARAM      = BiocParallel::MulticoreParam(n_cores)
)
saveRDS(mypcrel_r2, pcrel_r2_rds)
close(geno_data)
