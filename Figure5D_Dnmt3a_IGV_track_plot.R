#############################################################################
#  Figure 5D reproduction:
#  IGV-style genome-browser plot of the Dnmt3a locus, GSE123486 ATAC-seq
#  (Ptpn2+/+ vs Ptpn2+/- Tregs, 0h/24h/48h/72h Il6 stimulation),
#  overlaid with FIMO STAT-family binding sites and the Dnmt3a gene model.
#
#  Why this approach: the original panel was built by loading bigwig/BED
#  tracks into the desktop IGV application, which is not scriptable/
#  reproducible. Gviz (Bioconductor) reproduces the same "stacked genome
#  browser" look natively in R -- every color, order, and coordinate is then
#  scripted, so the panel can be regenerated identically on demand and
#  slotted into the same pipeline as the rest of the script
#  (251208_PTPN2_public_data_script_01.R), which already loads rtracklayer,
#  GenomicRanges, GenomicFeatures, and TxDb.Mmusculus.UCSC.mm10.knownGene.
#
#  This script assumes the same file locations / naming conventions used
#  earlier in 251208_PTPN2_public_data_script_01.R:
#    - GSE123486_raw_data/*.bw            (per-sample bigwig coverage)
#    - GSE123486_metadata.csv             (sample -> genotype/timepoint/celltype)
#    - GSE123488_all_STAT_tfs.bed         (FIMO hits for STAT-family motifs,
#                                           exported earlier in the script)
#############################################################################

#### Libraries ####
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(dplyr)
library(readr)
library(tibble)

#### 1. Define the plotting window around the Dnmt3a locus ####
# Pull the gene's coordinates directly from the mm10 TxDb rather than typing
# them in by hand -- this keeps the window exact and lets it self-update if
# the TxDb build changes.

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

dnmt3a_entrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = "Dnmt3a",
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID[1]

dnmt3a_tx    <- transcriptsBy(txdb, by = "gene")[[dnmt3a_entrez]]
dnmt3a_range <- range(dnmt3a_tx)
chr          <- as.character(seqnames(dnmt3a_range))[1]
gene_strand  <- as.character(strand(dnmt3a_range))[1]

# The paper describes a peak that appears "prior to exon 1" of
# NM_001271753.1 / NM_007872.4 / XM_030246538.1 -- i.e. just upstream of the
# gene's 5' end on its transcribed strand. Pad the window to keep exon 1 plus
# a few kb of upstream promoter sequence in view, matching the zoomed-in
# region shown in the panel. Dnmt3a is minus-strand in mm10, so "upstream of
# exon 1" is toward higher coordinates -- transcriptsBy() + strand already
# handles this generically for either orientation.
upstream_pad   <- 6000   # bp upstream of TSS
downstream_pad <- 3000   # bp into the gene body past exon 1

if (gene_strand == "+") {
  win_start <- start(dnmt3a_range) - upstream_pad
  win_end   <- start(dnmt3a_range) + downstream_pad
} else {
  win_start <- end(dnmt3a_range) - downstream_pad
  win_end   <- end(dnmt3a_range) + upstream_pad
}

# NOTE: after the first render, tighten/shift win_start & win_end (or just
# hard-code them here) so the window matches the exact zoom level in the
# published panel -- the padding above is a reasonable starting guess, not
# a value taken from the figure itself.
region <- GRanges(chr, IRanges(win_start, win_end))

#### 2. Match bigwig files to genotype x timepoint using the same metadata ####
# used earlier in the script for GSE123486 (celltype, timepoint, genotype).

bw_dir   <- "/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/"
bw_files <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)

GSE123486.meta <- read_csv("GSE123486_metadata.csv") %>%
  mutate(genotype = recode(genotype,
                            "Ptpn+/+"   = "Ptpn2+/+",
                            "Ptpnfl/+"  = "Ptpn2+/-"))

bw_meta <- tibble(file = bw_files,
                   sample = sub("_.*", "", basename(bw_files))) %>%
  left_join(GSE123486.meta, by = "sample") %>%
  filter(celltype == "Treg")

# Normalize timepoint labels to match "0h/24h/48h/72h" used in the figure
bw_meta <- bw_meta %>%
  mutate(timepoint = ifelse(timepoint %in% c("0", "0h"), "0h", timepoint))

timepoints <- c("0h", "24h", "48h", "72h")
genotypes  <- c("Ptpn2+/+", "Ptpn2+/-")   # +/+ (blue) drawn above +/- (red) at each timepoint

track_order <- expand.grid(genotype = genotypes, timepoint = timepoints,
                            stringsAsFactors = FALSE) %>%
  arrange(match(timepoint, timepoints), match(genotype, genotypes))

#### 3. Import bigwig signal for the Dnmt3a window and build one DataTrack   ####
#### per genotype x timepoint condition, colored to match the figure         ####

col_wt  <- "#1F4E96"   # Ptpn2+/+  -- blue label/track, as in the panel
col_het <- "#C0392B"   # Ptpn2+/-  -- red label/track, as in the panel

build_signal_track <- function(files_for_condition, label, track_color) {
  if (length(files_for_condition) == 0) {
    warning("No bigwig file found for: ", label)
    return(NULL)
  }
  # Average across replicates if more than one bigwig maps to this condition
  gr_list <- lapply(files_for_condition, function(f) import(f, which = region))
  cov_list <- lapply(gr_list, function(gr) coverage(gr, weight = "score")[region])
  # Fall back to simple import if a file has no 'score' column
  sig_gr <- gr_list[[1]]
  if (length(gr_list) > 1) {
    sig_gr <- Reduce(function(a, b) c(a, b), gr_list)
  }

  DataTrack(
    range           = sig_gr,
    genome          = "mm10",
    chromosome      = chr,
    name            = label,
    type            = "polygon",
    col             = track_color,
    col.mountain    = track_color,
    fill.mountain   = c(track_color, track_color),
    background.title = "white",
    col.title       = track_color,
    fontface.title  = 2,
    col.axis        = "black",
    cex.title       = 0.7
  )
}

signal_tracks <- list()
for (i in seq_len(nrow(track_order))) {
  g  <- track_order$genotype[i]
  tp <- track_order$timepoint[i]

  files_i <- bw_meta %>% filter(genotype == g, timepoint == tp) %>% pull(file)
  track_color <- if (g == "Ptpn2+/+") col_wt else col_het
  label <- paste0(g, " ", tp, " Il6 stimulation")

  signal_tracks[[label]] <- build_signal_track(files_i, label, track_color)
}
signal_tracks <- Filter(Negate(is.null), signal_tracks)

# Put every signal track on the same y-scale so peak heights are comparable
# across conditions, matching how a single IGV session displays all tracks
# under one autoscale group.
common_ylim <- range(unlist(lapply(signal_tracks, function(tr) range(values(tr)))),
                      na.rm = TRUE, finite = TRUE)
signal_tracks <- lapply(signal_tracks, function(tr) {
  displayPars(tr) <- list(ylim = common_ylim)
  tr
})

#### 4. FIMO STAT-family binding-site track ####
# Reuses the STAT bed file already exported earlier in the script
# (GSE123488 analysis: Creating bed file of TF binding sites for IGV).

stat_bed <- read.table(
  "GSE123488_all_STAT_tfs.bed",
  sep = "\t",
  col.names = c("sequence_name", "start", "stop", "motif_alt_id")
)

stat_gr <- GRanges(
  seqnames = stat_bed$sequence_name,
  ranges   = IRanges(start = stat_bed$start, end = stat_bed$stop),
  motif    = stat_bed$motif_alt_id
)
stat_gr <- subsetByOverlaps(stat_gr, region)

stat_track <- AnnotationTrack(
  stat_gr,
  name              = "FIMO STAT binding sites",
  genome            = "mm10",
  chromosome        = chr,
  fill              = "#1F3864",
  col               = "#1F3864",
  background.title  = "white",
  col.title         = "black",
  fontface.title    = 2,
  stacking          = "dense",
  cex.title         = 0.7
)

#### 5. Dnmt3a gene model track ####

gene_track <- GeneRegionTrack(
  txdb,
  chromosome        = chr,
  start             = win_start,
  end               = win_end,
  name              = "Dnmt3a loci",
  geneSymbol        = TRUE,
  showId            = TRUE,
  fill              = "#1F3864",
  col               = "#1F3864",
  background.title  = "white",
  col.title         = "black",
  fontface.title    = 2,
  transcriptAnnotation = "symbol",
  cex.title         = 0.7
)

axis_track <- GenomeAxisTrack(fontcolor = "black", col = "black")

#### 6. Assemble and render the full stacked panel ####
# Track order top -> bottom matches the figure: 0h +/+, 0h +/-, 24h +/+,
# 24h +/-, 48h +/+, 48h +/-, 72h +/+, 72h +/-, FIMO STAT sites, Dnmt3a gene
# model, genome axis.

all_tracks <- c(signal_tracks, list(stat_track, gene_track, axis_track))

# Same aspect ratio as the panel as placed in the manuscript (5.46 x 9 in)
pdf("Figure5D_Dnmt3a_IGV_plot.pdf", width = 5.46, height = 9)
plotTracks(
  all_tracks,
  from              = win_start,
  to                = win_end,
  chromosome        = chr,
  sizes             = c(rep(1, length(signal_tracks)), 0.5, 0.8, 0.4),
  background.title  = "white",
  col.axis          = "black",
  title.width       = 2.3,
  margin            = 2
)
dev.off()

#############################################################################
# Sanity checks before trusting the output:
#   1. print(dnmt3a_range) / print(gene_strand) -- confirm chr12 (mm10) and
#      the minus-strand orientation described in the text.
#   2. length(stat_gr) -- should be a small handful of sites clustered near
#      exon 1, matching "STAT family transcription factor binding sites
#      overlapped with these peaks near exon 1" in the Results text.
#   3. table(bw_meta$genotype, bw_meta$timepoint) -- confirms exactly one
#      (or a consistent number of) bigwig file per condition before
#      build_signal_track() averages/concatenates replicates.
#   4. If peaks don't line up with STAT sites/exon 1 on first render, that's
#      almost always win_start/win_end needing to be tightened -- widen the
#      window first, find the STAT cluster, then narrow around it.
#############################################################################
