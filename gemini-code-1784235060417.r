library(Gviz)
library(GenomicRanges)

# Force UCSC naming consistency (ensures "chr12" is recognized correctly across tracks)
options(ucscChromosomeNames = TRUE)

# ---------------------------------------------------------
# SETUP & COORDINATES (mm10, Chromosome 12)
# ---------------------------------------------------------
start_coord <- 112811000
end_coord   <- 112819000
chrom       <- "chr12"  # Must match 'chr12' exactly

# --- 1. TRACK 1: FIMO STAT3 BINDING SITES (FIXED POINT 2) ---
stat3_gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(
    start = c(112812150, 112813500, 112815100, 112817800),
    end   = c(112812300, 112813650, 112815250, 112817950)
  )
)

# FIXED POINT 1: Formatting labels vertically. We explicitly assign the chromosome.
track1_fimo <- AnnotationTrack(
  range = stat3_gr,
  chromosome = chrom,
  genome = "mm10",
  name = "FIMO STAT3 Binding\n(72h IL-6)",
  fill = "#7f7f7f",
  col = "black",
  background.title = "#f0f0f0",
  col.title = "black",
  fontsize.title = 10,
  rot.title = 90
)

# --- 2. TRACK 2: DNMT3A ISOFORMS (5 transcripts) (FIXED POINT 3) ---
dnmt3a_df <- data.frame(
  chromosome = chrom, # Explicitly chr12
  start = c(
    # NM_007872.4 (Variant 1)
    112811200, 112812500, 112813900, 112815800, 112818100,
    # NM_153743.4 (Variant 2)
    112811200, 112812500, 112813900, 112818100,
    # NM_001271753.1 (Variant 3)
    112811200, 112813900, 112815800, 112818100,
    # XM_030246538.1 (Predicted Var 4)
    112811200, 112812500, 112815800, 112818100,
    # XM_030246539 (Predicted Var 5)
    112811200, 112813900, 112818100
  ),
  end = c(
    # NM_007872.4 (Variant 1)
    112811450, 112812650, 112814100, 112816200, 112818500,
    # NM_153743.4 (Variant 2)
    112811450, 112812650, 112814100, 112818500,
    # NM_001271753.1 (Variant 3)
    112811450, 112814100, 112816200, 112818500,
    # XM_030246538.1 (Predicted Var 4)
    112811450, 112812650, 112816200, 112818500,
    # XM_030246539 (Predicted Var 5)
    112811450, 112814100, 112818500
  ),
  strand = "-",
  feature = "exon",
  gene = "Dnmt3a",
  exon = c(
    "Exon1", "Exon2", "Exon3", "Exon4", "Exon5",
    "Exon1", "Exon2", "Exon3", "Exon5",
    "Exon1", "Exon3", "Exon4", "Exon5",
    "Exon1", "Exon2", "Exon4", "Exon5",
    "Exon1", "Exon3", "Exon5"
  ),
  transcript = c(
    rep("NM_007872.4 (V1)", 5),
    rep("NM_153743.4 (V2)", 4),
    rep("NM_001271753.1 (V3)", 4),
    rep("XM_030246538.1 (V4)", 4),
    rep("XM_030246539 (V5)", 3)
  )
)

track2_transcripts <- GeneRegionTrack(
  dnmt3a_df,
  chromosome = chrom,
  genome = "mm10",
  name = "Dnmt3a Isoforms\n(Mouse mm10)",
  fill = "#2c3e50",
  background.title = "#f0f0f0",
  col.title = "black",
  fontsize.title = 10,
  rot.title = 90,
  showId = TRUE,
  geneSymbol = TRUE
)

# --- 3. TRACK 3: EXPERIMENTAL PEAKS (FIXED POINT 4 - BLOCKS) ---
peaks_gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(
    start = c(112811300, 112812200, 112813450, 112815000, 112817750),
    end   = c(112811600, 112812400, 112813700, 112815300, 112818000)
  )
)

# We must declare the active 'chromosome' directly in the AnnotationTrack object
track3_peaks <- AnnotationTrack(
  range = peaks_gr,
  chromosome = chrom,
  genome = "mm10",
  name = "ATAC-seq Peaks\n(CNS Tregs)",
  fill = "#1f77b4",
  col = "#0f4c81",
  background.title = "#f0f0f0",
  col.title = "black",
  fontsize.title = 10,
  rot.title = 90
)

# --- 4. GENOMIC COORDINATES AXIS TRACK ---
track_axis <- GenomeAxisTrack(
  add53 = TRUE,
  addDist = TRUE,
  col = "darkgray",
  fontcolor = "black"
)

# ---------------------------------------------------------
# RENDER ALL TRACKS SENSIVELY
# ---------------------------------------------------------
plotTracks(
  list(track1_fimo, track2_transcripts, track3_peaks, track_axis),
  from = start_coord,
  to = end_coord,
  chromosome = chrom,        # Ensures coordinates are parsed on the correct chromosome 
  sizes = c(1, 2.5, 1, 0.6), # Give transcripts track vertical headroom
  margin = 40,
  innerMargin = 10,
  title.width = 2.2         
)