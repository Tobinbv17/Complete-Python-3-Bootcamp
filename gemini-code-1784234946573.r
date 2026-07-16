library(Gviz)
library(GenomicRanges)

# Disable UCSC chromosome name strictness checks if necessary
options(ucscChromosomeNames = FALSE)

# ---------------------------------------------------------
# SETUP & COORDINATES (mm10, Chromosome 12)
# ---------------------------------------------------------
start_coord <- 112811000
end_coord   <- 112819000
chrom       <- "chr12"

# --- 1. TRACK 1: FIMO STAT3 BINDING SITES (FIXED POINT 2) ---
# Coordinates mapping across the promoter and intronic enhancers
stat3_gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(
    start = c(112812150, 112813500, 112815100, 112817800),
    end   = c(112812300, 112813650, 112815250, 112817950)
  ),
  id = c("Site 1", "Site 2", "Site 3", "Site 4")
)

# FIXED POINT 1: Set rot.title = 90 (or 0 if horizontal is preferred)
# Gviz handles the clean, single-line alignment automatically.
track1_fimo <- AnnotationTrack(
  range = stat3_gr,
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
# Constructing all 5 key transcript variants for Mouse Dnmt3a manually
dnmt3a_df <- data.frame(
  chromosome = chrom,
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

# Convert metadata to track structure
track2_transcripts <- GeneRegionTrack(
  dnmt3a_df,
  genome = "mm10",
  name = "Dnmt3a Isoforms\n(Mouse mm10)",
  fill = "#2c3e50",
  background.title = "#f0f0f0",
  col.title = "black",
  fontsize.title = 10,
  rot.title = 90,
  showId = TRUE,          # Displays the transcript IDs next to each model
  geneSymbol = TRUE
)

# --- 3. TRACK 3: EXPERIMENTAL PEAKS (FIXED POINT 4) ---
# Representing regions as discrete BED blocks rather than continuous lines
peaks_gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(
    start = c(112811300, 112812200, 112813450, 112815000, 112817750),
    end   = c(112811600, 112812400, 112813700, 112815300, 112818000)
  )
)

track3_peaks <- AnnotationTrack(
  range = peaks_gr,
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
# Set up a plotting canvas that fits all labels without overlapping margins
plotTracks(
  list(track1_fimo, track2_transcripts, track3_peaks, track_axis),
  from = start_coord,
  to = end_coord,
  chromosome = chrom,
  sizes = c(1, 2.5, 1, 0.6), # Give transcripts track the most vertical room
  margin = 40,
  innerMargin = 10,
  title.width = 2.2         # Widens the label area so text doesn't break/wrap ugly
)