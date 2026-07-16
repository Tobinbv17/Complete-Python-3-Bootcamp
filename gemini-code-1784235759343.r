# ------------------------------------------------------------------------------
# 1. Load Libraries & Setup Coordinate Space for Dnmt3a (mm10)
# ------------------------------------------------------------------------------
library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Set the genomic range around the Dnmt3a locus (mm10 coordinates)
# Dnmt3a is located on Chromosome 12: ~112,610,000 to 112,710,000 bp
chr <- "chr12"
start_coord <- 112610000
end_coord <- 112710000

# ------------------------------------------------------------------------------
# 2. Initialize Genomic Axis and Gene Annotation Tracks
# ------------------------------------------------------------------------------
# Genome axis track (Scale and coordinate indicator)
axis_track <- GenomeAxisTrack()

# Gene model track using UCSC mm10 annotations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gene_track <- GeneRegionTrack(
  txdb, 
  genome = "mm10", 
  chromosome = chr, 
  start = start_coord, 
  end = end_coord, 
  name = "Dnmt3a loci",
  transcriptAnnotation = "symbol",
  fill = "darkblue"
)

# Rename internal IDs to Gene Symbols
symbols <- mapIds(org.Mm.eg.db, keys = geneId(gene_track), column = "SYMBOL", keytype = "ENTREZID")
symbol(gene_track) <- symbols[geneId(gene_track)]

# ------------------------------------------------------------------------------
# 3. Create the FIMO STAT3 Binding Site Track
# ------------------------------------------------------------------------------
# Import the FIMO STAT3 bed file (ensure tabix/bgzip processed if necessary)
stat3_peaks_df <- read.table("GSE123488_all_STAT3.bed", sep="\t", header=FALSE)
colnames(stat3_peaks_df) <- c("chromosome", "start", "end", "motif")

# Filter for region of interest
stat3_peaks_roi <- stat3_peaks_df[
  stat3_peaks_df$chromosome == chr & 
  stat3_peaks_df$start >= start_coord & 
  stat3_peaks_df$end <= end_coord, 
]

stat3_gr <- GRanges(
  seqnames = stat3_peaks_roi$chromosome,
  ranges = IRanges(start = stat3_peaks_roi$start, end = stat3_peaks_roi$end)
)

fimo_track <- AnnotationTrack(
  stat3_gr, 
  name = "FIMO STAT\nbinding sites", 
  fill = "darkblue", 
  col = NULL
)

# ------------------------------------------------------------------------------
# 4. Construct ATAC-seq bigWig Signal Tracks
# ------------------------------------------------------------------------------
# Vector mapping the files to their respective condition labels in order
conditions <- c(
  "Ptpn2+/+ 0h"   = "Ptpn2_WT_0h.bw",
  "Ptpn2+/- 0h"   = "Ptpn2_HET_0h.bw",
  "Ptpn2+/+ 24h"  = "Ptpn2_WT_24h.bw",
  "Ptpn2+/- 24h"  = "Ptpn2_HET_24h.bw",
  "Ptpn2+/+ 48h"  = "Ptpn2_WT_48h.bw",
  "Ptpn2+/- 48h"  = "Ptpn2_HET_48h.bw",
  "Ptpn2+/+ 72h"  = "Ptpn2_WT_72h.bw",
  "Ptpn2+/- 72h"  = "Ptpn2_HET_72h.bw"
)

bw_dir <- "/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/"

# Programmatically generate tracks
data_tracks <- list()
for (i in seq_along(conditions)) {
  label <- names(conditions)[i]
  file_path <- file.path(bw_dir, conditions[i])
  
  # Set up track plotting properties matching paper visualization aesthetics
  data_tracks[[label]] <- DataTrack(
    range = file_path, 
    genome = "mm10", 
    type = "horizon",             # Polygon-filled look
    chromosome = chr, 
    name = label,
    fill.mountain = c("#F2A7A7", "#E27B7B"), # Replicates the subtle pink/violet gradient filling
    col.mountain = "#4B92DB",                # Soft blue boundary tracing
    ylim = c(0, 150)                         # Y-axis scaling to 150
  )
}

# ------------------------------------------------------------------------------
# 5. Compile and Render Panel D
# ------------------------------------------------------------------------------
# Combine all structural and genomic tracks into standard list
track_list <- c(
  axis_track,
  data_tracks[["Ptpn2+/+ 0h"]],
  data_tracks[["Ptpn2+/- 0h"]],
  data_tracks[["Ptpn2+/+ 24h"]],
  data_tracks[["Ptpn2+/- 24h"]],
  data_tracks[["Ptpn2+/+ 48h"]],
  data_tracks[["Ptpn2+/- 48h"]],
  data_tracks[["Ptpn2+/+ 72h"]],
  data_tracks[["Ptpn2+/- 72h"]],
  fimo_track,
  gene_track
)

# Plot track layout
plotTracks(
  track_list,
  from = start_coord,
  to = end_coord,
  chromosome = chr,
  sizes = c(1, rep(1.5, 8), 0.8, 1.2), # Control height proportions
  background.title = "white",
  col.title = "black",
  col.axis = "black",
  cex.title = 0.8,
  title.width = 3.5
)