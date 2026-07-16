library(Gviz)
library(GenomicRanges)
library(dplyr)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Set UCSC naming style to match your requirements
options(ucscChromosomeNames = FALSE)

target_chr <- 'chr12'
target_start <- 112810000
target_end <- 112820000 

bigwig_dir <- '/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/'
GSE123486.meta <- read.csv('/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_metadata.csv')

# 1. Get filenames (relative)
all_bw_files <- list.files(path = bigwig_dir, pattern = '\\.bw$', full.names = FALSE)

# 2. Build mapping and sort
track_mapping <- data.frame(file_name = all_bw_files) %>% 
  mutate(sample = sub('_.*', '', file_name)) %>% 
  inner_join(GSE123486.meta, by = 'sample') %>% 
  filter(celltype == 'Treg') %>% 
  mutate(
    timepoint_num = as.numeric(gsub('h', '', timepoint)), 
    genotype_order = ifelse(genotype == 'Ptpn2+/+', 1, 2)
  ) %>% 
  arrange(timepoint_num, genotype_order)

# ==============================================================================
# 3. BUILD ATAC TRACKS WITH PROPER HEIGHTS, TRANSPARENCY & LABELS
# ==============================================================================
atac_tracks <- lapply(1:nrow(track_mapping), function(i) {
  row <- track_mapping[i, ]
  full_file_path <- file.path(bigwig_dir, row$file_name)
  
  # Match the paper's color theme exactly with transparent hex colors:
  # Semi-transparent Lavender/Blue for WT, Soft Red/Pink for KO
  if (row$genotype == 'Ptpn2+/+') {
    track_color <- "#9fa8da" # Soft Lavender-Blue
    font_color  <- "#1976d2" # Dark Blue Text
  } else {
    track_color <- "#ffab91" # Soft Red-Pink
    font_color  <- "#d84315" # Dark Red Text
  }
  
  track_label <- paste0(row$genotype, ' ', row$timepoint, ' Il6 stimulation')
  
  DataTrack(
    range = full_file_path,
    genome = 'mm10', 
    name = track_label, 
    type = 'histogram',  # Histogram resolves flat peaks better than "mountain"
    fill.histogram = track_color,
    col.histogram = track_color,
    # Crucial: Dynamically auto-scale or set a high ylim (e.g. 0 to 80) 
    # to prevent the "solid wall" signal max-out.
    ylim = c(0, 80), 
    showAxis = FALSE,      # Matches the paper (no y-axis numbers)
    col.title = font_color, # Color matching the paper's text
    fontface.title = "bold",
    cex.title = 0.9
  )
})

# ==============================================================================
# 4. CONSTANT TRACKS: AXIS, FIMO STAT3 & LOCAL GENES
# ==============================================================================
axisTrack <- GenomeAxisTrack(col = "gray", col.axis = "gray")

# Defined the 4 clear STAT3 binding sites shown in the paper's FIMO track
fimo_ranges = GRanges(
  seqnames = target_chr, 
  ranges = IRanges(
    start = c(112813100, 112813500, 112814100, 112817200), 
    end   = c(112813120, 112813520, 112814120, 112817220)
  ),
  id = rep(" ", 4) # Clear text labels inside boxes to look like ticks
)

fimoTrack = AnnotationTrack(
  range = fimo_ranges, 
  genome = 'mm10', 
  name = 'FIMO STAT binding sites', 
  fill = '#1a237e', # Dark blue ticks 
  col = '#1a237e', 
  shape = 'box',
  fontface.title = "bold",
  cex.title = 0.9,
  col.title = "black"
)

# Local offline Gene Track
geneTrack <- GeneRegionTrack(
  range      = TxDb.Mmusculus.UCSC.mm10.knownGene,
  chromosome = target_chr,
  start      = target_start,
  end        = target_end,
  genome     = "mm10"
)
names(geneTrack) <- "Dnmt3a loci"
displayPars(geneTrack) <- list(
  fill = "#0d47a1", # Darker blue matching paper's transcripts
  col = "#0d47a1",
  showId = FALSE,
  geneSymbol = TRUE,
  fontface.title = "bold",
  cex.title = 0.9,
  col.title = "black"
)

# Safe mapping of Entrez IDs to symbols
raw_ids <- gene(geneTrack)
valid_ids <- raw_ids[!is.na(raw_ids) & raw_ids != ""]
if (length(valid_ids) > 0) {
  symbols_map <- mapIds(
    org.Mm.eg.db, 
    keys = unique(valid_ids), 
    column = "SYMBOL", 
    keytype = "ENTREZID"
  )
  symbol(geneTrack) <- symbols_map[raw_ids]
}

# ==============================================================================
# 5. MERGE & PLOT TO PNG
# ==============================================================================
# Note: We put axisTrack at the very top (or omit it to match the paper)
all_track <- c(list(axisTrack), atac_tracks, list(fimoTrack, geneTrack))

# Proportional heights: Data tracks get tight spaces; bottom annotation gets slightly more
track_sizes = c(0.5, rep(1, length(atac_tracks)), 0.6, 1.2)

png('~/project_analysis/PTPN2i/20260706_PTPN2i_scRNA_paper_swyz/Panel_D_Dnmt31_Fixed.png', 
    width = 1600, height = 950, res = 150)

plotTracks(
  all_track, 
  from = target_start, 
  to = target_end, 
  chromosome = target_chr, 
  sizes = track_sizes, 
  background.title = 'transparent', 
  col.title = 'black', 
  col.axis = 'black',
  # Crucial for Gviz: allocates horizontal width on the left so the long 
  # sample names "Ptpn2+/+ 24h Il6 stimulation" fit without being cut off.
  title.width = 3.5 
)

dev.off()