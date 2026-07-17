#### Plotting Figure 5 Panel D: IGV Track Plot for Dnmt3a Locus ####

# Load Gviz for genomic track visualization
library(Gviz)

# 1. Define the genomic region of interest (Dnmt3a locus on mm10)
# Dnmt3a is located on chromosome 12
chr <- "chr12"
start_coord <- 39000000  # Adjust specific coordinates based on your peak files
end_coord   <- 39050000

# 2. Create Genome Axis Track
axisTrack <- GenomeAxisTrack()

# 3. Create Gene Region Track using the mm10 TxDb object
geneTrack <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, 
                             chromosome = chr, 
                             start = start_coord, 
                             end = end_coord,
                             name = "Dnmt3a loci",
                             transcriptAnnotation = "symbol",
                             fill = "darkblue", 
                             col = "darkblue")

# 4. Create Annotation Track for FIMO STAT3 Binding Sites
# Using the test3 dataframe generated in your previous step
stat3_gr <- GRanges(seqnames = test3$sequence_name,
                    ranges = IRanges(start = test3$start, end = test3$stop),
                    id = test3$motif_alt_id)

stat3Track <- AnnotationTrack(stat3_gr, 
                              chromosome = chr,
                              name = "FIMO STAT binding sites",
                              fill = "darkblue", 
                              col = "darkblue",
                              shape = "box")

# 5. Create Data Tracks for BigWig files (0h, 24h, 48h, 72h for both Genotypes)
# Define your bigwig file paths dynamically based on your directory structure
bw_dir <- "/projects/users/wilsosx11/251208_PTPN2_public_datasets/GSE123486_raw_data/"

# List of files corresponding to the panels in Figure 5D
conditions <- list(
  "Ptpn2+/+ 0h Il6 stimulation"  = paste0(bw_dir, "Ptpn2_WT_0h.bw"),
  "Ptpn2+/- 0h Il6 stimulation"  = paste0(bw_dir, "Ptpn2_HET_0h.bw"),
  "Ptpn2+/+ 24h Il6 stimulation" = paste0(bw_dir, "Ptpn2_WT_24h.bw"),
  "Ptpn2+/- 24h Il6 stimulation" = paste0(bw_dir, "Ptpn2_HET_24h.bw"),
  "Ptpn2+/+ 48h Il6 stimulation" = paste0(bw_dir, "Ptpn2_WT_48h.bw"),
  "Ptpn2+/- 48h Il6 stimulation" = paste0(bw_dir, "Ptpn2_HET_48h.bw"),
  "Ptpn2+/+ 72h Il6 stimulation" = paste0(bw_dir, "Ptpn2_WT_72h.bw"),
  "Ptpn2+/- 72h Il6 stimulation" = paste0(bw_dir, "Ptpn2_HET_72h.bw")
)

# Generate DataTracks with alternating blue (WT) and red (HET) colors matching the figure
data_tracks <- list()
for (i in seq_along(conditions)) {
  track_name <- names(conditions)[i]
  bw_path <- conditions[[i]]
  
  # Determine color matching Panel D
  track_color <- ifelse(grepl("Ptpn2\\+\\/\\+", track_name), "steelblue", "darkred")
  
  data_tracks[[track_name]] <- DataTrack(range = bw_path, 
                                         chromosome = chr,
                                         type = "histogram", 
                                         name = track_name,
                                         fill.histogram = track_color, 
                                         col.histogram = track_color,
                                         baseline = 0,
                                         background.title = "white",
                                         col.axis = "black",
                                         fontcolor.title = track_color,
                                         cex.title = 0.8)
}

# 6. Plot the integrated tracks matching the stack order of Figure 5D
plotTracks(c(data_tracks, list(stat3Track, geneTrack, axisTrack)),
           from = start_coord, 
           to = end_coord,
           chromosome = chr,
           sizes = c(rep(1, 8), 0.5, 1, 0.5), # Relative track heights
           title.width = 3.5)