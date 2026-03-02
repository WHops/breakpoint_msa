library(DECIPHER)
library(ggplot2)
library(reshape2)
library(dplyr)


trio = 5

if (trio==1){
  seqs_fasta = '/Users/whops/projects/mini-programs/crosshair/breakpoint_lab/trio1/trio1_all.fa'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/trio1.pdf'
  xplotstart = 0
  xplotend = 10000
}

if (trio==2){
  seqs_fasta = '/Users/whops/projects/mini-programs/crosshair/breakpoint_lab/trio2/trio2_all.fa'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/trio2.pdf'
  xplotstart = 44000
  xplotend = 64000
}


if (trio==3){
  seqs_fasta = '/Users/whops/projects/mini-programs/crosshair/breakpoint_lab/trio3/trio3_all.fa'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/trio3.pdf'
  xplotstart = 0
  xplotend = 20000
}

if (trio==4){
  seqs_fasta = '/Users/whops/projects/mini-programs/crosshair/breakpoint_lab/trio4/trio4_all.fa'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/trio4.pdf'
  xplotstart = 0
  xplotend = 3000
}

if (trio==6){
  
  # Trio6
  seqs_fasta = '/Users/whops/projects/mini-programs/crosshair/breakpoint_lab/trio6/trio6_all.fa'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/trio6.pdf'
  xplotstart = 15000
  xplotend = 35000

}

if (trio==66){
  seqs_fasta = '/Users/whops/projects/mini-programs/breakpoint_msa/data/T6/T6_seqs.fasta'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/T6.pdf'
  xplotstart = 0
  xplotend = 40000
}
# Trio 2


if (trio==5){
  seqs_fasta = '/Users/whops/projects/mini-programs/breakpoint_msa/data/T5/T5_seqs.fasta'
  out_file = '/Users/whops/projects/mini-programs/breakpoint_msa/plots/sep25/T5.pdf'
  xplotstart = 0
  xplotend = 22000
}
# Trio 2



reference_seqname <- "child"



# Base-specific colors for Plot 1 
base_colors <- c(
  "A" = "#7fc97f",
  "T" = "#beaed4",
  "G" = "#fdc086",
  "C" = "#ffff99",
  "-" = "#386cb0",
  "N" = "gray"
)

# Define ref_colors as a named vector (not a list)
ref_colors <- c(
  match = "orange",
  mismatch = "black",
  gap = "darkgrey"
)

# Alignment parameters (simple for now)
alignment_params <- list(
  iterations = 2,
  refinements = 2,
  gapOpening = -1000,
  gapExtension = -500,
  terminalGap = -2
)

# ===== FUNCTIONS =====

all_same <- function(column) {
  length(unique(column)) == 1
}

assign_basematch_by_reference <- function(df, ref_name, color_scheme) {
  stopifnot(all(c("Position", "Sequence", "Base") %in% names(df)))
  
  ref_bases <- df %>%
    filter(.data$Sequence == ref_name) %>%
    select(Position, Base) 
  ref_bases$Ref_base = ref_bases$Base
  ref_bases$Base = NULL
  
  # df <- left_join(df, ref_bases, by = "Position") %>%
  #   mutate(Color = case_when(
  #     Base == Ref_base ~ color_scheme$reference_color,
  #     Base == "-" & Ref_base != "-" ~ color_scheme$gap_color,
  #     TRUE ~ color_scheme$mismatch_color
  #   ))
  
  df <- left_join(df, ref_bases, by = "Position") %>%
    mutate(match = case_when(
      Base == Ref_base ~ 'match',
      Base == "-" & Ref_base != "-" ~ 'gap',
      TRUE ~ 'mismatch'
    ))

  return(df)
}

trim_leading_and_ending_gaps <- function(plot_data_ref) {
  # Get only the relevant columns
  gap_check <- plot_data_ref %>%
    group_by(Position) %>%
    summarise(all_non_gaps = all(Base != "-"))

  # Find the first and last positions where all bases are non-gap
  valid_positions <- gap_check$Position[gap_check$all_non_gaps]
  trim_start <- min(valid_positions)
  trim_end <- max(valid_positions)

  plot_data_ref_trimmed <- plot_data_ref %>%
    filter(Position >= trim_start, Position <= trim_end)
  return(plot_data_ref_trimmed)
}

BrowseSeqsCustom <- function(aligned_seqs, start, end) {
  # Ensure it's a DNAStringSet
  if (!inherits(aligned_seqs, "DNAStringSet")) {
    stop("Input must be a DNAStringSet from DECIPHER::AlignSeqs")
  }
  
  # Convert to alignment matrix
  aln_matrix <- as.matrix(aligned_seqs)
  
  # Validate range
  aln_len <- ncol(aln_matrix)
  if (start < 1 || end > aln_len || start > end) {
    stop("Invalid start/end: must be within alignment length and start <= end.")
  }
  
  # Subset matrix by alignment positions
  subset_matrix <- aln_matrix[, start:end, drop = FALSE]
  
  # Reconstruct sequences
  subset_seqs <- DNAStringSet(apply(subset_matrix, 1, paste0, collapse = ""))
  names(subset_seqs) <- rownames(subset_matrix)
  
  # Launch browser
  BrowseSeqs(subset_seqs)
}


reorder_sequences <- function(df, ref_name, seqs) {
  others <- setdiff(names(seqs), ref_name)
  df$Sequence <- factor(df$Sequence, levels = c(others[1], ref_name, others[2]))
  return(df)
}

# ===== MAIN =====

# Load sequences
seqs <- readDNAStringSet(seqs_fasta)

# Align sequences
aligned <- AlignSeqs(
  seqs,
  iterations = alignment_params$iterations,
  refinements = alignment_params$refinements,
  gapOpening = alignment_params$gapOpening,
  gapExtension = alignment_params$gapExtension,
  terminalGap = alignment_params$terminalGap,
  useStructures = FALSE
)

# Filter non-informative columns
aligned_matrix <- as.matrix(aligned)

# Build plot data frame
plot_data <- data.frame(
  Position = rep(1:ncol(aligned_matrix), times = length(seqs)),
  Sequence = rep(names(seqs), each = ncol(aligned_matrix)),
  Base = as.vector(t(aligned_matrix))
)

# ===== Prepare PLOT 1: base identity colors =====

plot_data_plot1 <- reorder_sequences(plot_data, reference_seqname, seqs)

plot_base <- ggplot(plot_data_plot1, aes(
  x = Position,
  y = Sequence,
  fill = Base
)) +
  geom_tile() +
  scale_fill_manual(values = base_colors, na.value = "white") +
  theme_minimal() +
  labs(x = "Position in alignment", y = "Sequence", fill = "Base") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# ===== Prepare PLOT 2: reference comparison colors =====

plot_data_ref <- assign_basematch_by_reference(plot_data_plot1, reference_seqname, ref_colors)
plot_data_ref <- reorder_sequences(plot_data_ref, reference_seqname, seqs)

plot_data_ref = trim_leading_and_ending_gaps(plot_data_ref)



# Plot 

# Get the factor levels in the correct plotting order
y_levels <- levels(plot_data_ref$Sequence)

# Get their numeric positions
y_pos <- seq_along(y_levels)

# Compute positions between rows for orange lines
hline_positions <- head(y_pos, -1) + 0.5

# Split data
# Group and summarize data
alignment_intervals <- plot_data_ref %>%
  group_by(Sequence, match, grp = cumsum(match != lag(match, default = dplyr::first(match)) | Sequence != lag(Sequence, default = dplyr::first(Sequence)))) %>%
  summarize(start = dplyr::first(Position), end = dplyr::last(Position), .groups = "drop")

# Extract unique sequence names
sequence_names <- levels(alignment_intervals$Sequence)

plot_ref <- ggplot() +
  geom_rect(
    data = alignment_intervals %>% 
      filter(match == "match") %>% 
      filter(!is.na(start), !is.na(end)),
    aes(
      xmin = pmin(start, end),
      xmax = pmax(start, end),
      ymin = as.numeric(Sequence),
      ymax = as.numeric(Sequence) + 1
    ),
    fill = ref_colors["match"],
    na.rm = TRUE
  ) +
  geom_rect(
    data = alignment_intervals %>% 
      filter(match == "gap") %>%
      filter(!is.na(start), !is.na(end)),
    aes(
      xmin = pmin(start, end),
      xmax = pmax(start, end),
      ymin = as.numeric(Sequence),
      ymax = as.numeric(Sequence) + 1
    ),
    fill = ref_colors["gap"],
    color = ref_colors["gap"],
    size = 0.5,
    na.rm = TRUE
  ) +
  geom_rect(
    data = alignment_intervals %>% 
      filter(match == "mismatch") %>%
      filter(!is.na(start), !is.na(end)),
    aes(
      xmin = pmin(start, end),
      xmax = pmax(start, end),
      ymin = as.numeric(Sequence),
      ymax = as.numeric(Sequence) + 1
    ),
    fill = ref_colors["mismatch"],
    color = ref_colors["mismatch"],
    size = 0.5,
    na.rm = TRUE
  ) +
  geom_hline(yintercept = c(1,2,3,4), color = "white", size = 2.5) +
  scale_y_reverse(
    breaks = seq_along(sequence_names) + 0.5,
    labels = sequence_names
  ) +
  labs(x = "MSA position", y = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(xlim = c(xplotstart, xplotend), expand = FALSE)

plot_ref

ggsave(plot_ref, file=out_file, device = pdf, width=10, height=3, units='in')

# Define the custom order using indices
custom_order_indices <- c(3, 1, 2)

custom_order_indices <- c(2, 1, 3)

# Sort the sequences based on the custom order
sorted_aligned <- aligned[custom_order_indices]
BrowseSeqs(sorted_aligned)

