library(DECIPHER)
library(ggplot2)
library(reshape2)
library(dplyr)


# Function to check if all elements in a vector are the same
all_same <- function(column) {
  #return(column['CH'] == column['BP2'])
  return(length(unique(column)) == 1)
}

color_list_for_plot <- function(plot_data, colors){
  
  plot_data$Color = colors['mismatch_color']
  
  # Get a helper column indicating what the 'ref' sequence is per position
  ref_bases = plot_data[plot_data$Sequence == reference_seqname,c('Position','Base')]
  colnames(ref_bases)[2] = 'Ref_base'
  plot_data = left_join(plot_data, ref_bases, by="Position")
  
  # Per position: if base is the same base as the reference, then color = ref color. 
  plot_data[plot_data$Base == plot_data$Ref_base, 'Color'] <- colors['reference_color']
  
  # Per position: if ref base is not '-', but the seq base is '-', then color = grey
  plot_data$Color[plot_data$Base == '-' & plot_data$Ref_base != '-'] <- colors['gap_color']
  
  return(plot_data)
}


fill_empty_spaces_for_secondplot <- function(plot_data){
  # Step 1 & 2: Create a complete range of positions
  complete_positions <- seq(min(plot_data$Position), max(plot_data$Position))
  
  # Step 3: Expand your data frame
  expanded_data <- expand.grid(Position = complete_positions, Sequence = unique(plot_data$Sequence))
  
  # Step 4: Merge with the original data
  merged_data <- merge(expanded_data, plot_data, by = c("Position", "Sequence"), all.x = TRUE)
  
  # Step 5: Fill missing colors with orange
  merged_data$Color[is.na(merged_data$Color)] <- 'orange'
  
  return(merged_data)
}

plot_data_factorize_put_ref_in_middle <- function(plot_data, seqs){
  plot_data$Sequence <- factor(plot_data$Sequence, 
                               levels = c(names(seqs)[names(seqs) != reference_seqname][1],
                                          reference_seqname, 
                                          names(seqs)[names(seqs) != reference_seqname][2]))
  return(plot_data)
}

# Load the sequences from the file
seqs_fasta <- "../data/trio4/all.fa"

seqs_fasta = '/Users/Z364220/projects/mini-programs/breakpoint_msa/data/toy.fa'
seqs_fasta = '/Users/Z364220/projects/15q13/trio3/breakpoint_lab/15kbp/trio3_bps_all.fa'

seqs <- readDNAStringSet(seqs_fasta)

aln_custom = T
gapOpen = -10
gapExtension = -10000000
misMatch = -1
terminalGap = -0
iterations = 10
refinements = 10
x_label_dist <- 10  # for example, every 10th position


reference_seqname = 'child'
colors = list(
  'reference_color' = "orange",
  'mismatch_color' = 'darkred',
  'gap_color' = 'grey'
)


if (aln_custom){
  aligned <- AlignSeqs(
            seqs, 
            iterations=iterations, 
            refinements=refinements, 
            gapOpening = c(gapOpen,gapOpen), 
            gapExtension = c(gapExtension,gapExtension),
            useStructures = F,
            terminalGap = terminalGap,
            misMatch = misMatch
            )
} else {
  aligned <- AlignSeqs(seqs)

}




# Convert alignment to a matrix for easier manipulation
aligned_matrix <- as.matrix(aligned)


# Apply the function to each column to get a logical vector where TRUE means the column should be removed
columns_to_remove <- apply(aligned_matrix, 2, all_same)
positions_to_keep <- which(columns_to_remove==F)
filtered_aligned_matrix <- aligned_matrix[, positions_to_keep]


# Get the original positions of the kept columns
original_positions_filtered <- which(!columns_to_remove)


# Create a data frame with the original positions, sequence identifiers, and bases
plot_data <- data.frame(
  Position = rep(original_positions_filtered, times = 3),
  Sequence = factor(rep(names(seqs), each = length(original_positions_filtered))),
  Base = as.vector(t(filtered_aligned_matrix))
)

plot_data = color_list_for_plot(plot_data, colors)

### PLOT ###

# Sort and prepare plot_data
plot_data = plot_data[order(plot_data$Position, plot_data$Sequence),]
plot_data = plot_data_factorize_put_ref_in_middle(plot_data,seqs)

# Create vectors for x-axis breaks and labels (every nth position)
x_breaks <- unique(ceiling(order(plot_data$Position)/3))[seq(1, length(unique(plot_data$Position)), by = x_label_dist)]
x_labels <- unique(plot_data$Position)[seq(1, length(unique(plot_data$Position)), by = x_label_dist)]

plot_data[plot_data$Base == '-', 'Color'] = 'grey'

# Plot using geom_raster
plot_data2 = plot_data[plot_data$Position > start & plot_data$Position < end,]
plot_compress = ggplot(plot_data2, aes(x = ceiling(order(plot_data2$Position)/3), 
                                      y = Sequence, 
                                      fill = as.numeric(as.factor(Base)))) +
    geom_raster() +
    scale_fill_identity() +
    theme_minimal() +
    labs(x = "Position in original alignment", y = "Sequence") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
    guides(fill = FALSE) +  # Remove the legend for fill
    scale_x_continuous(breaks = x_breaks, labels = x_labels)  

plot_compress
  
plot_data_full = fill_empty_spaces_for_secondplot(plot_data)

start = 0
end = 100000
plot_data_full2 = plot_data_full[plot_data_full$Position > start & plot_data_full$Position < end,]

plot_full <- ggplot(plot_data_full2[plot_data_full2$Color != 'orange',],
                    aes(x = Position, y = Sequence, fill = Color)) +
  geom_tile(width=10) +
  scale_fill_identity() +
  theme_minimal() +
  labs(x = "Position in original alignment", y = "Sequence") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = FALSE)  
plot_full
plot_full <- ggplot(plot_data_full[plot_data_full$Color != 'orange',], 
                    aes(x = Position, y = Sequence, fill = Color)) +
  geom_tile(width=10) +
  scale_fill_identity() +
  theme_minimal() +
  labs(x = "Position in original alignment", y = "Sequence") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill = FALSE)  
# Remove the legend for fill
#plot_full = plot_full + 
#  ggplot(plot_data_full[plot_data_full$Color != 'orange',], 
#         aes(x = Position, y = Sequence, fill = Color)) + geom_raster()

#print(plot_compress)
print(plot_full)

