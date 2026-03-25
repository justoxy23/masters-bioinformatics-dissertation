# ============================
# 📂 1. Mutation Analysis Setup
# ============================
# ===========================================
# 📂 2. Compare Mutation Density: Full Genome vs PA1874 Region
# ===========================================
# Define genomes and sliding-window parameters
genome_ids <- c("A005", "A011", "A012", "A033", "A037", "A041", "A086", "A091", "A097", "A110")
genome_length <- 6385154
window_size <- 500
step_size <- 100
window_starts <- seq(0, genome_length - window_size, by = step_size)

# Define PA1874 region of interest
zoom_start <- 3566000
zoom_end <- 3574000

# Create output folder for comparison plots
if (!dir.exists("plots/zoom_vs_full")) dir.create("plots/zoom_vs_full", recursive = TRUE)

# Generate full-genome and zoomed-in mutation density plots for each genome
for (id in genome_ids) {
  message("Plotting full + zoom for: ", id)
  
  ## Load breseq mutation calls
  gd <- read.table(paste0("breseq_", id, ".gd"), header = FALSE, sep = "\t", comment.char = "#", fill = TRUE, quote = "")
  positions_breseq <- as.numeric(gd[gd$V1 %in% c("SNP", "DEL", "INS", "MOB", "SUB", "AMP", "UNAMP", "CNV"), 5])
  breseq_counts <- sapply(window_starts, function(start) {
    sum(positions_breseq >= start & positions_breseq < start + window_size)
  })
  
  ## Load marginal mutation calls
  html_lines <- readLines(paste0("marginal_", id, ".html"), warn = FALSE)
  positions_raw <- regmatches(html_lines, gregexpr("<td align=\"right\">[0-9,]+</td>", html_lines))
  positions_clean <- gsub("<.*?>", "", unlist(positions_raw))
  positions_marginal <- unique(as.numeric(gsub(",", "", positions_clean))[seq(1, length(positions_clean), 2)])
  marginal_counts <- sapply(window_starts, function(start) {
    sum(positions_marginal >= start & positions_marginal < start + window_size)
  })
  
  ## Load snippy mutation calls
  snippy <- read.csv(paste0("snippy_", id, ".csv"))
  positions_snippy <- snippy$POS
  snippy_counts <- sapply(window_starts, function(start) {
    sum(positions_snippy >= start & positions_snippy < start + window_size)
  })
  
  ## Prepare full-genome and zoomed-in plots
  in_zoom <- window_starts >= zoom_start & window_starts <= zoom_end
  y_max <- max(breseq_counts, marginal_counts, snippy_counts, na.rm = TRUE)
  y_max_zoom <- max(breseq_counts[in_zoom], marginal_counts[in_zoom], snippy_counts[in_zoom], na.rm = TRUE)
  
  png(filename = paste0("plots/zoom_vs_full/mutation_zoom_vs_full_", id, ".png"), width = 1200, height = 500)
  par(mfrow = c(1, 2))
  
  # Plot mutation density across the full genome
  plot(window_starts, breseq_counts, type = "l", col = "steelblue", lwd = 2,
       ylim = c(0, y_max), main = paste("Full Genome -", id),
       xlab = "Genome Position", ylab = "Mutation Count")
  lines(window_starts, marginal_counts, col = "firebrick", lwd = 2)
  lines(window_starts, snippy_counts, col = "darkgreen", lwd = 2)
  
  # Plot mutation density within the PA1874 region
  plot(window_starts[in_zoom], breseq_counts[in_zoom], type = "l", col = "steelblue", lwd = 2,
       ylim = c(0, y_max_zoom), main = paste("Zoom around PA1874 -", id),
       xlab = "Genome Position", ylab = "Mutation Count")
  lines(window_starts[in_zoom], marginal_counts[in_zoom], col = "firebrick", lwd = 2)
  lines(window_starts[in_zoom], snippy_counts[in_zoom], col = "darkgreen", lwd = 2)
  
  legend("topright", legend = c("Breseq", "Marginal", "Snippy"),
         col = c("steelblue", "firebrick", "darkgreen"), lty = 1, lwd = 2, bty = "n")
  
  dev.off()
}

# ===========================================
# 📂 3. Explore Breseq-Only Mutation Density Around PA1874
# ===========================================

# Define genome-wide sliding-window parameters
genome_length <- 6385154
window_size <- 500
step_size <- 100
window_starts <- seq(0, genome_length - window_size, by = step_size)

# Define zoomed region around PA1874
zoom_start <- 3560000
zoom_end   <- 3600000
in_zoom <- window_starts >= zoom_start & window_starts <= zoom_end

# Load breseq mutation calls for genome A005
gd <- read.table("breseq_A005.gd", header = FALSE, sep = "\t", comment.char = "#", fill = TRUE, quote = "")
positions_breseq <- as.numeric(gd[gd$V1 %in% c("SNP", "DEL", "INS", "MOB", "SUB", "AMP", "UNAMP", "CNV"), 5])

# Count breseq mutations within each sliding window
breseq_counts <- sapply(window_starts, function(start) {
  sum(positions_breseq >= start & positions_breseq < start + window_size)
})

# Plot zoomed-in breseq mutation density around PA1874
y_max_zoom <- max(breseq_counts[in_zoom], na.rm = TRUE)

plot(
  window_starts[in_zoom], breseq_counts[in_zoom], type = "l", col = "steelblue", lwd = 2,
  ylim = c(0, y_max_zoom),
  main = "Breseq Mutation Density Zoomed In (PA1874 ±500bp) - A005",
  xlab = "Genome Position", ylab = "Mutation Count"
)

# ===========================================
# 📂 4. Convert Marginal SNP Positions to Relative PA1874 Coordinates
# ===========================================

# Define genome IDs to process
genome_ids <- c("A005", "A011", "A012", "A025","A026","A027","A028A067","A032","A033", "A037", "A041", "A042","A043","A045","A047","A066","A078","A086", "A091", "A097", "A110", "A120", "A122", "A128", "A130", "A131", "A132", "A133", "A138", "A147", "A156", "A158", "A169", "A185", "A214", "A281", "A304", "A510")
# Define genomic coordinates for PA1874
gene_start <- 3566483
gene_end <- 3572893

# Create output folder for relative position tables
if (!dir.exists("marginal_relative_tables")) dir.create("marginal_relative_tables")

# Process marginal output files and extract SNP positions within PA1874
for (id in genome_ids) {
  cat("Processing genome:", id, "\n")
  
  # Define expected marginal HTML file path
  html_path <- file.path("marginal_htmls", paste0("marginal_", id, "_.*\\.html$"))
  
  if (!file.exists(html_path)) {
    warning("Missing file: ", html_path, " — skipping.")
    next
  }
  
  html_lines <- readLines(html_path, warn = FALSE)
  
  # Extract numeric mutation positions from HTML table entries
  matches <- regmatches(html_lines, gregexpr("<td align=\"right\">[0-9,]+</td>", html_lines))
  raw_positions <- unlist(matches)
  clean_positions <- gsub("<.*?>", "", raw_positions)
  numeric_positions <- as.numeric(gsub(",", "", clean_positions))
  
  # Keep only positions that fall within PA1874
  positions <- numeric_positions[seq(1, length(numeric_positions), by = 2)]
  
  # Filter to positions in gene range
  positions_in_gene <- positions[positions >= gene_start & positions <= gene_end]
  
  # Convert genomic positions to relative amino-acid positions
  # PA1874 lies on the reverse strand, so coordinates are reversed
  relative_positions <- (gene_end - positions_in_gene + 1)/3 # bc of reverse strand
  
  # Save table
  df <- data.frame(
    Genome = id,
    Genomic_Position = positions_in_gene,
    Relative_Position = relative_positions
  )
  
  # Write to CSV
  write.csv(df, file = paste0("marginal_relative_tables/marginal_", id, "_relative.csv"), row.names = FALSE)
}

# ===========================================
# 📂 5. Extract Relative SNP Positions from Multiple Marginal Files Per Genome
# ===========================================

# Define input directory containing marginal HTML files

input_dir <- "marginal_htmls"

for (id in genome_ids) {
  cat("Processing genome:", id, "\n")
  
  # Find all files for this genome ID
  files <- list.files(input_dir, pattern = paste0("marginal_", id, "_.*\\.html$"), full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No files found for genome: ", id)
    next
  }
  
  for (html_path in files) {
    cat("  Reading file:", basename(html_path), "\n")
    
    html_lines <- readLines(html_path, warn = FALSE)
    matches <- regmatches(html_lines, gregexpr("<td align=\"right\">[0-9,]+</td>", html_lines))
    raw_positions <- unlist(matches)
    clean_positions <- gsub("<.*?>", "", raw_positions)
    numeric_positions <- as.numeric(gsub(",", "", clean_positions))
    
    positions <- numeric_positions[seq(1, length(numeric_positions), by = 2)]
    positions_in_gene <- positions[positions >= gene_start & positions <= gene_end]
    
    relative_positions <- ceiling((gene_end - positions_in_gene + 1) / 3)
    
    df <- data.frame(
      Genome = id,
      File = basename(html_path),
      Genomic_Position = positions_in_gene,
      Relative_Position = relative_positions
    )
    
    out_file <- paste0("marginal_relative_tables/", tools::file_path_sans_ext(basename(html_path)), "_relative.csv")
    write.csv(df, file = out_file, row.names = FALSE)
  }
}

# ===========================================
# 📂 6. Combine Relative-to-AlphaFold Position Mappings
# ===========================================
library(dplyr)

# Define folder containing original mapping tables
og_path <- "marginal_relative_tables/og"

# Identify all mapping files
csv_files <- list.files(og_path, pattern = "\\.csv$", full.names = TRUE)

# Merge all unique relative-to-AlphaFold mappings into one table
all_mappings <- lapply(csv_files, function(file) {
  read.csv(file) %>%
    select(Relative_Position, Alphafold_Position) %>%  # adjust if column names differ
    distinct()
}) %>%
  bind_rows() %>%
  distinct() %>%
  arrange(Relative_Position)

# Export combined mapping table
write.csv(all_mappings, "combined_relative_to_alphafold_mapping.csv", row.names = FALSE)

# ===========================================
# 📂 7. Map Relative SNP Positions to AlphaFold Coordinates
# ===========================================
library(dplyr)

# Define input, mapping, and output paths
mapping_file <- "combined_relative_to_alphafold_mapping.csv"
marginal_folder <- "marginal_relative_tables"   
output_folder <- "marginal_with_alphafold"

# Load combined relative-to-AlphaFold mapping
mapping <- read.csv(mapping_file)

# Create output folder if not exists
if (!dir.exists(output_folder)) dir.create(output_folder)

# Identify marginal relative-position files to process
marginal_files <- list.files(marginal_folder, pattern = "\\.csv$", full.names = TRUE)

# Rerun loop but store warnings manually
warn_list <- c()

# Join AlphaFold coordinates onto each marginal SNP table
for (file in marginal_files) {
  
  df <- read.csv(file)
  
  # Append AlphaFold positions using relative-position mapping
  df_new <- df %>%
    left_join(mapping, by = "Relative_Position")
  
  # Check for unmapped positions
  missing_positions <- df_new %>%
    filter(is.na(Alphafold_Position)) %>%
    pull(Relative_Position) %>%
    unique()
  
  if (length(missing_positions) > 0) {
    warn_list <- c(warn_list, paste0(basename(file), ": ", paste(missing_positions, collapse = ", ")))
  }
  
  # Save updated file
  write.csv(df_new, file.path(output_folder, basename(file)), row.names = FALSE)
}

cat("✅ Processing complete. Updated files saved in:", output_folder, "\n")

# View all warnings collected
warn_list

# ===========================================
# 📂 8. Summarise Breseq vs Snippy SNP Calls Across Genomes
# ===========================================
library(dplyr)
library(readr)
library(ggplot2)

# Define genomes included in the comparison
genome_ids <- c("A005","A011","A012","A033","A037","A041","A086","A091","A097","A110")

# Combine per-genome comparison tables into one dataset
all_points <- bind_rows(lapply(genome_ids, function(id) {
  df <- read.csv(paste0("breseq_vs_snippy_", id, ".csv"), na.strings = "")
  df$Genome <- id
  df
}))

# Shared SNPs by Breseq and Snippy
shared <- all_points %>% filter(!is.na(Breseq) & !is.na(Snippy))
# Breseq-only SNPs
breseq_only <- all_points %>% filter(!is.na(Breseq) &  is.na(Snippy))
# Snippy-only SNPs
snippy_only <- all_points %>% filter(is.na(Breseq) & !is.na(Snippy))

# Visualise agreement and disagreement between callers across all genomes
gg <- ggplot() +
  geom_point(data = shared, aes(x = Breseq, y = Snippy), 
             colour = "#e59866", alpha = 0.6, size = 2) +
  geom_point(data = breseq_only, aes(x = Breseq, y = min(shared$Snippy, na.rm = TRUE) - 10000), 
             colour = "blue", shape = 17, alpha = 0.6) +
  geom_point(data = snippy_only, aes(x = min(shared$Breseq, na.rm = TRUE) - 10000, y = Snippy), 
             colour = "darkgreen", shape = 17, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(
    title = "Breseq vs Snippy SNP Calls (All Genomes Combined)",
    x = "Breseq Position",
    y = "Snippy Position"
  ) +
  theme_minimal()

#ggsave("breseq_vs_snippy_all_combined.png", gg, width = 8, height = 6, dpi = 300)
gg

# ===========================================
# 📂 9. Create Base R Summary Plot for Breseq vs Snippy Calls
# ===========================================
library(dplyr)
library(readr)

genome_ids <- c("A005","A011","A012","A033","A037","A041","A086","A091","A097","A110")

# Combine all CSVs into one big df
all_points <- bind_rows(lapply(genome_ids, function(id) {
  df <- read.csv(paste0("breseq_vs_snippy_", id, ".csv"), na.strings = "")
  df
}))

# Prepare subsets
shared <- all_points[!is.na(all_points$Breseq) & !is.na(all_points$Snippy), ]
breseq_only <- all_points[!is.na(all_points$Breseq) & is.na(all_points$Snippy), ]
snippy_only <- all_points[is.na(all_points$Breseq) & !is.na(all_points$Snippy), ]

# Axis limits
x_min <- min(c(all_points$Breseq, all_points$Snippy), na.rm = TRUE)
x_max <- max(c(all_points$Breseq, all_points$Snippy), na.rm = TRUE)

# Create scatter plot comparing breseq and snippy SNP positions
plot(NA, NA,
     xlim = c(x_min, x_max),
     ylim = c(x_min, x_max),
     xlab = "SNP Position Called by Breseq on PES Genome (bp)",
     ylab = "SNP Position Called by Snippy on PES Genome (bp)",
     main = "Breseq vs Snippy Mutation Calls (All Genomes)")

# Add shared and caller-specific SNP calls
points(shared$Breseq, shared$Snippy, pch = 19, col = "#e59866")        # Shared
points(breseq_only$Breseq, rep(x_min, nrow(breseq_only)), pch = 17, col = "blue")  # Breseq only
points(rep(x_min, nrow(snippy_only)), snippy_only$Snippy, pch = 17, col = "darkgreen")  # Snippy only

# Add y = x reference line
abline(a = 0, b = 1, col = "red", lty = 2)

# Add legend summarising overlap counts
legend(x = x_max * 0.6, y = x_max * 0.3,
       legend = c(paste0("Shared (", nrow(shared), ")"),
                  paste0("Breseq Only (", nrow(breseq_only), ")"),
                  paste0("Snippy Only (", nrow(snippy_only), ")")),
       col = c("#e59866", "blue", "darkgreen"),
       pch = c(19, 17, 17),
       bty = "n")

# ===========================================
# 📂 10. Generate Per-Genome Breseq vs Snippy Comparison Plots
# ===========================================
# List of genome IDs
genome_ids <- c("A005", "A011", "A012", "A033", "A037", "A041", "A086", "A091", "A097", "A110")

# Create output folder for plots
if (!dir.exists("plots/breseq_vs_snippy")) dir.create("plots/breseq_vs_snippy", recursive = TRUE)

for (id in genome_ids) {
  # Load SNP comparison table for the current genome
  filename <- paste0("breseq_vs_snippy_", id, ".csv")
  df <- read.csv(filename, na.strings = "")
  
  # Prepare subsets
  shared <- df[!is.na(df$Breseq) & !is.na(df$Snippy), ]
  breseq_only <- df[!is.na(df$Breseq) & is.na(df$Snippy), ]
  snippy_only <- df[is.na(df$Breseq) & !is.na(df$Snippy), ]
  
  # Axis limits
  x_min <- min(c(df$Breseq, df$Snippy), na.rm = TRUE) 
  x_max <- max(c(df$Breseq, df$Snippy), na.rm = TRUE) 
  
  # Open PNG device to export plot
  png(filename = paste0("plots/breseq_vs_snippy/", id, "_breseq_vs_snippy.png"), width = 1000, height = 800)
  
  # Create plot
  plot(NA, NA,
       xlim = c(x_min, x_max),
       ylim = c(x_min, x_max),
       xlab = "Breseq Position",
       ylab = "Snippy Position",
       main = paste("Breseq vs Snippy Mutation Calls -", id))
  
  points(shared$Breseq, shared$Snippy, pch = 19, col = "#e59866")  # Shared
  points(breseq_only$Breseq, rep(x_min, nrow(breseq_only)), pch = 17, col = "blue")  # Breseq only
  points(rep(x_min, nrow(snippy_only)), snippy_only$Snippy, pch = 17, col = "darkgreen")  # Snippy only
  
  # Diagonal
  abline(a = 0, b = 1, col = "red", lty = 2)
  
  # Add legend with counts
  legend(x = x_max * 0.6, y = x_max * 0.3,
         legend = c(paste0("Shared (", nrow(shared), ")"),
                    paste0("Breseq Only (", nrow(breseq_only), ")"),
                    paste0("Snippy Only (", nrow(snippy_only), ")")),
         col = c("#e59866", "blue", "darkgreen"),
         pch = c(19, 17, 17),
         bty = "n")
  
  dev.off()
}

# ===========================================
# 📂 11. Compare Full Genome vs Zoomed Region Across Selected Years
# ===========================================
# Setup
genome_id <- "A005"
years <- c("2002","2008")
genome_length <- 6385154
window_size <- 500
step_size <- 100
window_starts <- seq(0, genome_length - window_size, by = step_size)

# Region of interest
zoom_start <- 3523000
zoom_end <- 3532000

# Create output folder
if (!dir.exists("plots/zoom_vs_full")) dir.create("plots/zoom_vs_full", recursive = TRUE)

# Generate comparison plots for each selected year
for (year in years) {
  message("Plotting full + zoom for: ", genome_id, " - ", year)
  
  ## Load Breseq
  gd_path <- file.path("A005", paste0("breseq_", genome_id, "_", year, ".gd"))
  gd <- read.table(gd_path, header = FALSE, sep = "\t", comment.char = "#", fill = TRUE, quote = "")
  positions_breseq <- as.numeric(gd[gd$V1 %in% c("SNP", "DEL", "INS", "MOB", "SUB", "AMP", "UNAMP", "CNV"), 5])
  breseq_counts <- sapply(window_starts, function(start) {
    sum(positions_breseq >= start & positions_breseq < start + window_size)
  })
  
  ## Load Marginal
  html_lines <- readLines(file.path("A005", paste0("marginal_", genome_id, "_", year, ".html")), warn = FALSE)
  positions_raw <- regmatches(html_lines, gregexpr("<td align=\"right\">[0-9,]+</td>", html_lines))
  positions_clean <- gsub("<.*?>", "", unlist(positions_raw))
  positions_marginal <- unique(as.numeric(gsub(",", "", positions_clean))[seq(1, length(positions_clean), 2)])
  marginal_counts <- sapply(window_starts, function(start) {
    sum(positions_marginal >= start & positions_marginal < start + window_size)
  })
  
  ## Load Snippy
  html_lines_snippy <- readLines(paste0("A005/snps_", genome_id, "_", year, ".html"), warn = FALSE)
  
  # Extract lines that match a table row with a numeric POS value
  # Look for all TD entries and parse the second one (which is POS)
  td_lines <- grep("<TD>", html_lines_snippy, value = TRUE)
  td_matrix <- matrix(td_lines, ncol = 6, byrow = TRUE)  # assumes at least 6 fields per variant
  positions_snippy <- as.numeric(gsub("<.*?>", "", td_matrix[, 2]))  # extract POS column
  positions_snippy <- positions_snippy[!is.na(positions_snippy)]  # remove NAs
  
  snippy_counts <- sapply(window_starts, function(start) {
    sum(positions_snippy >= start & positions_snippy < start + window_size)
  })
  
  ## Plotting
  in_zoom <- window_starts >= zoom_start & window_starts <= zoom_end
  # y_max <- max(breseq_counts, marginal_counts, snippy_counts, na.rm = TRUE)
  # y_max_zoom <- max(breseq_counts[in_zoom], marginal_counts[in_zoom], snippy_counts[in_zoom], na.rm = TRUE)
  y_max <- 15
  y_max_zoom <- 15
  
  png(filename = paste0("plots/zoom_vs_full/mutation_zoom_vs_full_", genome_id, "_", year, ".png"), width = 1200, height = 500)
  par(mfrow = c(1, 2))
  
  # Full genome
  plot(window_starts, breseq_counts, type = "l", col = "steelblue", lwd = 2,
       ylim = c(0, y_max), main = paste("Full Genome -", genome_id, year),
       xlab = "Genome Position", ylab = "Mutation Count")
  lines(window_starts, marginal_counts, col = "firebrick", lwd = 2)
  lines(window_starts, snippy_counts, col = "darkgreen", lwd = 2)
  
  # Zoom
  plot(window_starts[in_zoom], breseq_counts[in_zoom], type = "l", col = "steelblue", lwd = 2,
       ylim = c(0, y_max_zoom), main = paste("Zoom around PA1874 -", genome_id, year),
       xlab = "Genome Position", ylab = "Mutation Count")
  lines(window_starts[in_zoom], marginal_counts[in_zoom], col = "firebrick", lwd = 2)
  lines(window_starts[in_zoom], snippy_counts[in_zoom], col = "darkgreen", lwd = 2)
  
  legend("topright", legend = c("Breseq", "Marginal", "Snippy"),
         col = c("steelblue", "firebrick", "darkgreen"), lty = 1, lwd = 2, bty = "n")
  
  dev.off()
}

# ===========================================
# 📂 12. Compare Full Genome vs Zoomed Region for Isolates with P-Style IDs
# ===========================================
# Setup
genome_id <- "A138"
years <- c("P2375_2004")
# "P2189_2004", "P2376_2004", "P2377_2004", "P2378_2004"
genome_length <- 6385154
window_size <- 500
step_size <- 100
window_starts <- seq(0, genome_length - window_size, by = step_size)

# Region of interest
zoom_start <- 3566000
zoom_end <- 3574000

# Create output folder
if (!dir.exists("plots/zoom_vs_full")) dir.create("plots/zoom_vs_full", recursive = TRUE)

# Loop over all A012 genomes
for (year in years) {
  message("Plotting full + zoom for: ", genome_id, " - ", year)
  
  ## ---- Load Breseq ----
  gd_path <- file.path("A138", paste0("breseq_", genome_id, "_", year, ".gd"))
  gd <- read.table(gd_path, header = FALSE, sep = "\t", comment.char = "#", fill = TRUE, quote = "")
  positions_breseq <- as.numeric(gd[gd$V1 %in% c("SNP", "DEL", "INS", "MOB", "SUB", "AMP", "UNAMP", "CNV"), 5])
  breseq_counts <- sapply(window_starts, function(start) {
    sum(positions_breseq >= start & positions_breseq < start + window_size)
  })
  
  ## ---- Load Marginal ----
  html_lines <- readLines(file.path("A138", paste0("marginal_", genome_id, "_", year, ".html")), warn = FALSE)
  positions_raw <- regmatches(html_lines, gregexpr("<td align=\"right\">[0-9,]+</td>", html_lines))
  positions_clean <- gsub("<.*?>", "", unlist(positions_raw))
  positions_marginal <- unique(as.numeric(gsub(",", "", positions_clean))[seq(1, length(positions_clean), 2)])
  marginal_counts <- sapply(window_starts, function(start) {
    sum(positions_marginal >= start & positions_marginal < start + window_size)
  })
  
  ## ---- Load Snippy ----
  html_lines_snippy <- readLines(file.path("A138", paste0("snps_", genome_id, "_", year, ".html")), warn = FALSE)
  td_lines <- grep("<TD>", html_lines_snippy, value = TRUE)
  td_matrix <- matrix(td_lines, ncol = 6, byrow = TRUE)
  positions_snippy <- as.numeric(gsub("<.*?>", "", td_matrix[, 2]))
  positions_snippy <- positions_snippy[!is.na(positions_snippy)]
  
  snippy_counts <- sapply(window_starts, function(start) {
    sum(positions_snippy >= start & positions_snippy < start + window_size)
  })
  
  ## ---- Plotting ----
  in_zoom <- window_starts >= zoom_start & window_starts <= zoom_end
  y_max <- 15
  y_max_zoom <- 15
  
  png(filename = paste0("plots/zoom_vs_full/mutation_zoom_vs_full_", genome_id, "_", year, ".png"), width = 1200, height = 500)
  par(mfrow = c(1, 2))
  
  # Full genome
  plot(window_starts, breseq_counts, type = "l", col = "steelblue", lwd = 2,
       ylim = c(0, y_max), main = paste("Full Genome -", genome_id, year),
       xlab = "Genome Position", ylab = "Mutation Count")
  lines(window_starts, marginal_counts, col = "firebrick", lwd = 2)
  lines(window_starts, snippy_counts, col = "darkgreen", lwd = 2)
  
  # Zoom
  plot(window_starts[in_zoom], breseq_counts[in_zoom], type = "l", col = "steelblue", lwd = 2,
       ylim = c(0, y_max_zoom), main = paste("Zoom around PA1874 -", genome_id, year),
       xlab = "Genome Position", ylab = "Mutation Count")
  lines(window_starts[in_zoom], marginal_counts[in_zoom], col = "firebrick", lwd = 2)
  lines(window_starts[in_zoom], snippy_counts[in_zoom], col = "darkgreen", lwd = 2)
  
  legend("topright", legend = c("Breseq", "Marginal", "Snippy"),
         col = c("steelblue", "firebrick", "darkgreen"), lty = 1, lwd = 2, bty = "n")
  
  dev.off()
}

# ===========================================
# 📂 13. Count Frequency of AlphaFold-Mapped Positions Across Files
# ===========================================
# Set working directory to where the CSVs are stored
setwd("marginal_relative_tables/")

# Get all CSV file names
files <- list.files(pattern = "*.csv")

# Initialize a vector to collect all Alphafold positions
all_positions <- c()

# Loop through files and extract Alphafold_Position column
for (file in files) {
  data <- read.csv(file)
  
  # If the column has any naming issues, adjust here:
  if ("Alphafold_Position" %in% colnames(data)) {
    all_positions <- c(all_positions, data$Alphafold_Position)
  } else {
    warning(paste("No Alphafold_Position column in", file))
  }
}

# Count frequency of each position
position_counts <- table(all_positions)

# Convert to a sorted dataframe
position_df <- as.data.frame(position_counts)
colnames(position_df) <- c("Alphafold_Position", "Count")
position_df$Alphafold_Position <- as.integer(as.character(position_df$Alphafold_Position))
position_df <- position_df[order(-position_df$Count), ]

# View result
print(position_df)

# Optional: save to CSV
write.csv(position_df, "../alphafold_position_counts.csv", row.names = FALSE)

# Set working directory to where the CSVs are stored
setwd("marginal_with_alphafold/")

# Get all CSV file names
files <- list.files(pattern = "*.csv")

# Initialize a list to store unique positions per file
all_positions <- c()

# Loop through files and extract unique Alphafold_Position values
for (file in files) {
  data <- read.csv(file)
  
  # Check the column exists
  if ("Alphafold_Position" %in% colnames(data)) {
    unique_positions <- unique(data$Alphafold_Position)
    all_positions <- c(all_positions, unique_positions)
  } else {
    warning(paste("No Alphafold_Position column in", file))
  }
}

# Count in how many files each position appears
position_counts <- table(all_positions)

# Convert to a sorted dataframe
position_df <- as.data.frame(position_counts)
colnames(position_df) <- c("Alphafold_Position", "File_Count")
position_df$Alphafold_Position <- as.integer(as.character(position_df$Alphafold_Position))
position_df <- position_df[order(-position_df$File_Count), ]

# View result
print(position_df)

# Optional: save to CSV
write.csv(position_df, "../alphafold_position_counts.csv", row.names = FALSE)

# ===========================================
# 📂 14. Final Docking Scatter Plot with Cluster and Energy Thresholds
# ===========================================
library(dplyr)
library(ggplot2)
library(ggrepel)


#  Load data
docking_data <- data.frame(
  Interaction = c(
    "PA1875–PA1874", "PA1876–PA1874", "PA1877–PA1874",
    "4844ctpL–5369pstS", "3477rhIR–1000pqsE", "1713exsA–1714exsD", "1706pcrV–1708popB", "1706pcrV–1709popD",
    "2360T6SS–PA1874", "4964parC–PA1874", "4270rpoB–PA1874", "3477rhIR–5369pstS", "4844ctpL–1000pqsE",
    "0085T6SS–1713exsA", "0090clpv1–1706pcrV", "1714exsD–4525pilA"
  ),
  Best_Energy_Score = c(
    -1477.5, -1706.3, -1320.6,
    -1176.2, -1041.1, -1198.6, -1140.7, -1198.6,
    -1192.0, -1005.8, -1102.3, -929.2, -1142.7,
    -898.5, -989.5, -982.5
  ),
  Largest_Cluster_Size = c(
    77, 26, 31,
    56, 91, 100, 49, 60,
    44, 44, 36, 92, 53,
    70, 50, 110
  ),
  Group = c(
    "PA187X-PA1874", "PA187X-PA1874", "PA187X-PA1874",
    "Positive control", "Positive control", "Positive control", "Positive control", "Positive control",
    "Negative control", "Negative control", "Negative control", "Negative control", "Negative control",
    "Negative control", "Negative control", "Negative control"
  ),
  stringsAsFactors = FALSE
)

# Remove "0085T6SS–1713exsA" AND ambiguous controls
docking_data <- docking_data %>%
  filter(
    !(Interaction %in% c("0085T6SS–1713exsA",
                         "1714exsD–4525pilA",
                         "3477rhIR–5369pstS"))
  )

# Color/shape maps (only 3 groups now)
group_colors <- c(
  "PA187X-PA1874" = "#2F6BFF",
  "Positive control" = "seagreen3",
  "Negative control" = "grey40"
)
shape_map <- c(
  "PA187X-PA1874" = 18,   # diamond
  "Positive control" = 17,# triangle
  "Negative control" = 16 # circle
)

# Compute thresholds using CLEAN negatives
neg_clean <- docking_data %>% filter(Group == "Negative control")

p75_neg <- quantile(neg_clean$Largest_Cluster_Size, 0.75, na.rm = TRUE)
p95_neg <- quantile(neg_clean$Largest_Cluster_Size, 0.95, na.rm = TRUE)
mad_cut <- median(neg_clean$Largest_Cluster_Size, na.rm = TRUE) +
  2 * mad(neg_clean$Largest_Cluster_Size, constant = 1, na.rm = TRUE)

neg_mean <- mean(neg_clean$Largest_Cluster_Size, na.rm = TRUE)
neg_sd   <- sd(neg_clean$Largest_Cluster_Size, na.rm = TRUE)
cdf_neg  <- ecdf(neg_clean$Largest_Cluster_Size)

# Calculate energy-based thresholds from clean negative controls
p25_neg_E <- quantile(neg_clean$Best_Energy_Score, 0.25, na.rm = TRUE)
p05_neg_E <- quantile(neg_clean$Best_Energy_Score, 0.05, na.rm = TRUE)
mad_cut_E <- median(neg_clean$Best_Energy_Score, na.rm = TRUE) -
  2 * mad(neg_clean$Best_Energy_Score, constant = 1, na.rm = TRUE)


docking_annot <- docking_data %>%
  mutate(
    z_vs_neg   = (Largest_Cluster_Size - neg_mean)/neg_sd,
    pct_vs_neg = cdf_neg(Largest_Cluster_Size)
  )

# Plot with vertical lines at clean thresholds
ggplot(docking_annot,
       aes(x = Largest_Cluster_Size, y = Best_Energy_Score,
           color = Group, shape = Group, label = Interaction)) +
  geom_vline(xintercept = as.numeric(p75_neg), linetype = 2, linewidth = 0.5) +
  geom_vline(xintercept = as.numeric(p95_neg), linetype = 2, linewidth = 0.5) +
  geom_vline(xintercept = as.numeric(mad_cut),  linetype = 3, linewidth = 0.5) +
  annotate("text", x = as.numeric(p75_neg), y = min(docking_annot$Best_Energy_Score)+40,
           label = "       75th %ile (clean negatives)", angle = 90, vjust = -0.4, size = 3) +
  annotate("text", x = as.numeric(p95_neg), y = min(docking_annot$Best_Energy_Score)+40,
           label = "       95th %ile (clean negatives)", angle = 90, vjust = -0.4, size = 3) +
  annotate("text", x = as.numeric(mad_cut), y = min(docking_annot$Best_Energy_Score)+40,
           label = "               Robust cutoff (median+2×MAD)", angle = 90, vjust = -0.4, size = 3) +
  # Horizontal threshold lines (energy)
  geom_hline(yintercept = as.numeric(p25_neg_E), linetype = 2, linewidth = 0.7) +
  geom_hline(yintercept = as.numeric(p05_neg_E), linetype = 2, linewidth = 0.7) +
  geom_hline(yintercept = as.numeric(mad_cut_E),  linetype = 3, linewidth = 0.7) +
  
  # Horizontal threshold labels
  annotate("text",
           x = max(docking_annot$Largest_Cluster_Size) + 2,
           y = as.numeric(p25_neg_E),
           label = "25th %ile (negatives, energy)",
           hjust = 1, vjust = -0.4, size = 5) +
  annotate("text",
           x = max(docking_annot$Largest_Cluster_Size) + 2,
           y = as.numeric(p05_neg_E),
           label = "5th %ile (negatives, energy)",
           hjust = 1, vjust = -0.4, size = 5) +
  annotate("text",
           x = max(docking_annot$Largest_Cluster_Size) + 2,
           y = as.numeric(mad_cut_E),
           label = "Robust cutoff (median − 2×MAD)",
           hjust = 1, vjust = -0.4, size = 5) +
  
  geom_point(size = ifelse(docking_annot$Group=="PA187X-PA1874", 4.5, 3)) +
  geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.3,
                  min.segment.length = 0, seed = 1) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shape_map) +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(15, 60, 15, 15),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.97, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.box.background = element_blank()
  ) +
  labs(
    title = "Docking convergence: cluster size vs energy",
    subtitle = "Vertical lines computed from CLEAN negatives only; cluster size = primary reliability metric",
    x = "Largest cluster size (higher = more reliable)",
    y = "Docking score (more negative = lower energy)"
  )

ggplot(docking_annot,
       aes(x = Largest_Cluster_Size, y = Best_Energy_Score,
           color = Group, shape = Group, label = Interaction)) +
  # Vertical threshold lines
  geom_vline(xintercept = as.numeric(p75_neg), linetype = 2, linewidth = 0.7) +
  geom_vline(xintercept = as.numeric(p95_neg), linetype = 2, linewidth = 0.7) +
  geom_vline(xintercept = as.numeric(mad_cut),  linetype = 3, linewidth = 0.7) +
  
  # Threshold labels
  annotate("text", x = as.numeric(p75_neg), y = min(docking_annot$Best_Energy_Score) + 40,
           label = "                  75th %ile (negatives)", angle = 90, vjust = -0.4, size = 5) +
  annotate("text", x = as.numeric(p95_neg), y = min(docking_annot$Best_Energy_Score) + 40,
           label = "                  95th %ile (negatives)", angle = 90, vjust = -0.4, size = 5) +
  annotate("text", x = as.numeric(mad_cut), y = min(docking_annot$Best_Energy_Score) + 40,
           label = "                               Robust cutoff (median+2×MAD)", angle = 90, vjust = -0.4, size = 5) +
  
  # Horizontal threshold lines (energy)
  geom_hline(yintercept = as.numeric(p25_neg_E), linetype = 2, linewidth = 0.7) +
  geom_hline(yintercept = as.numeric(p05_neg_E), linetype = 2, linewidth = 0.7) +
  geom_hline(yintercept = as.numeric(mad_cut_E),  linetype = 3, linewidth = 0.7) +
  
  # Horizontal threshold labels
  annotate("text",
           x = max(docking_annot$Largest_Cluster_Size) + 2,
           y = as.numeric(p25_neg_E),
           label = "        25th %ile (negatives)",
           hjust = 1, vjust = -0.4, size = 5) +
  annotate("text",
           x = max(docking_annot$Largest_Cluster_Size) + 2,
           y = as.numeric(p05_neg_E),
           label = "        5th %ile (negatives)",
           hjust = 1, vjust = -0.4, size = 5) +
  annotate("text",
           x = max(docking_annot$Largest_Cluster_Size) + 2,
           y = as.numeric(mad_cut_E),
           label = "        Robust cutoff (median − 2×MAD)",
           hjust = 1, vjust = -0.4, size = 5) +
  
  
  # Points and point labels
  geom_point(size = ifelse(docking_annot$Group == "PA187X-PA1874", 6, 4.5)) +  # BIGGER POINTS
  geom_text_repel(size = 5,                          # BIGGER POINT LABELS
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  min.segment.length = 0,
                  seed = 1) +
  
  # Colors and shapes
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shape_map) +
  
  # Theme adjustments
  theme_minimal(base_size = 18) +  # boost base font size
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(15, 60, 15, 15),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = c(0.97, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 16),
    plot.title   = element_text(size = 22, face = "bold"),
    plot.subtitle= element_text(size = 16)
  ) +
  
  labs(
    title = "Docking convergence: cluster size vs energy",
    subtitle = "Vertical lines computed from CLEAN negatives only; cluster size = primary reliability metric",
    x = "Largest cluster size (higher = more reliable)",
    y = "Docking score\n(more negative = lower energy = stronger predicted binding)"
  )

ggsave(
  filename = "docking_convergence.png",
  plot = last_plot(),
  width = 16,       # inches
  height = 8,      # inches
  dpi = 600        # high resolution
)

# ===========================================
# 📂 15. Plot PA1874 SNP Appearances Across Patients and Years
# ===========================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(ggplot2); library(purrr)
})

patient_map <- c(
  "A005"="A005", "A011"="A011", "A012"="A012", "A025"="A025","A026"="A026","A027"="A027","A028A067"="A028A067","A032"="A032","A033"="A033","A037"="A037","A041"="A041","A042"="A042","A043"="A043","A045"="A045","A047"="A047","A066"="A066","A078"="A078","A086"="A086","A091"="A091","A097"="A097","A110"="A110","A120"="A120","A122"="A122","A128"="A128","A130"="A130","A131"="A131","A132"="A132","A133"="A133","A138"="A138","A147"="A147","A156"="A156","A158"="A158","A169"="A169","A185"="A185","A214"="A214","A281"="A281","A304"="A304","A510"="A510"
)

# PA1874 coordinates 
pa1874_start <- 3566483
pa1874_end   <- 3572893

# guide lines 
hotspots <- c(3567195, 3571943)

# 38 fixed colours 
patient_cols <- c(
  "#b8604c", "#4cc75a", "#b45fd2", "#87bb37", "#6f5fd2",
  "#c4b93d", "#497adb", "#e09332", "#53a2d7", "#df6834",
  "#47caca", "#d8423d", "#5cb764", "#d14fb1", "#4c902b",
  "#e7488b", "#5dc391", "#b13b78", "#448b50", "#cc405d",
  "#369780", "#a3421d", "#a595e0", "#878f2a", "#8d539d",
  "#af892a", "#5667a5", "#9fb86f", "#de88be", "#52712a",
  "#954c6f", "#2a6a45", "#e88883", "#666020", "#a7525f",
  "#938a49", "#d9a26a", "#95632e"
)

# ---------- HELPERS ----------
# 1) Find all years present in a genome folder, robust to different filename styles
discover_years <- function(genome_id) {
  files <- list.files(genome_id, full.names = FALSE)
  files <- files[grepl("^(breseq|marginal|snps)_", files, ignore.case = TRUE)]
  # Grab a 4-digit year anywhere after an underscore, before extension
  yrs <- str_match(files, "_(\\d{4})\\.(?:gd|html)$")[,2]
  yrs <- sort(unique(na.omit(as.integer(yrs))))
  as.character(yrs)
}

# 2) Breseq "real" SNPs from .gd
extract_breseq_real <- function(gd_path) {
  if (!file.exists(gd_path)) return(integer(0))
  gd <- read.table(gd_path, header = FALSE, sep = "\t",
                   comment.char = "#", fill = TRUE, quote = "")
  pos <- suppressWarnings(as.numeric(gd[gd$V1 %in% c("SNP"), 5]))
  pos[!is.na(pos)]
}

# 3) Breseq "marginal" SNPs from marginal_*.html
extract_breseq_marginal <- function(html_path) {
  if (!file.exists(html_path)) return(integer(0))
  html_lines <- readLines(html_path, warn = FALSE)
  tds <- regmatches(html_lines, gregexpr("<td[^>]*>[0-9,]+</td>", html_lines))
  tds <- unlist(tds)
  nums <- as.numeric(gsub("[^0-9]", "", tds))
  pos <- unique(nums[seq(1, length(nums), 2)])
  pos[!is.na(pos)]
}

# 4) Snippy SNPs from snps_*.html (POS is 2nd <TD> in 6-col table)
extract_snippy <- function(html_path) {
  if (!file.exists(html_path)) return(integer(0))
  html_lines <- readLines(html_path, warn = FALSE)
  td_lines <- grep("<TD>", html_lines, value = TRUE)
  if (length(td_lines) == 0) return(integer(0))
  ncol <- 6
  nrow <- length(td_lines) %/% ncol
  if (nrow == 0) return(integer(0))
  td_matrix <- matrix(td_lines[seq_len(nrow * ncol)], ncol = ncol, byrow = TRUE)
  pos <- as.numeric(gsub("<.*?>", "", td_matrix[, 2]))
  pos[!is.na(pos)]
}

# 5) Process one genome folder across all its years
process_one_genome <- function(genome_id, patient_id) {
  years <- discover_years(genome_id)
  if (length(years) == 0) {
    message("No years found for ", genome_id, " — skipping.")
    return(tibble())
  }
  bind_rows(lapply(years, function(year) {
    gd_path        <- file.path(genome_id, paste0("breseq_",  year, ".gd"))
    if (!file.exists(gd_path)) {
      # try alternative naming with genome or sample in the middle
      cand <- list.files(genome_id, pattern = paste0("^breseq_.*_", year, "\\.gd$"), full.names = TRUE)
      if (length(cand)) gd_path <- cand[1]
    }
    marg_html_path <- file.path(genome_id, paste0("marginal_", year, ".html"))
    if (!file.exists(marg_html_path)) {
      cand <- list.files(genome_id, pattern = paste0("^marginal_.*_", year, "\\.html$"), full.names = TRUE)
      if (length(cand)) marg_html_path <- cand[1]
    }
    snippy_html    <- list.files(genome_id, pattern = paste0("^snps_.*_", year, "\\.html$"),
                                 full.names = TRUE)
    snippy_html    <- if (length(snippy_html)) snippy_html[1] else file.path(genome_id, paste0("snps_", year, ".html"))
    
    b_real <- extract_breseq_real(gd_path)
    b_marg <- extract_breseq_marginal(marg_html_path)
    s_pos  <- extract_snippy(snippy_html)
    
    # restrict to PA1874 region
    b_real <- b_real[b_real >= pa1874_start & b_real <= pa1874_end]
    b_marg <- b_marg[b_marg >= pa1874_start & b_marg <= pa1874_end]
    s_pos  <- s_pos [s_pos  >= pa1874_start & s_pos  <= pa1874_end]
    
    run_id <- paste0(genome_id, "_", year)
    
    bind_rows(
      tibble(Run = run_id, Patient = patient_id, Year = as.integer(year),
             Position = b_real, CallType = "Breseq_real"),
      tibble(Run = run_id, Patient = patient_id, Year = as.integer(year),
             Position = b_marg, CallType = "Breseq_marginal"),
      tibble(Run = run_id, Patient = patient_id, Year = as.integer(year),
             Position = s_pos,  CallType = "Snippy")
    )
  }))
}

# Build the full table
all_df <- imap_dfr(patient_map, ~ process_one_genome(genome_id = .y, patient_id = .x))

# Order Patient factor to match colour vector order
# (Assumes names(patient_map) length == length(patient_cols))
patients_in_order <- unname(patient_map)
all_df <- all_df %>%
  filter(!is.na(Position)) %>%
  mutate(
    Patient = factor(Patient, levels = patients_in_order),
    CallType = factor(CallType, levels = c("Breseq_real","Breseq_marginal","Snippy"))
  )

write_csv(all_df, "pa1874_snps_all_patients.csv")
message("Wrote: pa1874_snps_all_patients.csv with ", nrow(all_df), " rows")

# Plot
p <- ggplot(all_df, aes(x = Position, y = Year, colour = Patient, shape = CallType)) +
  geom_point(alpha = 0.9, size = 2.5) +
  scale_shape_manual(values = c(Breseq_real = 16, Breseq_marginal = 2, Snippy = 15)) +
  scale_colour_manual(values = patient_cols, drop = FALSE) +
  labs(title = "PA1874 SNP appearances — patients coloured, callers shaped",
       x = "PA1874 genomic position", y = "Collection year") +
  theme_minimal(base_size = 12) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 3),
      nrow = 2  # Patient legend in 2 rows
    ),
    shape = guide_legend(
      nrow = 1  # CallType legend in 1 row
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p

p <- p +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(colour = "grey95", linewidth = 0.2)
  )

p

p <- p +
  scale_y_continuous(
    breaks = c(1990, 1995, 2000, 2005, 2010, 2015, 2020),  
    limits = NULL)+                         # show full range incl. 1995
      scale_x_continuous(
        breaks = c(3567000, 3568000, 3569000, 3570000, 3571000, 3572000)
      )+
  theme(
    panel.grid.minor = element_blank(),  # remove minor grid lines
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5) # keep major grid lines
  )
    
p

p <- p +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 20, b = 40, l = 10))

p

ggsave("pa1874_snps_scatter_all_patients_2.png", p, width = 20, height = 6, dpi = 300)


# ===========================================
# 📂 16. Explore Marginal Density of PA1874 SNP Positions and Collection Years
# ===========================================
# Density plot of SNP collection years
p2 <- ggplot(all_df, aes(x = Year)) +
  geom_density(fill = "#C6DBEF", alpha = 0.3) +
  scale_x_continuous(
    breaks = c(1990, 1995, 2000, 2005, 2010, 2015, 2020)
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5) ,
    axis.title.y = element_blank(),  # remove redundant axis label
    axis.text.y = element_blank(),   # remove ticks so it’s clean beside plot
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),   # remove ticks so it’s clean beside plot
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
  )

p2

# Alternative density visualisation for year distribution
p2 <- ggplot(all_df, aes(x = Year)) +
  coord_flip() +
  geom_density(aes(y = ..density.. / max(..density..)), 
               fill = "#C6DBEF", alpha = 0.3) +
  scale_y_continuous(name = "Normalised Density (max = 1)") +
  scale_x_continuous(breaks = seq(1990, 2020, 5)) +
  theme_minimal()

p2

p2 <- ggplot(all_df, aes(x = Year)) +
  coord_flip() +
  geom_density(fill = "#C6DBEF", alpha = 0.3) +
  scale_x_continuous(breaks = seq(1990, 2020, 5)) +
  theme_minimal()
p2

# ===========================================
# 📂 17. Plot RTX Motif Scan Results Across Adhesins and PA1874
# ===========================================
library(ggplot2)

# Define RTX motif counts from ScanProsite results
rtx <- data.frame(
  protein = c("LapA (P. putida)", "FrhA (V. cholerae)", "RtxA (L. pneumophila)",
              "BpfA (S. oneidensis)", "LapF (P. putida)", "PA1874 (P. aeruginosa)"),
  hits    = c(4, 3, 3, 3, 1, 0)
)

# Order proteins by motif count and highlight PA1874 separately
rtx$protein <- factor(rtx$protein,
                      levels = c("LapA (P. putida)", "FrhA (V. cholerae)", "RtxA (L. pneumophila)",
                                 "BpfA (S. oneidensis)", "LapF (P. putida)", "PA1874 (P. aeruginosa)"))
rtx$group <- ifelse(grepl("^PA1874", rtx$protein), "PA1874", "Other")

p <- ggplot(rtx, aes(x = protein, y = hits, fill = group)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.3) +
  geom_text(aes(label = hits), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = c(Other = "#A7C7E7", PA1874 = "#2E6FBE"), guide = "none") +
  labs(x = NULL, y = "Number of RTX motifs detected",
       title = "RTX motif scan across adhesins vs PA1874") +
  coord_cartesian(ylim = c(0, max(rtx$hits) + 0.8)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
R.version

p