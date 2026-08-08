#!/bin/bash

# ==========================================
# MitoCarta RPA Pipeline v1.0.0
# ==========================================

# ==========================================
# Help Display Function
# ==========================================
show_help() {
    cat << EOF
[Usage]
  bash $0 <input_file> <human|mouse>

[Description]
  Performs mitochondrial pathway enrichment analysis using the Broad Institute's MitoCarta3.0 database.
  
  Workflow & Features:
  - Automatically downloads the specified species dataset (Human or Mouse) directly from the Broad Institute.
  - Maps input genes (requires a 'Symbol' or 'GeneSymbol' column) to MitoCarta pathways.
  - Calculates statistical significance (P-value via Hypergeometric test) and False Discovery Rate (FDR via Benjamini-Hochberg method).
  - If a 'Dominant_Direction' column (e.g., UP/DOWN) is detected in the input, the analysis is automatically split by expression direction.
  - Runs entirely inside a Docker container, ensuring high reproducibility and eliminating local R dependencies.
  
  [Output Files (2)]
  1. Mito_hit_pathway.txt : Pathway aggregation and statistical testing results (Universe = ~20,000 total genes).
  2. Mito_hit_list.txt    : Detailed information of hit genes (Maintains the exact row order of your input data).

[Arguments]
  <input_file> File to be analyzed (CSV, TSV, TXT, etc.)
  <human|mouse> Target species ('human' or 'mouse')
EOF
}

# Display help if requested or if arguments are missing
if [[ "$1" == "-h" || "$1" == "--help" || $# -lt 2 ]]; then
    show_help
    exit 0
fi

INPUT_FILE="$1"
SPECIES=$(echo "$2" | tr '[:upper:]' '[:lower:]')

# Validate species argument
if [[ "$SPECIES" != "human" && "$SPECIES" != "mouse" ]]; then
    echo "Error: Please specify 'human' or 'mouse' for the second argument."
    exit 1
fi

INPUT_ABS=$(realpath "$INPUT_FILE")
WORK_DIR=$(pwd)

echo "========================================"
echo " MitoCarta RPA Pipeline (2-Files Output)"
echo " Version     : v1.0.0"
echo " Target File : $INPUT_FILE"
echo " Species     : $SPECIES"
echo "========================================"

docker run -i --rm \
    -v "$INPUT_ABS":/input_file:ro \
    -v "$WORK_DIR":/work \
    ezojika7713/rnaseq-meta:latest \
    Rscript - "$SPECIES" "/input_file" "/work" << 'EOF'

args <- commandArgs(trailingOnly = TRUE)
species    <- args[1]
input_path <- args[2]
work_dir   <- args[3]

options(timeout = 600)

library(readxl)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(tools)

# --- Define species-specific parameters ---
if (species == "human") {
  url <- "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"
  xls_filename <- "Human.MitoCarta3.0.xls"
  sheet_name <- "A Human MitoCarta3.0"
} else {
  url <- "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls"
  xls_filename <- "Mouse.MitoCarta3.0.xls"
  sheet_name <- "A Mouse MitoCarta3.0"
}

xls_path <- file.path(work_dir, xls_filename)

# Download MitoCarta database if it does not exist locally
if (!file.exists(xls_path)) {
  cat(sprintf("Downloading %s ...\n", xls_filename))
  download.file(url, destfile = xls_path, mode = "wb", quiet = FALSE, headers = c("User-Agent" = "Mozilla/5.0"))
} else {
  cat(sprintf("Using local file: %s\n", xls_filename))
}

# --- Load MitoCarta Database ---
mito_db <- read_excel(xls_path, sheet = sheet_name)

# Set the total number of genes (Universe) according to the species
if (species == "human") {
  total_universe_genes <- 19248
} else {
  total_universe_genes <- 20632
}

pathway_sizes <- mito_db %>%
  separate_rows(MitoCarta3.0_MitoPathways, sep = "\\s*\\|\\s*") %>%
  filter(!is.na(MitoCarta3.0_MitoPathways) & MitoCarta3.0_MitoPathways != "0") %>%
  group_by(MitoCarta3.0_MitoPathways) %>%
  summarise(Pathway_Total_Genes = n_distinct(Symbol), .groups = "drop")

# --- Load Input File and Record Original Row Order ---
ext <- tolower(file_ext(input_path))
if (ext == "csv") {
  deg_df <- read.csv(input_path, comment.char="#", stringsAsFactors=FALSE, check.names=FALSE)
} else {
  deg_df <- read.delim(input_path, sep="\t", comment.char="#", stringsAsFactors=FALSE, check.names=FALSE)
  if (ncol(deg_df) == 1 && any(grepl(",", deg_df[[1]]))) {
    deg_df <- read.csv(input_path, comment.char="#", stringsAsFactors=FALSE, check.names=FALSE)
  }
}

col_names <- colnames(deg_df)
if ("GeneSymbol" %in% col_names) {
  deg_df <- deg_df %>% rename(Symbol = GeneSymbol)
} else if (!"Symbol" %in% col_names) {
  if (ncol(deg_df) == 1) {
    first_gene <- colnames(deg_df)[1]
    deg_df <- data.frame(Symbol = c(first_gene, deg_df[[1]]), stringsAsFactors=FALSE)
  } else {
    stop("Error: 'Symbol' or 'GeneSymbol' column not found.")
  }
}

# Save the original row order of the input file
deg_df <- deg_df %>% mutate(Input_Order = row_number())

# --- Absorb Notation Variations in Dominant_Direction ---
if ("Dominant_Direction" %in% colnames(deg_df)) {
  deg_df$Dominant_Direction <- toupper(deg_df$Dominant_Direction)
  cat("Detected 'Dominant_Direction' column. Splitting analysis by UP/DOWN...\n")
} else {
  cat("No 'Dominant_Direction' column detected. Treating as a single gene list...\n")
}

total_input_genes <- n_distinct(deg_df$Symbol)

# ==========================================
# 1. Detailed Gene List (Sorted by Input)
# ==========================================
cat("Mapping genes and generating sorted list...\n")

# Join both datasets and sort by the original input file order
out_list_input <- deg_df %>%
  inner_join(mito_db, by = "Symbol") %>%
  arrange(Input_Order) %>%
  select(-Input_Order) # Remove temporary tracking column

# ==========================================
# 2. Pathway Enrichment Analysis
# ==========================================
cat("Performing Pathway Enrichment Analysis...\n")

pathway_mapped <- out_list_input %>%
  separate_rows(MitoCarta3.0_MitoPathways, sep = "\\s*\\|\\s*") %>%
  filter(!is.na(MitoCarta3.0_MitoPathways) & MitoCarta3.0_MitoPathways != "0")

if ("Dominant_Direction" %in% colnames(pathway_mapped)) {
  
  direction_totals_input <- deg_df %>%
    group_by(Dominant_Direction) %>%
    summarise(Sample_Size = n_distinct(Symbol), .groups = "drop")

  pathway_summary <- pathway_mapped %>%
    group_by(MitoCarta3.0_MitoPathways, Dominant_Direction) %>%
    summarise(
      Hit_Count = n_distinct(Symbol),
      Genes = paste(unique(Symbol), collapse = ", "),
      .groups = "drop"
    ) %>%
    left_join(pathway_sizes, by = "MitoCarta3.0_MitoPathways") %>%
    left_join(direction_totals_input, by = "Dominant_Direction") %>%
    mutate(
      P_value = phyper(Hit_Count - 1, Pathway_Total_Genes, total_universe_genes - Pathway_Total_Genes, Sample_Size, lower.tail = FALSE),
      FDR = p.adjust(P_value, method = "BH")
    ) %>%
    select(MitoCarta3.0_MitoPathways, Dominant_Direction, Hit_Count, Pathway_Total_Genes, P_value, FDR, Genes) %>%
    arrange(P_value)

} else {
  pathway_summary <- pathway_mapped %>%
    group_by(MitoCarta3.0_MitoPathways) %>%
    summarise(
      Hit_Count = n_distinct(Symbol),
      Genes = paste(unique(Symbol), collapse = ", "),
      .groups = "drop"
    ) %>%
    left_join(pathway_sizes, by = "MitoCarta3.0_MitoPathways") %>%
    mutate(
      P_value = phyper(Hit_Count - 1, Pathway_Total_Genes, total_universe_genes - Pathway_Total_Genes, total_input_genes, lower.tail = FALSE),
      FDR = p.adjust(P_value, method = "BH")
    ) %>%
    select(MitoCarta3.0_MitoPathways, Hit_Count, Pathway_Total_Genes, P_value, FDR, Genes) %>%
    arrange(P_value)
}

# --- Save Results ---
file_pathway    <- file.path(work_dir, "Mito_hit_pathway.txt")
file_out_input  <- file.path(work_dir, "Mito_hit_list.txt")

write.table(pathway_summary, file_pathway, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(out_list_input, file_out_input, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Success! Generated the following 2 files:\n")
cat("  1. Mito_hit_pathway.txt\n")
cat("  2. Mito_hit_list.txt\n")

EOF

echo "========================================"
echo "Analysis completed."