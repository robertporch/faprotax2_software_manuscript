#### /usr/bin/Rscript
#
#
#
# This is an R script for performing cross-validation between FAPROTAX v1 and v2.
# There is a command-line argument for specifying whether to perform taxonomy- or sequence-based placement
# in FAPROTAX v2.
#
#
# Rob Porch
# last updated 2026 March 17


###################
# Load packages
library(data.table)
library(seqinr)
library(stringr)
library(tidyverse)
library(seqinr)
library(optparse)
###################

####################
# Set up command-line arguments

option_list <- list(
	make_option(c("-i", "--in_faprotax2_table"), type="character", default="", help="Path to FAPROTAX2 phenotypic database."),
	make_option("--in_faprotax1_table", type="character", default="", help="Path to FAPROTAX v1 phenotypic database, at species-level."),
	make_option("--in_shared_functions", type="character", default="", help="Path to .txt file listing shared functions between FAPROTAX v1 and v2, in order to limit functions for cross-validation."),
	make_option("--minimum_absent", type="integer", default=100, help="Minimum number of verified absence records for which phenotypic data exist in order to consider the function for cross-validation. (default %default)"),
	make_option("--minimum_present", type="integer", default=100, help="Minimum number of verified presence records for which phenotypic data exist in order to consider the function for cross-validation. (default %default)"),
	make_option("--test_pool_fraction", type="double", default=0.05, help="Fraction of verified presence/absence records to sample as a test pool for cross-validation. (default %default)"),
	make_option(c("-m", "--placement_method"), type="character", default="taxonomy", help="Method of query-to-reference placement for FAPROTAX2. Options are either 'taxonomy' or 'sequence'. (default %taxonomy)"),
	make_option(c("-s", "--in_sequences"), type="character", default="", help="Path to file containing representative SILVA sequences, for the purposes of extracting sequences from test-set species."),
	make_option(c("-n", "--n_rounds"), type="integer", default=10, help="Number of rounds of cross-validation to perform, per function. (default %default)"),
	make_option(c("-o", "--out_final_summary"), type="character", default="", help="Path to store output figure showing CV results.")
)

# Parse command-line arguments
opt <- parse_args(optparse::OptionParser(option_list = option_list))


# --- Setup Paths ---
set.seed(543)
if (opt$placement_method == "taxonomy") {
	project_dir <- "output/cross_validation/taxonomy"
} else if (opt$placement_method == "sequencee") {
	project_dir <- "output/cross_validation/sequence"
} else {
	stop("Error: must specify desired placement method ('taxonomy' or 'placement') for FAPROTAX2, using --placement_method")
}

results_dir <- file.path(project_dir, "results")
dir.create(results_dir, showWarnings = FALSE)

# --- Load Shared Functions ---
# Load the list of shared functions and clean the V2 database names
shared_traits_list <- readLines(opt$in_shared_functions)

# Load V2 Database
FAPROTAX2_data <- data.table::fread(opt$in_faprotax2_table) %>%
distinct(accession, .keep_all = TRUE)

FAPROTAX1_data <- data.table::fread(opt$in_faprotax1_table) %>%
distinct(accession, .keep_all = TRUE)

# Subset to only include novel species records from FAPROTAX2 (i.e. remove species records present in FAPROTAX1).
# This step ensures that the test set species are not present in the FAPROTAX v1 database and therefore will not be
# predicted by FAPROTAX v1 with 100% true positive rate (because we are not removing taxon records from FAPROTAX v1 in this analysis).

FAPROTAX2_data <- FAPROTAX2_data %>%
filter(!(species %in% FAPROTAX1_data$species))

# Map database columns to clean names (removing .value)
actual_trait_cols <- colnames(FAPROTAX2_data)[grepl("\\.value$", colnames(FAPROTAX2_data))]
trait_mapping <- actual_trait_cols
names(trait_mapping) <- sub("\\.value$", "", actual_trait_cols)


# Only keep traits that exist in the shared list
valid_traits <- intersect(names(trait_mapping), shared_traits_list)

# --- Data Pre-processing ---
# Convert P/A/U to 1/0/NA
FAPROTAX2_data <- FAPROTAX2_data %>%
mutate(across(all_of(actual_trait_cols), ~ case_when(
  . == "P" ~ 1,
  . == "A" ~ 0,
  TRUE ~ NA_real_
)))


#### Load in SILVA sequences from FAPROTAX2-db
if (opt$placement_method == "sequence") {
	message("Loading master SILVA sequences..")
	all_sequences <- seqinr::read.fasta(file=opt$in_sequences,
                                    	forceDNAtolower = FALSE)
} else {}

# --- Helper Functions ---
evaluate_metrics <- function(truth, pred, n_present_total, n_absent_total) {
truth <- as.numeric(truth)
pred <- as.numeric(pred)
min_present <- opt$minimum_present
min_absent <- opt$minimum_absent

tp <- sum(pred == 1 & truth == 1, na.rm = TRUE)
tn <- sum(pred == 0 & truth == 0, na.rm = TRUE)
fp <- sum(pred == 1 & truth == 0, na.rm = TRUE)
fn <- sum(pred == 0 & truth == 1, na.rm = TRUE)

# STRICT FILTERING: 
# Only calculate TPR if the total pool of known presences is >= 50
tpr <- if(n_present_total >= min_present) tp / sum(truth == 1, na.rm = TRUE) else NA

# Only calculate TNR if the total pool of known absences is >= 50
tnr <- if(n_absent_total >= min_absent) tn / sum(truth == 0, na.rm = TRUE) else NA

# Balanced Accuracy is ONLY valid if both metrics met the threshold
ba <- if(!is.na(tpr) && !is.na(tnr)) (tpr + tnr) / 2 else NA

return(list(TPR = tpr, TNR = tnr, BA = ba))
}

# --- Initialize Result Containers ---
round_result <- data.table()
true_vs_pred <- data.table()
n_rounds <- opt$n_rounds
min_present <- opt$minimum_present
min_absent <- opt$minimum_absent

# --- Cross-Validation Loop ---
for (round in 1:n_rounds) {
set.seed(54 + round)
message("Starting round ", round, "...")

for (trait in valid_traits) {
  real_col_name <- trait_mapping[trait]
  
  # 1. Filter Logic: Use real_col_name for subsetting
  subdb <- FAPROTAX2_data[, .SD, .SDcols = c("accession", "taxonomy", real_col_name)]
  subdb <- subdb[!is.na(get(real_col_name))]
  
  n_present_total <- sum(subdb[[real_col_name]] == 1, na.rm = TRUE)
  n_absent_total <- sum(subdb[[real_col_name]] == 0, na.rm = TRUE)
  
  if (n_present_total < min_present && n_absent_total < min_absent) {
    message("Trait: ", trait, " has insufficient records (P:", n_present_total, " A:", n_absent_total, "). Skipping..")
    next
  }
  
  # 2. Sampling Logic (10% or up to 500)
  present_pool <- subdb[subdb[[real_col_name]] == 1, ]
  absent_pool <- subdb[subdb[[real_col_name]] == 0, ]
  
  n_sample_pres <- if (nrow(present_pool) > 0) {
    floor(opt$test_pool_fraction * nrow(present_pool))
  } else 0
  
  n_sample_abs <- if (nrow(absent_pool) > 0) {
    floor(opt$test_pool_fraction * nrow(absent_pool))
  } else 0
  
  sampled_pres <- if (n_sample_pres > 0) {
    present_pool[sample(.N, n_sample_pres)]
  } else data.table()
  
  sampled_abs <- if (n_sample_abs > 0) {
    absent_pool[sample(.N, n_sample_abs)]
  } else data.table()
  
  sampled <- rbind(sampled_pres, sampled_abs)
  sampled[, OTU := as.character(seq_len(.N))]
  setnames(sampled, real_col_name, trait)
  
  test_accessions <- as.character(sampled$accession)
  
  if (opt$placement_method == "sequence") {
  	test_sequences <- all_sequences[test_accessions]
  	names(test_sequences) <- sampled$OTU
  } else {}
  
  # --- Create Folders ---
  trait_dir <- file.path(results_dir, trait)
  input_dir <- file.path(trait_dir, "input")
  output_dir_v1 <- file.path(trait_dir, "output/v1")
  output_dir_v2 <- file.path(trait_dir, "output/v2")
  
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_v1, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_v2, recursive = TRUE, showWarnings = FALSE)
  
  # --- Prepare Tool-Specific Data Tables ---
  
  # 1. FAPROTAX v1 full table (Taxonomy + Sample dummy column)
  # .SDcols excludes the 'accession' and 'OTU' columns to match your original flow
  sampled_faprotax_v1 <- sampled[, .(OTU, taxonomy, Sample_01 = 1)]
  # Note: If your v1 script needs the trait column in the input, use:
  # sampled_faprotax_v1 <- sampled[, .(taxonomy, get(trait), Sample_01 = 1)]
  
  # 2. FAPROTAX v2 table (OTU ID + Sample dummy column)
  sampled_faprotax_v2 <- data.table(OTU = sampled$OTU, Sample_01 = 1)
  
  # --- Save Input Files (The "Thorough" part) ---
  
  # Save the simple OTU table (OTU ID and dummy abundance)
  data.table::fwrite(sampled_faprotax_v2, file.path(input_dir, "otu_table.tsv"), sep = "\t")
  
  # Save the full table (usually used for reference or specific v1 parsing)
  data.table::fwrite(sampled_faprotax_v1, file.path(input_dir, "otu_table_faprotaxv1.tsv"), sep = "\t")
  
  # Save the taxonomy lookup
  data.table::fwrite(sampled[, .(OTU, taxonomy)], file.path(input_dir, "otu_taxonomy_table.tsv"), sep = "\t")
  
  # Save the accession lookup
  data.table::fwrite(sampled[, .(accession)], file.path(input_dir, "otu_accession_table.tsv"), sep = "\t")
  
  # Save the SILVA sequences corresponding to each accession
  if (opt$placement_method == "sequence") {
  	seqinr::write.fasta(sequences = test_sequences,
                      	names = names(test_sequences),
                      	file.out = file.path(input_dir, "otu_alignments.fasta"))
  }
  
  # Save Ground Truth (OTU, accession, taxonomy, and the focal trait)
  # Using 'trait' here because we renamed the column above
  truth_to_save <- sampled[, .(OTU, accession, taxonomy, get(trait))]
  setnames(truth_to_save, "V4", trait)
  data.table::fwrite(truth_to_save, file.path(input_dir, "true_traits.tsv"), sep = "\t")
  
  # Save accessions for removal (the "Leave-one-out" equivalent)
  writeLines(as.character(sampled$accession), file.path(input_dir, "accessions_to_remove.txt"))
  
  # --- Command Line Execution ---
  
  # V1 COMMAND
  # Note: Fill in your specific flags here as per your Version 1 requirements
  message("Running FAPROTAX v1 on trait: ", trait, " in round ", round, " of ", n_rounds, "..")
  v1_cmd <- paste(
    c("dependencies/collapse_table.py",
      "-f",
      "-v",
 	    "-i", file.path(input_dir, "otu_table_faprotaxv1.tsv"),
 	    "--row_names_are_in_column", "taxonomy",
 	    "-g", "dependencies/FAPROTAX.txt",
 	    "--out_groups2records_table", file.path(output_dir_v1, "ASV2function_mapping.tsv"),
      "-n", "columns_before_collapsing"),
    collapse = " "
  )
  
  # Dynamic logic for FAPROTAX v2, depending on placement method
  v2_input_flag <- if(opt$placement_method == "sequence") "-a" else "-t"
  v2_input_file <- if(opt$placement_method == "sequence") {
  	file.path(input_dir, "otu_alignments.fasta")
  } else if (opt$placement_method == "taxonomy") {
  	paste0("file:", file.path(input_dir, "otu_taxonomy_table.tsv"))
  } else {
  	stop("Error: must specify desired placement method ('taxonomy' or 'placement') for FAPROTAX2, using --placement_method")
  }
  
  # V2 COMMAND
  message("Running FAPROTAX v2 on trait: ", trait, " in round ", round, " of ", n_rounds, "..")
  v2_cmd <- paste(
    c("Rscript", "dependencies/faprotax.R",
      "-d", "dependencies/faprotax",
      "-c", "dependencies/faprotax/default_config.tsv",
      "-i", file.path(input_dir, "otu_table.tsv"),
      v2_input_flag, v2_input_file,
      "--omit_references_in_file", file.path(input_dir, "otu_accession_table.tsv"),
      "--only_functions", trait,
      "-r", file.path(output_dir_v2, "report.txt"),
      "-l", file.path(output_dir_v2, "log.txt"),
      "-o", file.path(output_dir_v2, "function_table.tsv"),
      "--out_otu_function_table", file.path(output_dir_v2, "ASV2function_mapping.tsv"),
      "--otu_names_are_in_colum", "OTU",
      "-f",
      "-v",
      "--out_intermediates", file.path(output_dir_v2, "intermediates"),
      "--include_summary_comments",
      "--out_nearest_refs", file.path(output_dir_v2, "nearest_references.tsv"),
      "--out_annotation_probabilities", file.path(output_dir_v2, "annotation_probabilities.tsv"),
      "--weigh_by_probability",
      "--normalize", "binarize_before_collapsing",
      "--Nthreads", "1"),
    collapse = " "
  )
  
  system(v1_cmd)
  system(v2_cmd)
  
  # --- Parsing Results ---
  v1_res <- fread(file.path(output_dir_v1, "ASV2function_mapping.tsv"))
  v2_res <- fread(file.path(output_dir_v2, "ASV2function_mapping.tsv"))
  # Check if the trait column exists in both (safety check)
  if (!(trait %in% colnames(v1_res)) | !(trait %in% colnames(v2_res))) {
    message("Trait ", trait, " missing from output columns. Skipping..")
    next
  }
  
  # V1: Ignore 'record' column, use row index to create "1", "2", "3"...
  # We convert to character to match the 'OTU' column in our 'sampled' data.table
  v1_preds <- v1_res[, .(OTU = as.character(1:.N), pred_v1 = get(trait))]
  
  # V2: Use same row-index logic to ensure perfect alignment 
  v2_preds <- v2_res[, .(OTU = as.character(1:.N), pred_v2 = get(trait))]
  
  # Join predictions with the ground truth
  # We use 'sampled' which we already prepared with OTU = seq_len(.N)
  compare_df <- sampled %>%
    mutate(OTU = as.character(OTU)) %>%
    left_join(v1_preds, by = "OTU") %>%
    left_join(v2_preds, by = "OTU")
  
  # Fill NAs in predictions as 0
  compare_df$pred_v1[is.na(compare_df$pred_v1)] <- 0
  compare_df$pred_v2[is.na(compare_df$pred_v2)] <- 0
  
  # --- Store Raw Predictions ---
  true_vs_pred <- rbind(true_vs_pred, compare_df[, .(
    OTU, accession, taxonomy, trait = trait, round = round, 
    true_state = get(trait), pred_v1, pred_v2
  )])
  
  # Calculate round metrics
  m1 <- evaluate_metrics(compare_df[[trait]], compare_df$pred_v1, n_present_total, n_absent_total)
  m2 <- evaluate_metrics(compare_df[[trait]], compare_df$pred_v2, n_present_total, n_absent_total)
  
  message("Storing cross-validation results for round ", round, " of ", n_rounds, "..")
  round_result <- rbind(round_result, data.table(
    trait = trait, round = round,
    v1_TPR = m1$TPR, v1_TNR = m1$TNR, v1_BA = m1$BA,
    v2_TPR = m2$TPR, v2_TNR = m2$TNR, v2_BA = m2$BA
  ))
}
}

# --- Final Aggregation ---
final_summary <- round_result %>%
group_by(trait) %>%
summarise(
  mean_v1_TPR = mean(v1_TPR, na.rm = TRUE),
  mean_v1_TNR = mean(v1_TNR, na.rm = TRUE),
  mean_v1_BA  = mean(v1_BA, na.rm = TRUE),
  mean_v2_TPR = mean(v2_TPR, na.rm = TRUE),
  mean_v2_TNR = mean(v2_TNR, na.rm = TRUE),
  mean_v2_BA  = mean(v2_BA, na.rm = TRUE)
)

fwrite(final_summary, opt$out_final_summary, sep = "\t")