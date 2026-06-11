library(GeomxTools)
library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BiocParallel)
library(NanoStringNCTools)
library(lme4)
library(pbapply)
library(lmerTest)
library(readxl) 

# --- 1. Load Data ---
raw_expression <- fread("/Users/nikamikhailava/Desktop/Lab_stuff/batch_merged_new/fibroblasts_gene_expression_data_new_batch_three_merged_transposed.csv")
metadata <- read_excel("/Volumes/h30492/farkkilab2/9_EyeMT/Data/geomx/batch123/metadata/dcc_metadata_batch123_cleaned.xlsx")

# --- 2. Format Data ---
all_genes <- setdiff(colnames(raw_expression), "dcc_filename")

# Merge expression with metadata directly 
raw_merged <- raw_expression %>%
  left_join(metadata, by = "dcc_filename")

# --- 3. Filter for Stroma & TSI ---
raw_merged <- raw_merged %>%
  filter(Segment == "stroma", 
         Segment_geomx %in% c("stroma", "tsi"))

# Set the factor levels for region (stroma as baseline)
raw_merged$Segment_geomx <- factor(
  raw_merged$Segment_geomx,
  levels = c("stroma", "tsi")
)

# Set factor levels for NACT status
raw_merged$NACT_status <- factor(
  raw_merged$NACT_status,
  levels = c("pre", "post")
)

# --- 4. Define LMM Function ---
run_gene_test_satterthwaite <- function(gene, df) {
  tryCatch({
    y_values <- as.numeric(trimws(as.character(df[[gene]])))
    
    if (sum(!is.na(y_values)) < 3) {
      return(list(gene = gene, status = "error", message = "Not enough non-NA values", result = NA))
    }
    if (var(y_values, na.rm = TRUE) == 0) {
      return(list(gene = gene, status = "error", message = "Zero variance", result = NA))
    }
    
    # Fit model with NACT_status covariate
    mod <- lmer(y_values ~ Segment_geomx + NACT_status + (1 | Sample), data = df)
    summary_mod <- summary(mod)$coefficients
    
    # Extract contrast for 'tsi'
    target_row <- grep("Segment_geomxtsi", rownames(summary_mod), value = TRUE)
    
    if (length(target_row) == 0) {
      return(list(gene = gene, status = "error", message = "Contrast not found", result = NA))
    }
    
    stats <- summary_mod[target_row, ]
    return(list(
      gene = gene, status = "success", message = NA,
      result = c(Estimate = stats["Estimate"], Std_Error = stats["Std. Error"],
                 df = stats["df"], t_value = stats["t value"], P_Value = stats["Pr(>|t|)"])
    ))
  }, error = function(e) {
    return(list(gene = gene, status = "error", message = e$message, result = NA))
  })
}

# --- 5. Run Analysis ---
cat("Starting adjusted LMM analysis for", length(all_genes), "genes...\n")
results_list <- pblapply(all_genes, run_gene_test_satterthwaite, df = raw_merged)

# --- 6. Process Results ---
names(results_list) <- all_genes
clean_results <- results_list[sapply(results_list, function(x) !is.null(x) && x$status == "success")]

final_df <- as.data.frame(do.call(rbind, lapply(clean_results, `[[`, "result")))
final_df$Gene <- rownames(final_df)
final_df[, 1:5] <- lapply(final_df[, 1:5], as.numeric)
colnames(final_df)[1:5] <- c("Estimate", "Std_Error", "df", "t_value", "P_Value")

final_df$FDR <- p.adjust(final_df$P_Value, method = "fdr")
final_df <- final_df[, c("Gene", "Estimate", "Std_Error", "t_value", "P_Value", "FDR")]
final_df <- final_df[order(final_df$P_Value), ]

# --- 7. Save ---
write.csv(final_df, "/Users/nikamikhailava/Desktop/Fibroblasts_signature_derivation/fibroblast_DEG_tsi_stroma_corrected_for_chemo_and_sample.csv", row.names = FALSE)