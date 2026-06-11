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

raw_merged <- raw_expression %>%
  left_join(metadata, by = "dcc_filename")

# --- 3. Filter for Stroma, TSI, AND ONLY POST-NACT ---
raw_merged <- raw_merged %>%
  filter(Segment == "stroma", 
         Segment_geomx %in% c("stroma", "tsi"),
         NACT_status == "post") # Added filter for post-samples only

raw_merged$Segment_geomx <- factor(
  raw_merged$Segment_geomx,
  levels = c("stroma", "tsi")
)

# --- 4. Define LMM Function ---
run_gene_test_satterthwaite_post <- function(gene, df) {
  tryCatch({
    y_values <- as.numeric(trimws(as.character(df[[gene]])))
    
    if (sum(!is.na(y_values)) < 3) {
      return(list(gene = gene, status = "error", message = "Not enough non-NA values", result = NA))
    }
    if (var(y_values, na.rm = TRUE) == 0) {
      return(list(gene = gene, status = "error", message = "Zero variance", result = NA))
    }
    
    # Fit model WITHOUT NACT_status (since all samples are now "post")
    mod <- lmer(y_values ~ Segment_geomx + (1 | Sample), data = df)
    summary_mod <- summary(mod)$coefficients
    
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
cat("Starting Post-NACT only LMM analysis for", length(all_genes), "genes...\n")
results_list <- pblapply(all_genes, run_gene_test_satterthwaite_post, df = raw_merged)

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
# Updated file name to reflect this is the post-only subset
write.csv(final_df, "/Users/nikamikhailava/Desktop/Fibroblasts_signature_derivation/fibroblast_tsi_stroma_POST_ONLY.csv", row.names = FALSE)