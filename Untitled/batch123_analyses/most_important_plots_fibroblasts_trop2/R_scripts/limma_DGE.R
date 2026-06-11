library(GeomxTools)
library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(BiocParallel)
library(NanoStringNCTools)
library(lme4)
library(pbapply)
library(lmerTest)

raw_expression   <- fread("/Volumes/h30492/farkkilab2/4_CellCycle/geomx/RNA_Expression_for_Analysis/batch123_gene_expression_not_deconvoluted_harmony_batch_corr_deseq2_vst_3003.csv")

metadata <- fread("/Volumes/h30492/farkkilab2/9_EyeMT/Data/geomx/batch123/metadata/dcc_metadata_batch123_cleaned.csv")

aoi_summary <- fread("/Volumes/h30492/farkkilab2/4_CellCycle/geomx/Data/proliferative_aois.csv")

# Fix column name
setnames(metadata, "dcc_filename", "dcc_filename_1")

transpose_expr <- function(df) {
  df_t <- as.data.frame(t(df))
  colnames(df_t) <- df_t[1, ]
  df_t <- df_t[-1, ]
  df_t$dcc_filename_1 <- rownames(df_t)
  rownames(df_t) <- NULL
  return(df_t)
}

raw_expr_t   <- transpose_expr(raw_expression)


raw_merged <- aoi_summary %>%
  inner_join(raw_expr_t, by = "dcc_filename_1") %>%
  left_join(metadata %>% select(dcc_filename_1, Sample), by = "dcc_filename_1")

raw_merged <- raw_merged %>%
  filter(AOI_category %in% c("Low P / Low C", "High P / High C"))

raw_merged$AOI_category <- factor(
  raw_merged$AOI_category,
  levels = c("Low P / Low C", "High P / High C")
)

vars_to_exclude <- c("dcc_filename_1", "Sample", "AOI_category", "Tumor_total_count", 
                     "Tumor_Ki67_count", "Tumor_Ki67_pRb_count", "Tumor_PH3_count", 
                     "Ki67_pRb_niche_count", "Ki67_pRb_niche_prop", "Ki67_prop", "Ki67_pRb_prop", "PH3_prop", "V1")

all_genes <- setdiff(colnames(raw_merged), vars_to_exclude)

############## function for the lmm running #################
run_gene_test_satterthwaite <- function(gene, df) {
  tryCatch({
    y_values <- as.numeric(trimws(as.character(df[[gene]])))
    # remove the genes where the number of available values is 0
    if (sum(!is.na(y_values)) < 3) {
      return(list(
        gene = gene,
        status = "error",
        message = "Not enough non-NA values",
        result = NA
      ))
    }
    #### removing the genes where the variance is 0
    if (var(y_values, na.rm = TRUE) == 0) {
      return(list(
        gene = gene,
        status = "error",
        message = "Zero variance",
        result = NA
      ))
    }

    # Fit model
    mod <- lmer(y_values ~ AOI_category + (1 | Sample), data = df)
    summary_mod <- summary(mod)$coefficients
    
    # Extract contrast
    target_row <- grep("AOI_categoryHigh P / High C", rownames(summary_mod), value = TRUE)
    
    if (length(target_row) == 0) {
      return(list(
        gene = gene,
        status = "error",
        message = "Contrast not found",
        result = NA
      ))
    }
    
    stats <- summary_mod[target_row, ]
    
    return(list(
      gene = gene,
      status = "success",
      message =NA ,
      result = c(
        Estimate = stats["Estimate"],
        Std_Error = stats["Std. Error"],
        df = stats["df"],
        t_value = stats["t value"],
        P_Value = stats["Pr(>|t|)"]
      )
    ))

  }, error = function(e) {
    return(list(
      gene = gene,
      status = "error",
      message = e$message,
      result = NA
    ))
  })
}

# --- 5. Run Analysis ---
cat("Starting LMM analysis for", length(all_genes), "genes...\n")


# Using pbapply for progress bar
results_list <- pblapply(all_genes, run_gene_test_satterthwaite, df = raw_merged)


# 1. Assign gene names to the list (ensure all_genes matches results_list length)
names(results_list) <- all_genes

# 2. Remove NULLs and genes that returned an 'error' status
# This ensures we only try to bind successful models
clean_results <- results_list[sapply(results_list, function(x) !is.null(x) && x$status == "success")]

# 3. Extract the 'result' vector from each entry and stack them
matrix_data <- do.call(rbind, lapply(clean_results, `[[`, "result"))

# 4. Convert to a clean Data Frame
final_df <- as.data.frame(matrix_data)

# 5. Add the Gene names as a proper column (instead of just row names)
final_df$Gene <- rownames(final_df)

# 6. Ensure all statistical columns are numeric (rbind sometimes makes them 'any')
final_df[, 1:5] <- lapply(final_df[, 1:5], as.numeric)

# 7. Add FDR correction (Essential for 14,000+ genes!)

colnames(final_df) <- c("Estimate", "Std_Error", "df", "t_value", "P_Value", "Gene")

final_df$FDR <- p.adjust(final_df$P_Value, method = "fdr")
# 8. Reorder columns to put Gene first and sort by P-value
final_df <- final_df[, c("Gene", "Estimate", "Std_Error", "t_value", "P_Value", "FDR")]
final_df <- final_df[order(final_df$P_Value), ]

# --- 7. Save ---
write.csv(final_df, "/Volumes/h30492/farkkilab2/4_CellCycle/geomx/GeoMx_DGE_LMM_Satterthwaite_Results_Bulk.csv", row.names = FALSE)

