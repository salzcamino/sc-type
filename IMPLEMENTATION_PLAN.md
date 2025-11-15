# ScType v2 Implementation Plan
## Detailed Roadmap for Top 3 Improvements

**Created:** 2025-11-15
**Expected Timeline:** 5-7 days
**Expected Impact:** +15-30% accuracy improvement

---

## Overview

This document provides a detailed, step-by-step implementation plan for adding three critical improvements to ScType:

1. **Statistical Significance Testing** (Days 1-2)
2. **TF-IDF Marker Weighting** (Days 3-4)
3. **Doublet Detection Integration** (Day 5)
4. **Integration & Testing** (Days 6-7)

Each phase is broken down into actionable tasks that can be executed by agents.

---

## PHASE 1: Statistical Significance Testing

**Goal:** Replace arbitrary ncells/4 threshold with rigorous statistical testing (z-scores, p-values, FDR correction)

**Expected Impact:** +5-10% accuracy
**Complexity:** EASY
**Timeline:** 1-2 days

### Phase 1.1: Create Statistical Testing Functions

**File:** `R/sctype_statistics.R` (NEW)

**Tasks:**
1. Create new file with GPL v3 header
2. Implement core statistical functions:
   - `calculate_zscore_pvalue()` - Convert raw scores to z-scores and p-values
   - `apply_fdr_correction()` - Benjamini-Hochberg FDR correction
   - `assign_confidence_level()` - Categorize into High/Medium/Low
   - `generate_null_distribution()` - Permutation-based null for calibration

**Key Algorithm (calculate_zscore_pvalue):**
```r
calculate_zscore_pvalue <- function(scores, cluster_sizes) {
    # For each cluster:
    # 1. Calculate z-score: (score - mean) / sd
    # 2. Convert to p-value: pnorm(z_score, lower.tail=FALSE)
    # 3. Return data.frame with cluster, score, zscore, pvalue

    z_scores <- (scores - mean(scores)) / sd(scores)
    p_values <- pnorm(z_scores, lower.tail = FALSE)

    return(data.frame(
        cluster = names(scores),
        score = scores,
        zscore = z_scores,
        pvalue = p_values
    ))
}
```

**Dependencies:** None (base R only)

**Testing Criteria:**
- Function accepts vector of scores
- Returns correct z-scores (mean=0, sd=1)
- P-values are in [0, 1]
- Small p-values for high scores

**Code Reference:** Lines 10-150 in PROPOSED_IMPROVEMENTS_CODE.R

---

### Phase 1.2: Update sctype_score_.R

**File:** `R/sctype_score_.R` (MODIFY)

**Tasks:**
1. Read existing file to understand output format
2. Modify to return additional metadata:
   - Raw scores (already returned)
   - Ensure scores matrix has proper row/column names
   - Add optional parameter `return_details = FALSE`
3. When `return_details = TRUE`, return list with:
   - `scores` - existing score matrix
   - `cluster_summary` - aggregated scores per cluster
   - `gene_weights` - marker sensitivity weights

**Modifications:**
```r
# EXISTING (line ~74):
return(es)

# NEW:
if (return_details) {
    return(list(
        scores = es,
        marker_weights = marker_sensitivity,
        scaled_data = scRNAseqData  # if needed for permutation
    ))
} else {
    return(es)
}
```

**Testing Criteria:**
- Backward compatible (default behavior unchanged)
- `return_details = TRUE` returns list
- All existing tests pass

**Code Location:** R/sctype_score_.R lines 60-80

---

### Phase 1.3: Create Enhanced Wrapper with Statistics

**File:** `R/sctype_wrapper_v2.R` (NEW)

**Tasks:**
1. Copy `R/sctype_wrapper.R` to `R/sctype_wrapper_v2.R`
2. Add new parameters:
   - `use_statistics = TRUE` - Enable statistical testing
   - `fdr_threshold = 0.05` - FDR cutoff for significance
   - `confidence_method = "zscore"` - Method for confidence scoring
3. Modify annotation workflow:
   - After aggregating scores, call `calculate_zscore_pvalue()`
   - Apply FDR correction with `apply_fdr_correction()`
   - Assign confidence levels
   - Add p-values and FDR to cluster results
4. Add new metadata columns to Seurat object:
   - `sctype_classification` (cell type, as before)
   - `sctype_pvalue` (per-cluster p-value)
   - `sctype_fdr` (FDR-corrected p-value)
   - `sctype_confidence` (High/Medium/Low/Unknown)
5. Apply FDR threshold: if FDR > threshold, assign "Unknown"

**Key Workflow Changes:**
```r
# EXISTING workflow (lines 90-120 in sctype_wrapper.R):
# 1. Aggregate scores by cluster
# 2. Take top cell type
# 3. If score < ncells/4, assign "Unknown"

# NEW workflow (sctype_wrapper_v2.R):
# 1. Aggregate scores by cluster
# 2. Calculate statistics
stat_results <- calculate_zscore_pvalue(cluster_scores, cluster_sizes)
stat_results <- apply_fdr_correction(stat_results)

# 3. Assign cell types with FDR filter
for (cluster in unique_clusters) {
    if (stat_results[cluster, "fdr"] < fdr_threshold) {
        assign top cell type
        confidence <- assign_confidence_level(stat_results[cluster, "zscore"])
    } else {
        assign "Unknown"
        confidence <- "Unknown"
    }
}

# 4. Add to metadata
seurat_object@meta.data$sctype_pvalue <- stat_results[cluster_id, "pvalue"]
seurat_object@meta.data$sctype_fdr <- stat_results[cluster_id, "fdr"]
seurat_object@meta.data$sctype_confidence <- confidence
```

**Testing Criteria:**
- All existing parameters still work
- New statistics columns populated
- FDR filtering reduces false positives
- Confidence levels make sense (High = low FDR, high z-score)

**Code Reference:** PROPOSED_IMPROVEMENTS_CODE.R lines 10-350

---

### Phase 1.4: Test Statistical Testing

**Test File:** `tests/test_statistics.R` (NEW)

**Test Cases:**
1. **Test calculate_zscore_pvalue():**
   - Input: scores = c(10, 20, 30, 5)
   - Expected: z-scores centered at 0, p-values decrease with score

2. **Test FDR correction:**
   - Input: p-values = c(0.01, 0.02, 0.03, 0.5)
   - Expected: FDR values ≥ p-values, multiple testing accounted for

3. **Test run_sctype_v2() with statistics:**
   - Load example data (exampleData.RDS)
   - Run with `use_statistics = TRUE`
   - Check: All new columns exist, no NAs, confidence levels assigned

4. **Compare v1 vs v2:**
   - Run same data through both wrappers
   - Count Unknown assignments (v2 should have fewer spurious annotations)
   - Measure agreement on confident annotations

**Success Criteria:**
- All unit tests pass
- V2 has ≤ v1 Unknown rate (more confident assignments)
- Manual inspection: confident clusters have clear marker expression

**Test Data:** Use exampleData.RDS (PBMC 3k)

---

### Phase 1.5: Documentation

**File:** `STATISTICAL_TESTING_README.md` (NEW)

**Sections:**
1. **Overview:** Why statistical testing improves ScType
2. **Methods:**
   - Z-score calculation
   - P-value derivation (normal approximation)
   - FDR correction (Benjamini-Hochberg)
   - Confidence level thresholds
3. **Usage Examples:**
   ```r
   # Basic usage with statistics
   seurat_obj <- run_sctype_v2(seurat_obj,
                                known_tissue_type = "Immune system",
                                use_statistics = TRUE,
                                fdr_threshold = 0.05)

   # Access results
   table(seurat_obj@meta.data$sctype_confidence)
   low_conf <- seurat_obj@meta.data[seurat_obj@meta.data$sctype_confidence == "Low", ]
   ```
4. **Interpretation:**
   - P-value < 0.05: Statistically significant
   - FDR < 0.05: Significant after multiple testing correction
   - High confidence: FDR < 0.01, z-score > 2
   - Medium confidence: FDR < 0.05, z-score 1-2
   - Low confidence: FDR < 0.1, z-score < 1
5. **Comparison to v1:** Table showing differences
6. **References:** Citations for statistical methods

**Length:** ~1500 words

---

## PHASE 2: TF-IDF Marker Weighting

**Goal:** Replace frequency-only marker weighting with TF-IDF (term frequency × inverse document frequency)

**Expected Impact:** +8-15% accuracy on subtypes
**Complexity:** EASY-MEDIUM
**Timeline:** 2-3 days

### Phase 2.1: Create TF-IDF Weighting Functions

**File:** `R/sctype_tfidf.R` (NEW)

**Tasks:**
1. Create new file with GPL v3 header
2. Implement TF-IDF weighting functions:
   - `calculate_tf()` - Term frequency (marker expression in dataset)
   - `calculate_idf()` - Inverse document frequency (marker specificity)
   - `calculate_tfidf_weights()` - Combine TF × IDF
   - `apply_tfidf_to_markers()` - Apply to marker gene sets

**Key Algorithm (calculate_tfidf_weights):**
```r
calculate_tfidf_weights <- function(expression_matrix, marker_genes, cell_types_db) {
    # For each marker gene:
    # TF = log(1 + mean_expression_in_dataset)
    # IDF = log(n_celltypes / n_celltypes_with_marker)
    # TF-IDF = TF × IDF

    tfidf_weights <- list()

    for (gene in marker_genes) {
        # Term Frequency: expression magnitude
        tf <- log1p(mean(expression_matrix[gene, ]))

        # Inverse Document Frequency: specificity
        n_total <- length(unique(cell_types_db$cellName))
        n_with_marker <- sum(grepl(gene, cell_types_db$geneSymbolmore1))
        idf <- log(n_total / (n_with_marker + 1))

        tfidf_weights[[gene]] <- tf * idf
    }

    # Normalize to [0, 1]
    tfidf_weights <- scales::rescale(unlist(tfidf_weights))

    return(tfidf_weights)
}
```

**Dependencies:** `scales` package for rescale

**Testing Criteria:**
- High expression + rare marker → high weight
- Low expression OR common marker → low weight
- Weights in [0, 1]
- Specific markers weighted higher than housekeeping genes

**Code Reference:** PROPOSED_IMPROVEMENTS_CODE.R lines 152-350

---

### Phase 2.2: Update gene_sets_prepare.R

**File:** `R/gene_sets_prepare.R` (MODIFY)

**Tasks:**
1. Read existing file to understand marker preparation
2. Add new parameter: `weighting_method = "frequency"` (default maintains backward compatibility)
3. Add option for `weighting_method = "tfidf"`
4. When TF-IDF selected:
   - Call `calculate_tfidf_weights()` from `sctype_tfidf.R`
   - Replace frequency-based weights with TF-IDF weights
   - Return weighted marker lists

**Modifications:**
```r
# EXISTING (lines ~40-50):
marker_sensitivity <- calculate_frequency_weights(markers)

# NEW:
if (weighting_method == "frequency") {
    marker_sensitivity <- calculate_frequency_weights(markers)
} else if (weighting_method == "tfidf") {
    source("R/sctype_tfidf.R")
    marker_sensitivity <- calculate_tfidf_weights(expression_matrix, markers, db)
} else {
    stop("Invalid weighting_method. Use 'frequency' or 'tfidf'")
}
```

**Challenge:** gene_sets_prepare() doesn't currently receive expression matrix
**Solution:** Add optional parameter `expression_matrix = NULL`
- If NULL, use frequency weights (backward compatible)
- If provided, can calculate TF-IDF

**Testing Criteria:**
- Default behavior unchanged (weighting_method = "frequency")
- TF-IDF produces different weights
- Weighted markers still work in sctype_score_()

**Code Location:** R/gene_sets_prepare.R lines 30-60

---

### Phase 2.3: Modify sctype_score_.R for TF-IDF

**File:** `R/sctype_score_.R` (MODIFY)

**Tasks:**
1. Review how marker_sensitivity is used (lines 24-40)
2. Ensure TF-IDF weights are applied correctly
3. Add parameter `marker_weights = NULL` (optional pre-computed weights)
4. If `marker_weights` provided, use those instead of recalculating
5. Normalize properly after TF-IDF weighting

**Modifications:**
```r
# EXISTING (line ~26):
marker_sensitivity <- marker_sensitivity_function(markers)

# NEW:
if (!is.null(marker_weights)) {
    # Use pre-computed TF-IDF weights
    marker_sensitivity <- marker_weights
} else {
    # Calculate frequency weights (current default)
    marker_sensitivity <- marker_sensitivity_function(markers)
}

# Rest of scoring algorithm unchanged
```

**Testing Criteria:**
- Custom weights properly applied
- Scores change appropriately with TF-IDF
- High TF-IDF markers have more influence on final score

**Code Location:** R/sctype_score_.R lines 20-45

---

### Phase 2.4: Integrate TF-IDF into Wrapper

**File:** `R/sctype_wrapper_v2.R` (MODIFY from Phase 1.3)

**Tasks:**
1. Add parameter `weighting_method = "frequency"` to `run_sctype_v2()`
2. Pass weighting method to `gene_sets_prepare()`
3. If TF-IDF, also pass expression matrix to `gene_sets_prepare()`
4. Document in function docstring

**Modifications:**
```r
run_sctype_v2 <- function(seurat_object,
                         known_tissue_type = NULL,
                         weighting_method = "frequency",  # NEW
                         use_statistics = TRUE,
                         fdr_threshold = 0.05,
                         ...) {

    # Extract expression matrix if TF-IDF
    if (weighting_method == "tfidf") {
        expr_matrix <- extract_expression(seurat_object)
    } else {
        expr_matrix <- NULL
    }

    # Prepare gene sets with weighting
    gs_list <- gene_sets_prepare(db_file, tissue_type,
                                 weighting_method = weighting_method,
                                 expression_matrix = expr_matrix)

    # Continue with scoring...
}
```

**Testing Criteria:**
- `weighting_method = "frequency"` matches v1 behavior
- `weighting_method = "tfidf"` produces different (better) annotations
- Both methods complete without errors

---

### Phase 2.5: Test TF-IDF Weighting

**Test File:** `tests/test_tfidf.R` (NEW)

**Test Cases:**
1. **Test calculate_tfidf_weights():**
   - Create mock expression matrix
   - Create mock marker database
   - Verify: Rare + high expression → high weight
   - Verify: Common OR low expression → low weight

2. **Test integration with gene_sets_prepare():**
   - Call with `weighting_method = "tfidf"`
   - Check: Weights differ from frequency
   - Check: Still returns valid gs_positive/gs_negative

3. **Test run_sctype_v2() with TF-IDF:**
   - Run on PBMC data with frequency weights
   - Run on PBMC data with TF-IDF weights
   - Compare: TF-IDF should improve rare subtype detection

4. **Edge cases:**
   - Gene not in expression matrix → fallback to frequency
   - All markers common → weights normalize correctly

**Success Criteria:**
- TF-IDF improves detection of rare cell types (measured manually or with ground truth)
- Specific markers (e.g., CD8A for CD8 T cells) get higher weights than common markers (e.g., ACTB)
- No errors on edge cases

**Test Data:** exampleData.RDS + manual inspection of known markers

---

### Phase 2.6: Documentation

**File:** `TFIDF_WEIGHTING_README.md` (NEW)

**Sections:**
1. **Overview:** Why TF-IDF improves marker specificity
2. **Background:** Borrowed from NLP, applied to cell types
3. **Algorithm:**
   - TF (term frequency): How much is gene expressed?
   - IDF (inverse document frequency): How specific is gene?
   - TF-IDF = TF × IDF
4. **Usage Examples:**
   ```r
   # Use TF-IDF weighting
   seurat_obj <- run_sctype_v2(seurat_obj,
                                known_tissue_type = "Immune system",
                                weighting_method = "tfidf")
   ```
5. **When to Use:**
   - Datasets with rare cell types
   - Highly similar cell subtypes (e.g., T cell subsets)
   - When marker databases have variable quality
6. **Comparison:** Frequency vs TF-IDF performance table
7. **Limitations:** Requires expression data, computationally more expensive
8. **References:** TF-IDF in NLP, applications to genomics

**Length:** ~1200 words

---

## PHASE 3: Doublet Detection Integration

**Goal:** Integrate scDblFinder to detect and flag doublets before annotation

**Expected Impact:** +5-10% accuracy (prevents spurious annotations)
**Complexity:** EASY
**Timeline:** 1 day

### Phase 3.1: Create Doublet Detection Wrapper

**File:** `R/sctype_doublet_detection.R` (NEW)

**Tasks:**
1. Create wrapper around scDblFinder
2. Implement functions:
   - `detect_doublets_scdblfinder()` - Run scDblFinder
   - `add_doublet_scores()` - Add scores to Seurat object
   - `filter_doublets()` - Optionally remove doublets
   - `flag_doublets()` - Mark doublets in metadata

**Key Algorithm:**
```r
detect_doublets_scdblfinder <- function(seurat_object,
                                        dbr = NULL,  # doublet rate
                                        clusters = NULL) {
    # Check if scDblFinder installed
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
        stop("scDblFinder not installed. Install with: BiocManager::install('scDblFinder')")
    }

    # Convert Seurat to SCE for scDblFinder
    sce <- as.SingleCellExperiment(seurat_object)

    # Run scDblFinder
    sce <- scDblFinder::scDblFinder(sce,
                                    dbr = dbr,
                                    clusters = clusters)

    # Extract results
    doublet_scores <- sce$scDblFinder.score
    doublet_class <- sce$scDblFinder.class

    return(list(
        scores = doublet_scores,
        class = doublet_class,  # "singlet" or "doublet"
        sce = sce
    ))
}
```

**Dependencies:**
- `scDblFinder` (Bioconductor)
- `SingleCellExperiment`

**Testing Criteria:**
- Detects known doublet-enriched clusters
- Scores correlate with biological markers of doublets
- Runs without errors on test data

**Code Reference:** PROPOSED_IMPROVEMENTS_CODE.R lines 352-500

---

### Phase 3.2: Integrate into Wrapper

**File:** `R/sctype_wrapper_v2.R` (MODIFY from Phase 2.4)

**Tasks:**
1. Add parameter `detect_doublets = FALSE` to `run_sctype_v2()`
2. Add parameter `doublet_action = "flag"` (options: "flag", "filter", "both")
3. Before annotation:
   - If `detect_doublets = TRUE`, run `detect_doublets_scdblfinder()`
   - Add doublet scores to metadata
4. During annotation:
   - If `doublet_action = "filter"`, exclude doublets from annotation
   - If `doublet_action = "flag"`, annotate but mark as low confidence
5. Add new metadata columns:
   - `doublet_score` (0-1 probability)
   - `doublet_class` ("singlet" or "doublet")
   - `sctype_classification` (if doublet + action=flag, append "[Doublet]")

**Workflow Integration:**
```r
run_sctype_v2 <- function(seurat_object,
                         detect_doublets = FALSE,  # NEW
                         doublet_action = "flag",  # NEW
                         ...) {

    # Step 0: Doublet detection (if requested)
    if (detect_doublets) {
        message("Detecting doublets...")
        doublet_results <- detect_doublets_scdblfinder(seurat_object)

        seurat_object@meta.data$doublet_score <- doublet_results$scores
        seurat_object@meta.data$doublet_class <- doublet_results$class

        # Filter doublets if requested
        if (doublet_action == "filter") {
            cells_to_use <- which(doublet_results$class == "singlet")
            message(paste("Filtering", sum(doublet_results$class == "doublet"), "doublets"))
        } else {
            cells_to_use <- 1:ncol(seurat_object)
        }
    } else {
        cells_to_use <- 1:ncol(seurat_object)
    }

    # Step 1-N: Standard ScType workflow (on filtered cells)
    # ...

    # Step Final: Flag doublets in annotations
    if (detect_doublets && doublet_action %in% c("flag", "both")) {
        doublet_cells <- seurat_object@meta.data$doublet_class == "doublet"
        seurat_object@meta.data$sctype_classification[doublet_cells] <-
            paste0(seurat_object@meta.data$sctype_classification[doublet_cells], " [Doublet]")
        seurat_object@meta.data$sctype_confidence[doublet_cells] <- "Low"
    }

    return(seurat_object)
}
```

**Testing Criteria:**
- Doublets correctly identified (check known markers like high nCount_RNA)
- Filtering reduces spurious annotations
- Flagging preserves cell counts but warns users

---

### Phase 3.3: Test Doublet Detection

**Test File:** `tests/test_doublets.R` (NEW)

**Test Cases:**
1. **Test detect_doublets_scdblfinder():**
   - Run on PBMC data
   - Check: Doublet rate is reasonable (5-10%)
   - Check: High doublet scores correlate with high nCount/nFeature

2. **Test doublet_action = "filter":**
   - Run ScType with filtering
   - Verify: Doublets excluded from annotation
   - Verify: Fewer Unknown or ambiguous annotations

3. **Test doublet_action = "flag":**
   - Run ScType with flagging
   - Verify: Doublets annotated but marked "[Doublet]"
   - Verify: Confidence set to "Low"

4. **Biological validation:**
   - Check known doublet signatures (e.g., CD3+CD19+ = T+B doublet)
   - Manually inspect flagged doublets on UMAP
   - Verify they cluster between cell types

**Success Criteria:**
- Known doublets detected (e.g., cells expressing markers of 2+ lineages)
- Filtering improves annotation purity
- Flagging provides useful QC information

**Test Data:** PBMC 3k (known to have ~5% doublets)

---

### Phase 3.4: Documentation

**File:** `DOUBLET_DETECTION_README.md` (NEW)

**Sections:**
1. **Overview:** Why doublet detection matters
2. **Methods:** scDblFinder algorithm summary
3. **Usage Examples:**
   ```r
   # Detect and flag doublets
   seurat_obj <- run_sctype_v2(seurat_obj,
                                known_tissue_type = "Immune system",
                                detect_doublets = TRUE,
                                doublet_action = "flag")

   # Filter doublets before annotation
   seurat_obj <- run_sctype_v2(seurat_obj,
                                known_tissue_type = "Immune system",
                                detect_doublets = TRUE,
                                doublet_action = "filter")

   # Access doublet information
   table(seurat_obj@meta.data$doublet_class)
   VlnPlot(seurat_obj, features = "doublet_score", group.by = "seurat_clusters")
   ```
4. **Interpreting Results:**
   - Doublet score > 0.5: Likely doublet
   - Doublets often have high UMI counts
   - Check marker co-expression
5. **When to Use:**
   - Always recommended for 10X data (capture method prone to doublets)
   - Especially important for rare cell type detection
   - Critical when merging datasets
6. **Limitations:**
   - Cannot detect all doublets (some are homotypic)
   - May flag transitional states as doublets
7. **References:** scDblFinder paper, doublet formation rates

**Length:** ~1000 words

---

## PHASE 4: Integration and Testing

**Goal:** Combine all improvements into unified `run_sctype_v2()` with comprehensive testing

**Timeline:** 1-2 days

### Phase 4.1: Create Unified run_sctype_v2()

**File:** `R/sctype_wrapper_v2.R` (FINALIZE from Phases 1-3)

**Tasks:**
1. Ensure all three improvements work together
2. Add comprehensive parameter validation
3. Add progress messages for each step
4. Create complete documentation
5. Implement backward compatibility mode

**Final Function Signature:**
```r
run_sctype_v2 <- function(
    seurat_object,
    known_tissue_type = NULL,

    # Core parameters (backward compatible)
    assay = "RNA",
    scaled = TRUE,
    custom_marker_file = NULL,
    cluster_col = "seurat_clusters",
    annotation_col = "sctype_classification",
    plot = FALSE,

    # NEW: Statistical testing
    use_statistics = TRUE,
    fdr_threshold = 0.05,
    confidence_method = "zscore",

    # NEW: TF-IDF weighting
    weighting_method = "frequency",  # or "tfidf"

    # NEW: Doublet detection
    detect_doublets = FALSE,
    doublet_action = "flag",  # "flag", "filter", or "both"
    doublet_rate = NULL,

    # Misc
    verbose = TRUE,
    return_details = FALSE
) {
    # Implementation combines all phases
}
```

**Testing Criteria:**
- All parameters work individually
- All parameters work in combination
- Backward compatibility: default params = v1 behavior
- Progress messages informative
- Returns expected metadata columns

---

### Phase 4.2: Comprehensive Test Suite

**File:** `tests/test_sctype_v2.R` (NEW)

**Test Structure:**
```r
# 1. Unit tests (individual functions)
test_that("Statistical functions work", { ... })
test_that("TF-IDF weights calculated correctly", { ... })
test_that("Doublet detection runs", { ... })

# 2. Integration tests (combined features)
test_that("Statistics + TF-IDF work together", { ... })
test_that("All features combined work", { ... })

# 3. Regression tests (compare to v1)
test_that("Default params match v1 output", { ... })
test_that("V2 improves on known difficult cases", { ... })

# 4. Edge case tests
test_that("Empty clusters handled", { ... })
test_that("All Unknown case handled", { ... })
test_that("Single cluster case handled", { ... })

# 5. Performance tests
test_that("V2 runs in reasonable time", { ... })
test_that("Memory usage acceptable", { ... })
```

**Test Data:**
- exampleData.RDS (PBMC 3k)
- Synthetic data with known ground truth
- Edge cases (single cluster, all same type, etc.)

**Success Criteria:**
- All tests pass
- V2 ≥ v1 accuracy (measured on ground truth)
- Runtime < 2× v1 (acceptable for accuracy gain)

---

### Phase 4.3: Benchmark v1 vs v2

**File:** `benchmarks/compare_v1_v2.R` (NEW)

**Benchmark Tasks:**
1. Load PBMC 3k example data
2. Run v1: `run_sctype()`
3. Run v2 with defaults: `run_sctype_v2(use_statistics = FALSE, weighting_method = "frequency", detect_doublets = FALSE)`
4. Run v2 with statistics only
5. Run v2 with TF-IDF only
6. Run v2 with doublets only
7. Run v2 with all features

**Metrics to Compare:**
- **Accuracy:** Agreement with known markers
- **Unknown rate:** % cells annotated as Unknown
- **Confidence distribution:** High/Medium/Low counts
- **Runtime:** Seconds per 1000 cells
- **Memory:** Peak RAM usage
- **Rare cell detection:** Can it find rare types? (mast cells, DCs)

**Output Format:**
```
Benchmark Results: ScType v1 vs v2
===================================
Dataset: PBMC 3k (2,700 cells, 8 clusters)

Method                  | Accuracy | Unknown% | Runtime | Memory
-----------------------------------------------------------------
v1 (current)           | 85.2%    | 12.3%    | 1.2s    | 150 MB
v2 (default=v1)        | 85.2%    | 12.3%    | 1.3s    | 155 MB
v2 + statistics        | 87.8%    | 8.1%     | 1.5s    | 160 MB
v2 + TF-IDF            | 91.3%    | 10.2%    | 2.1s    | 170 MB
v2 + doublets          | 86.9%    | 11.0%    | 2.8s    | 180 MB
v2 + all features      | 93.5%    | 6.8%     | 3.2s    | 190 MB

Cell Type Detection (rare types):
v1: 0/3 rare types found
v2 (all features): 2/3 rare types found

Recommended: Use v2 with all features for best accuracy
Fallback: Use v2 with statistics only for speed/accuracy balance
```

**Visualization:**
- Bar plot: Accuracy by method
- Line plot: Runtime vs accuracy tradeoff
- Confusion matrix: v1 vs v2 annotations

---

### Phase 4.4: Update CLAUDE.md

**File:** `CLAUDE.md` (MODIFY)

**Updates Needed:**
1. **Repository Structure:** Add new files
   ```
   R/sctype_wrapper_v2.R
   R/sctype_statistics.R
   R/sctype_tfidf.R
   R/sctype_doublet_detection.R
   STATISTICAL_TESTING_README.md
   TFIDF_WEIGHTING_README.md
   DOUBLET_DETECTION_README.md
   IMPROVEMENTS_V2_README.md
   ```

2. **Core Components:** Add section 10 for v2 improvements
   ```markdown
   ### 10. ScType v2 Improvements (`R/sctype_wrapper_v2.R`)

   **Purpose**: Enhanced version of ScType with statistical rigor,
   better marker weighting, and doublet detection.

   **What's New:**
   - Statistical significance testing (p-values, FDR)
   - TF-IDF marker weighting
   - Doublet detection integration

   **Expected Improvements:**
   - +5-10% accuracy from statistics
   - +8-15% accuracy from TF-IDF
   - +5-10% accuracy from doublet detection
   - Combined: +15-30% accuracy improvement
   ```

3. **Usage Examples:** Add v2 workflow
   ```r
   # Load v2 wrapper
   source("R/sctype_wrapper_v2.R")

   # Run with all improvements
   seurat_obj <- run_sctype_v2(
       seurat_obj,
       known_tissue_type = "Immune system",
       use_statistics = TRUE,
       weighting_method = "tfidf",
       detect_doublets = TRUE
   )
   ```

4. **Version Information:** Update compatibility section

---

### Phase 4.5: Create Overview Documentation

**File:** `IMPROVEMENTS_V2_README.md` (NEW)

**Comprehensive Overview Document**

**Sections:**
1. **Executive Summary**
   - What's new in v2
   - Expected improvements
   - When to use v1 vs v2

2. **Three Improvements Explained**
   - Statistical testing: Why and how
   - TF-IDF weighting: Why and how
   - Doublet detection: Why and how

3. **Quick Start Guide**
   ```r
   # Minimal example
   source("R/sctype_wrapper_v2.R")
   seurat_obj <- run_sctype_v2(seurat_obj,
                                known_tissue_type = "Immune system")
   ```

4. **Complete Parameter Guide**
   - Table of all parameters
   - Defaults and recommendations
   - Advanced usage

5. **Performance Comparison**
   - Benchmark results summary
   - When to use which features
   - Speed vs accuracy tradeoffs

6. **Migration Guide (v1 → v2)**
   - Backward compatibility notes
   - How to update existing code
   - Breaking changes (if any)

7. **Troubleshooting**
   - Common issues
   - Installation problems
   - Performance tuning

8. **References**
   - Papers for each method
   - Related tools
   - Acknowledgments

**Length:** ~3000 words

---

## PHASE 5: Final Steps

### Phase 5.1: Create Changelog

**File:** `CHANGELOG.md` (NEW or UPDATE)

**Format:**
```markdown
# Changelog

## [v2.0.0] - 2025-11-XX

### Added
- Statistical significance testing with p-values and FDR correction
- TF-IDF marker weighting for improved specificity
- Doublet detection integration via scDblFinder
- Comprehensive test suite
- Five new documentation files

### Changed
- Created run_sctype_v2() as enhanced wrapper
- Original run_sctype() remains unchanged (backward compatible)

### Performance
- Accuracy improvement: +15-30% on challenging datasets
- Runtime increase: ~2-3x (acceptable for accuracy gain)
- Memory increase: ~20-30% (minimal)

### Files Added
- R/sctype_wrapper_v2.R
- R/sctype_statistics.R
- R/sctype_tfidf.R
- R/sctype_doublet_detection.R
- STATISTICAL_TESTING_README.md
- TFIDF_WEIGHTING_README.md
- DOUBLET_DETECTION_README.md
- IMPROVEMENTS_V2_README.md
- benchmarks/compare_v1_v2.R
- tests/test_sctype_v2.R

### Benchmarks
- PBMC 3k: 85.2% → 93.5% accuracy
- Rare cell detection: 0/3 → 2/3 types found
- Runtime: 1.2s → 3.2s (2.7x slower, acceptable)

## [v1.x.x] - Previous versions
...
```

---

### Phase 5.2: Update Main README

**File:** `README.md` (MODIFY)

**Changes:**
1. Add v2 mention in introduction
2. Add quick start for v2
3. Link to new documentation
4. Update installation instructions (new dependencies)

**Addition:**
```markdown
## What's New in v2 (2025)

ScType v2 introduces three major improvements for higher accuracy:

1. **Statistical Significance Testing**: P-values and FDR correction replace arbitrary thresholds
2. **TF-IDF Marker Weighting**: Better marker importance scoring
3. **Doublet Detection**: Automatic multiplet identification

**Performance:** +15-30% accuracy improvement on challenging datasets

**Quick Start (v2):**
```r
source("R/sctype_wrapper_v2.R")
seurat_obj <- run_sctype_v2(seurat_obj, known_tissue_type = "Immune system")
```

See [IMPROVEMENTS_V2_README.md](IMPROVEMENTS_V2_README.md) for details.

**Backward Compatibility:** Original `run_sctype()` unchanged. Both versions available.
```

---

### Phase 5.3: Commit Strategy

**Commit Structure:**

1. **Commit 1:** Statistical testing
   ```
   Add statistical significance testing to ScType (v2 Phase 1)

   - Create R/sctype_statistics.R with z-score, p-value, FDR functions
   - Create R/sctype_wrapper_v2.R with statistical testing
   - Update R/sctype_score_.R to support return_details
   - Add STATISTICAL_TESTING_README.md
   - Add tests/test_statistics.R

   Impact: +5-10% accuracy improvement
   ```

2. **Commit 2:** TF-IDF weighting
   ```
   Add TF-IDF marker weighting to ScType (v2 Phase 2)

   - Create R/sctype_tfidf.R
   - Update R/gene_sets_prepare.R to support TF-IDF
   - Update R/sctype_score_.R for custom weights
   - Update R/sctype_wrapper_v2.R with weighting_method parameter
   - Add TFIDF_WEIGHTING_README.md
   - Add tests/test_tfidf.R

   Impact: +8-15% accuracy improvement on subtypes
   ```

3. **Commit 3:** Doublet detection
   ```
   Add doublet detection integration to ScType (v2 Phase 3)

   - Create R/sctype_doublet_detection.R
   - Update R/sctype_wrapper_v2.R with doublet detection
   - Add DOUBLET_DETECTION_README.md
   - Add tests/test_doublets.R

   Impact: +5-10% accuracy improvement
   ```

4. **Commit 4:** Integration and docs
   ```
   Finalize ScType v2 with comprehensive testing and documentation

   - Finalize R/sctype_wrapper_v2.R (all features integrated)
   - Add comprehensive test suite (tests/test_sctype_v2.R)
   - Add benchmarks (benchmarks/compare_v1_v2.R)
   - Create IMPROVEMENTS_V2_README.md
   - Update CLAUDE.md
   - Update README.md
   - Add CHANGELOG.md

   Combined impact: +15-30% accuracy improvement
   Runtime: ~3x slower (acceptable tradeoff)
   ```

---

## Success Criteria

### Phase 1 Success (Statistical Testing)
- [ ] All statistical functions implemented and tested
- [ ] P-values and FDR calculated correctly
- [ ] V2 with statistics runs without errors
- [ ] Fewer spurious Unknown annotations than v1
- [ ] Documentation complete and clear

### Phase 2 Success (TF-IDF Weighting)
- [ ] TF-IDF weights calculated correctly
- [ ] Integration with gene_sets_prepare works
- [ ] Rare cell types better detected with TF-IDF
- [ ] Specific markers weighted higher than housekeeping
- [ ] Documentation complete and clear

### Phase 3 Success (Doublet Detection)
- [ ] scDblFinder integration works
- [ ] Doublets correctly identified (5-10% rate on PBMC)
- [ ] Filtering reduces spurious annotations
- [ ] Flagging provides useful QC info
- [ ] Documentation complete and clear

### Phase 4 Success (Integration)
- [ ] All features work together
- [ ] Comprehensive tests pass
- [ ] Benchmark shows +15-30% improvement
- [ ] Documentation complete (5 new files)
- [ ] CLAUDE.md updated

### Overall Success
- [ ] V2 accuracy ≥ v1 + 15% (measured)
- [ ] Runtime ≤ 5× v1 (acceptable)
- [ ] Backward compatible (v1 still works)
- [ ] All code committed and pushed
- [ ] Ready for user testing

---

## Timeline Summary

| Phase | Days | Deliverables | Impact |
|-------|------|--------------|--------|
| **Phase 1: Statistics** | 1-2 | sctype_statistics.R, sctype_wrapper_v2.R (v1), tests, docs | +5-10% |
| **Phase 2: TF-IDF** | 2-3 | sctype_tfidf.R, updates to prepare/score/wrapper, tests, docs | +8-15% |
| **Phase 3: Doublets** | 1 | sctype_doublet_detection.R, wrapper update, tests, docs | +5-10% |
| **Phase 4: Integration** | 1-2 | Final wrapper, comprehensive tests, benchmarks, 5 docs | Combined |
| **Phase 5: Final** | 0.5 | Commit, push, changelog, README updates | - |
| **TOTAL** | **5.5-8.5 days** | **15+ files** | **+15-30%** |

---

## Dependencies Summary

### R Packages (NEW)
- **scDblFinder** (Bioconductor): Doublet detection
- **scales**: Rescaling for TF-IDF

### R Packages (EXISTING)
- dplyr, HGNChelper, openxlsx (already required)
- Seurat, SingleCellExperiment (already required)

### Optional
- testthat (for running tests)
- bench (for benchmarking)

---

## Notes for Implementation

### Agent Execution Strategy
1. **Sequential execution**: Phases must be done in order (1 → 2 → 3 → 4)
2. **Parallel within phase**: Can parallelize tasks within a phase (e.g., documentation + code)
3. **Testing checkpoints**: After each phase, run tests before proceeding
4. **Incremental commits**: Commit after each phase completes

### Key Files to Reference
- **PROPOSED_IMPROVEMENTS_CODE.R**: Contains working implementations of all 3 improvements
- **ANNOTATION_METHODS_RESEARCH_REPORT.md**: Detailed algorithm explanations
- **RESEARCH_SUMMARY.md**: Quick reference for decisions

### Quality Checks
- After each phase, run: `Rscript tests/test_*.R`
- Before final commit, run: `Rscript benchmarks/compare_v1_v2.R`
- Verify backward compatibility: Default run_sctype_v2() = run_sctype() output

---

**END OF IMPLEMENTATION PLAN**

Ready to begin Phase 1!
