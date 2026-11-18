#' Arguments shared by the llomics payload builders
#'
#' Descriptions of arguments that several llomics functions share. Each function's own help page
#' lists the arguments it takes.
#'
#' @param analysis_id Short identifier for the analysis. Stored in the run card and used in file
#'   names.
#' @param question The scientific question, in one sentence. It is sent to the model with every
#'   unit.
#' @param tool Name and version of the tool that ran the test.
#' @param steps Character vector describing the steps that produced the result, such as
#'   normalisation, model and filters.
#' @param observation_type What one observation is, for example `"sample"` or `"cell"`.
#' @param n_features Features loaded per cluster, or per direction for a contrast. Lists are never
#'   padded with ineligible features. `NULL` uses the contract default (15).
#' @param fdr Adjusted p-value below which a feature is eligible for the top list. `NULL` uses the
#'   contract default (0.05).
#' @param min_effect Absolute effect size a feature must exceed to be eligible for the top list.
#' @param effect_cutoff Effect size for the second set of hit counts (`n_up_effect`,
#'   `n_down_effect`). `NULL` uses the contract default (1).
#' @param exclude Regular expression. Features whose readable name matches are left out of the
#'   top list, not out of the test. `NULL` keeps all features.
#' @param panel Features reported for every unit whether or not they rank, matched on feature ID
#'   or symbol. A feature that is not in the result is written as `null`.
#' @param filters List of filters applied before the test, each `list(rule = , why = )`. They are
#'   recorded in the run card.
#' @param source_tables Named list of paths to the full result tables, for looking up features the
#'   payload leaves out.
#' @param run_extra Named list of further fields for the run card.
#' @param trust Declarations for quality checks the code cannot compute, as a named list of
#'   `c(status, note)`, for example
#'   `list(cluster_stability = c("pass", "stable across resolutions 0.4 to 1.0"))`. A declaration
#'   replaces the computed check of the same name. Statuses are `"pass"`, `"warn"`, `"fail"`,
#'   `"not_run"` and `"not_applicable"`.
#' @param technical_panels `TRUE` scores the shipped technical panels (see [llomics_panels()])
#'   against all significant features of each unit. `FALSE` skips this. A named list in the same
#'   shape replaces the shipped panels.
#' @param gene_sets Gene sets for over-representation of the hits: a named list of character
#'   vectors, as returned by [llomics_read_gmt()]. The tool wrappers also accept the path to a GMT
#'   file. `NULL` skips enrichment.
#' @param gene_sets_name Where the gene sets came from, recorded in the run card. The tool wrappers
#'   fill it from the GMT file name.
#' @param p_floor Smallest p-value written. Smaller values, including p-values that underflowed to
#'   zero, are clamped to it and counted in the run card. `NULL` uses the contract default (1e-300).
#'
#' @name llomics-params
#' @keywords internal
NULL
