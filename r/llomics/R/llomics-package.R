#' llomics: turn omics pipeline results into payloads for LLMs
#'
#' Builds small, self-describing JSON payloads from the results of statistical steps in omics
#' pipelines. A payload holds a run card (what was run and with which filters), quality checks,
#' and one compact state per unit of judgment: a cluster with its markers, or a contrast with its
#' up- and down-regulated features. A language model can be asked about each unit, and code
#' decides what its answers may settle.
#'
#' Wrappers that run the test and build the payload in one call:
#' * [llomics_FindAllMarkers()] and [llomics_FindMarkers()] for Seurat,
#' * [llomics_limma()] for a limma fit.
#'
#' For any other tool, rename its result columns and call [llomics_dge_payload()] or
#' [llomics_clusters_payload()]. Write a payload with [llomics_write()] and check it with
#' [llomics_validate()].
#'
#' The rules a payload follows are in [llomics_contract()]; the technical gene panels scored
#' against every hit list are in [llomics_panels()].
#'
#' @keywords internal
"_PACKAGE"
