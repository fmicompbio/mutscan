#' Generate QC report
#'
#' @param se A \code{SummarizedExperiment} object, typically generated with
#'     \code{summarizeExperiment()}.
#' @param outFile Character string providing the name of the output file.
#'     Should have the extension \code{.html}.
#' @param reportTitle Character string specifying the title of the QC report.
#' @param forceOverwrite Logical scalar, indicating whether an existing file
#'     with the same name as \code{outFile} should be overwritten.
#' @param ... Additional parameters to be forwarded to
#'     \code{\link[rmarkdown]{render}}, for example \code{quiet = TRUE}.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns Invisibly, the path to the generated html file.
#'
#' @seealso \code{\link[rmarkdown]{render}} used to render the html output file.
#'
#' @importFrom rmarkdown render
#' @importFrom xfun Rscript_call
#' @importFrom DT datatable
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr bind_rows full_join
#' @importFrom S4Vectors metadata
#'
#' @examples
#' ## Load SummarizedExperiment object
#' se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds",
#'                           package = "mutscan"))
#' ## Define output file
#' outfile <- tempfile(fileext = ".html")
#'
#' ## Generate QC report
#' generateQCReport(se, outfile)
#'
generateQCReport <- function(se, outFile, reportTitle = "mutscan QC report",
                             forceOverwrite = FALSE, ...) {

    ## --------------------------------------------------------------------- ##
    ## Check that input arguments are valid
    ## --------------------------------------------------------------------- ##
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = outFile, type = "character")
    .assertScalar(x = forceOverwrite, type = "logical")
    if (tools::file_ext(outFile) != "html") {
        stop("'outFile' must have the file extension '.html'.")
    }
    outDir <- dirname(outFile)
    if (!dir.exists(outDir)) {
        dir.create(outDir, recursive = TRUE)
    }
    if (file.exists(outFile)) {
        if (!forceOverwrite) {
            stop(outFile, " already exists and forceOverwrite = FALSE, stopping.")
        } else {
            message(outFile, " already exists but forceOverwrite = TRUE, overwriting.")
        }
    }

    ## --------------------------------------------------------------------- ##
    ## Render the Rmd file
    ## --------------------------------------------------------------------- ##
    args <- list()
    args$input <- system.file("templates", "qc_report_template.Rmd",
                              package = "mutscan")
    args$output_format <- "html_document"
    args$output_file <- basename(outFile)
    args$output_dir <- outDir
    args$intermediates_dir <- outDir
    args$run_pandoc <- TRUE
    args$params <- list(se = se, title = reportTitle)
    if (any(names(args) %in% ...names())) {
        warning("Will ignore arguments provided in '...' that are set ",
                "internally by generateQCReport: ",
                paste(intersect(names(args), ...names()), collapse = ", "))
    }
    args <- c(args, list(...)[setdiff(...names(), names(args))])

    outputReport <- xfun::Rscript_call(
        rmarkdown::render,
        args
    )

    ## --------------------------------------------------------------------- ##
    ## Return (invisibly) the path to the rendered html file
    ## --------------------------------------------------------------------- ##
    invisible(outFile)
}