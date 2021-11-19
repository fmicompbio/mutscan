#' Generate QC report
#' 
#' @param se A \code{SummarizedExperiment} object, typically generated with 
#'     \code{summarizeExperiment()}.
#' @param outFile Character string providing the name of the output file. 
#'     Should have the extension \code{.html}.
#' @param forceOverwrite Logical scalar, indicating whether an existing file
#'     with the same name as \code{outFile} should be overwritten. 
#' 
#' @export
#' @author Charlotte Soneson
#' 
#' @returns Invisibly, the path to the generated html file. 
#' 
#' @importFrom rmarkdown render
#' @importFrom xfun Rscript_call
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
generateQCReport <- function(se, outFile, forceOverwrite = FALSE) {
    
    ## --------------------------------------------------------------------- ##
    ## Check that input arguments are valid
    ## --------------------------------------------------------------------- ##
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot(length(outFile) == 1 && is.character(outFile))
    stopifnot(length(forceOverwrite) == 1 && is.logical(forceOverwrite))
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
    args$quiet <- FALSE
    args$run_pandoc <- TRUE
    args$params <- list(se = se)
    
    outputReport <- xfun::Rscript_call(
        rmarkdown::render,
        args
    )
    
    ## --------------------------------------------------------------------- ##
    ## Return (invisibly) the path to the rendered html file
    ## --------------------------------------------------------------------- ##
    invisible(outFile)
}