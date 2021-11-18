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
#' @author Charlotte 
#' 
#' @returns Invisibly, the path to the generated html file. 
#' 
#' @importFrom rmarkdown render
generateQCReport <- function(se, outFile, forceOverwrite = FALSE) {
    
    if (tools::file_ext(outFile) != "html") {
        stop("'outFile' must have the file extension '.html'. ")
    }
    outDir <- dirname(outFile)
    if (!dir.exists(outDir)) {
        dir.create(outDir, recursive = TRUE)
    }
    if (file.exists(outFile) && !forceOverwrite) {
        stop(outFile, " already exists and forceOverwrite = FALSE, stopping.")
    } else if (file.exists(outFile) && forceOverwrite) {
        message(outFile, " already exists but forceOverwrite = TRUE, overwriting.")
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