# This script is provided as a utility via the swissknife package 
# (https://github.com/fmicompbio/swissknife). This script is provided under 
# the MIT license, and package authors are permitted to 
# include the code as-is in other packages, as long as this note and the 
# information provided below crediting the authors of the respective 
# functions is retained. 

#' Utility function to check validity of scalar variable values.
#' 
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria. 
#'  
#' @param x The variable to be checked.
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x}, 
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x}, 
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' 
#' @author Michael Stadler
#' @noRd
#' @keywords internal
.assertScalar <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL) {
    args <- lapply(sys.call()[-1], as.character)
    xname <- if ("x" %in% names(args)) args$x else "argument"
    
    ## Check arguments
    stopifnot(is.null(type) || (length(type) == 1 && is.character(type)))
    stopifnot(is.null(rngIncl) || (length(rngIncl) == 2 && is.numeric(rngIncl)))
    stopifnot(is.null(rngExcl) || (length(rngExcl) == 2 && is.numeric(rngExcl)))
    
    if (length(x) != 1L) {
        stop("'", xname, "' must be a scalar value (length one)")
    }
    
    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }
    
    if (!is.null(type) && !is(x, type)) {
        stop("'", xname, "' must be of type '", type, "'")
    }
    
    if (!is.null(rngIncl) && is.numeric(rngIncl) && length(rngIncl) == 2L &&
        (x < rngIncl[1] || x > rngIncl[2])) {
        stop("'", xname, "' must be within [", rngIncl[1], ",", 
             rngIncl[2], "] (inclusive)")
    }
    
    if (!is.null(rngExcl) && is.numeric(rngExcl) && length(rngExcl) == 2L &&
        (x <= rngExcl[1] || x >= rngExcl[2])) {
        stop("'", xname, "' must be within (", rngExcl[1], ",", 
             rngExcl[2], ") (exclusive)")
    }
    
    if (!is.null(validValues) && !(x %in% validValues)) {
        stop("'", xname, "' must be one of: ", paste(validValues, 
                                                     collapse = ", "))
    }
    
    return(invisible(TRUE))
}

#' Utility function to check validity of vector variable values.
#' 
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria. 
#'  
#' @param x The variable to be checked
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x}, 
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x}, 
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' @param len The required length of \code{x}.
#' @param rngLen The allowed range for the length of \code{x}.
#' 
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
.assertVector <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          len = NULL, 
                          rngLen = NULL) {
    args <- lapply(sys.call()[-1], as.character)
    xname <- if ("x" %in% names(args)) args$x else "argument"
    
    ## Check arguments
    stopifnot(is.null(type) || (length(type) == 1 && is.character(type)))
    stopifnot(is.null(rngIncl) || (length(rngIncl) == 2 && is.numeric(rngIncl)))
    stopifnot(is.null(rngExcl) || (length(rngExcl) == 2 && is.numeric(rngExcl)))
    stopifnot(is.null(len) || (length(len) == 1 && is.numeric(len)))
    stopifnot(is.null(rngLen) || (length(rngLen) == 2 && is.numeric(rngLen)))
    
    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }
    
    if (!is.null(type) && !is(x, type)) {
        stop("'", xname, "' must be of class '", type, "'")
    }
    
    if (!is.null(rngIncl) && is.numeric(rngIncl) && length(rngIncl) == 2L &&
        any(x < rngIncl[1] | x > rngIncl[2])) {
        stop("values in '", xname, "' must be within [", rngIncl[1], ",", 
             rngIncl[2], "] (inclusive)")
    }
    
    if (!is.null(rngExcl) && is.numeric(rngExcl) && length(rngExcl) == 2L &&
        any(x <= rngExcl[1] | x >= rngExcl[2])) {
        stop("values in '", xname, "' must be within (", rngExcl[1], ",", 
             rngExcl[2], ") (exclusive)")
    }
    
    if (!is.null(len) && is.numeric(len) && length(len) == 1L && 
        length(x) != len) {
        stop("'", xname, "' must have length ", len)
    }
    
    if (!is.null(rngLen) && is.numeric(rngLen) && length(rngLen) == 2L &&
        any(length(x) < rngLen[1] | length(x) > rngLen[2])) {
        stop("length of '", xname, "' must be within [", rngLen[1], ",", 
             rngLen[2], "] (inclusive)")
    }
    
    if (!is.null(validValues) && !all(x %in% validValues)) {
        if (length(validValues) > 15) {
            vvPrint <- paste(c(validValues[seq_len(15)], 
                               "...(truncated)"),
                             collapse = ", ")
        } else {
            vvPrint <- paste(validValues, collapse = ", ")
        }
        stop("All values in '", xname, "' must be one of: ", vvPrint)
    }
    
    return(invisible(TRUE))
}
