#' @importFrom ggplot2 ggplot theme_minimal coord_cartesian theme labs 
#'     element_text geom_point 
#' @importFrom rlang .data
.makeScatterPlot <- function(res, xCol, yCol, xLabel, yLabel, 
                             padjCol, padjLabel, padjThreshold, 
                             labelCol = NULL, centerAxis = "x", 
                             pointSize = "small", interactivePlot = FALSE) {
    stopifnot(is(res, "data.frame"))
    stopifnot(length(padjThreshold) == 1 && is.numeric(padjThreshold) && 
                  padjThreshold >= 0)
    stopifnot(length(pointSize) == 1 && is.character(pointSize) && 
                  pointSize %in% c("small", "large"))
    stopifnot(length(interactivePlot) == 1 && is.logical(interactivePlot))
    if (interactivePlot && !requireNamespace("plotly", quietly = TRUE)) {
        stop("The 'plotly' package is required to make interactive plots.")
    }
    xr <- range(res[[xCol]], na.rm = TRUE)
    yr <- range(res[[yCol]], na.rm = TRUE)
    if ("x" %in% centerAxis) {
        xr <- c(-max(abs(xr), na.rm = TRUE), max(abs(xr), na.rm = TRUE))
    } 
    if ("y" %in% centerAxis) {
        yr <- c(-max(abs(yr), na.rm = TRUE), max(abs(yr), na.rm = TRUE))
    }
    gg <- ggplot2::ggplot(res, aes(x = .data[[xCol]], y = .data[[yCol]], 
                                   label = .data[[labelCol]])) + 
        ggplot2::theme_minimal() + ggplot2::coord_cartesian(xlim = xr, ylim = yr) + 
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                       axis.title = ggplot2::element_text(size = 14),
                       title = ggplot2::element_text(size = 14)) + 
        ggplot2::labs(x = xLabel, y = yLabel)
    if (pointSize == "large") {
        gg <- gg + 
            ggplot2::geom_point(fill = "lightgrey", color = "grey", pch = 21, size = 1.5)
    } else {
        gg <- gg + 
            ggplot2::geom_point(color = "black", size = 0.25, alpha = 0.5)
    }
    if (any(res[[padjCol]] < padjThreshold)) {
        gg <- gg + 
            ggplot2::labs(subtitle = paste0("Points are indicated in red if ", 
                                            padjCol, "<", padjThreshold))
        if (pointSize == "large") {
            gg <- gg + 
                ggplot2::geom_point(
                    data = res[res[[padjCol]] < padjThreshold, , drop = FALSE],
                    fill = "red", color = "grey", pch = 21, size = 1.5
                )
        } else {
            gg <- gg + 
                ggplot2::geom_point(
                    data = res[res[[padjCol]] < padjThreshold, , drop = FALSE],
                    color = "red", size = 0.25
                )
                
        }
    }
    if (interactivePlot) {
        plotly::ggplotly(gg)
    } else {
        gg
    }
}

.getColName <- function(res, validValues, aspect) {
    stopifnot(is(res, "data.frame"))
    colName = grep(paste(paste0("^", validValues, "$"), collapse = "|"), 
                   colnames(res), value = TRUE)
    if (length(colName) == 0) {
        stop("No suitable column found for ", aspect, ", one of the ",
             "following expected: ", paste(validValues, collapse = ", "))
    }
    if (length(colName) > 1) {
        warning("Multiple suitable columns found for ", aspect, 
                "; choosing the first: ", colName[1])
        colName <- colName[1]
    }
    colName
}

#' Construct an MA (mean-difference) plot
#' 
#' @param res \code{data.frame} (typically output from 
#'     \code{calculateRelativeFC()}) with columns corresponding to the 
#'     average abundance (\code{logCPM} or \code{AveExpr}), log-fold 
#'     change (\code{logFC}) and significance (\code{FDR} or 
#'     \code{adj.P.Val}). 
#' @param padjThreshold Numeric scalar indicating the adjusted p-value 
#'     threshold to use for coloring the points. All features with 
#'     adjusted p-value below the treshold will be shown in red.
#' @param pointSize Either \code{"small"} or \code{"large"}, indicating 
#'     which of the two available plot styles that will be used.
#' @param interactivePlot Logical scalar, indicating whether an 
#'     interactive plot should be returned, in which one can hover 
#'     over the individual points and obtain further information. 
#' 
#' @return If \code{interactivePlot} is \code{TRUE}, a \code{plotly} 
#'     object. If \code{interactivePlot} is \code{FALSE}, a \code{ggplot2}
#'     object. 
#' 
#' @author Charlotte Soneson
#' @export
#' 
#' @examples 
#' se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds", 
#'                           package = "mutscan"))[1:200, ]
#' design <- model.matrix(~ Replicate + Condition, 
#'                        data = SummarizedExperiment::colData(se))
#' res <- calculateRelativeFC(se, design, coef = "Conditioncis_output")
#' plotMeanDiff(res, pointSize = "large")
#' 
plotMeanDiff <- function(res, padjThreshold = 0.05, pointSize = "small",
                         interactivePlot = FALSE) {
    ## Get the columns to use
    xCol <- .getColName(res, validValues = c("logCPM", "AveExpr"), 
                        aspect = "x-axis")
    yCol <- .getColName(res, validValues = c("logFC"), 
                        aspect = "y-axis")
    padjCol <- .getColName(res, validValues = c("FDR", "adj.P.Val"),
                           aspect = "highlighting")
    res$feature <- rownames(res)
    
    .makeScatterPlot(res = res, xCol = xCol, yCol = yCol, 
                     xLabel = xCol, yLabel = yCol,
                     padjCol = padjCol, padjLabel = padjCol, 
                     padjThreshold = padjThreshold, 
                     labelCol = "feature", centerAxis = "y", 
                     pointSize = pointSize,
                     interactivePlot = interactivePlot)
}

#' Construct a volcano plot
#' 
#' @param res \code{data.frame} (typically output from 
#'     \code{calculateRelativeFC()}) with columns corresponding to the 
#'     log-fold change (\code{logFC}), p-value (\code{PValue} or 
#'     \code{P.Value}) and significance (\code{FDR} or \code{adj.P.Val}). 
#' @param padjThreshold Numeric scalar indicating the adjusted p-value 
#'     threshold to use for coloring the points. All features with 
#'     adjusted p-value below the treshold will be shown in red.
#' @param pointSize Either \code{"small"} or \code{"large"}, indicating 
#'     which of the two available plot styles that will be used.
#' @param interactivePlot Logical scalar, indicating whether an 
#'     interactive plot should be returned, in which one can hover 
#'     over the individual points and obtain further information. 
#' 
#' @return If \code{interactivePlot} is \code{TRUE}, a \code{plotly} 
#'     object. If \code{interactivePlot} is \code{FALSE}, a \code{ggplot2}
#'     object. 
#' 
#' @author Charlotte Soneson
#' @export
#' 
#' @examples 
#' se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds", 
#'                           package = "mutscan"))[1:200, ]
#' design <- model.matrix(~ Replicate + Condition, 
#'                        data = SummarizedExperiment::colData(se))
#' res <- calculateRelativeFC(se, design, coef = "Conditioncis_output")
#' plotVolcano(res, pointSize = "large")
#' 
plotVolcano <- function(res, padjThreshold = 0.05, pointSize = "small",
                        interactivePlot = FALSE) {
    ## Get the columns to use
    xCol <- .getColName(res, validValues = c("logFC"), 
                        aspect = "x-axis")
    yCol <- .getColName(res, validValues = c("PValue", "P.Value"), 
                        aspect = "y-axis")
    padjCol <- .getColName(res, validValues = c("FDR", "adj.P.Val"),
                           aspect = "highlighting")
    res$negLog10P <- -log10(res[[yCol]])
    res$feature <- rownames(res)
    
    .makeScatterPlot(res = res, xCol = xCol, yCol = "negLog10P", 
                     xLabel = xCol, yLabel = paste0("-log10(", yCol, ")"),
                     padjCol = padjCol, padjLabel = padjCol, 
                     padjThreshold = padjThreshold, 
                     labelCol = "feature", centerAxis = "x", 
                     pointSize = pointSize,
                     interactivePlot = interactivePlot)
}
