#' Construct an (interactive) scatterplot
#' 
#' This function provides a convenience wrapper for constructing a scatter 
#' plot of two columns in a data.frame, and color the points by whether 
#' the values in a third column are below a provided value. 
#' Users will interact with the function via the higher-level wrappers 
#' \code{plotVolcano()} and \code{plotMeanDiff()}, which both call this 
#' function internally. 
#' 
#' @param res A \code{data.frame} with values to plot.
#' @param xCol,yCol Character scalars defining the columns to plot on the
#'     x- and y-axis, respectively.
#' @param xLabel,yLabel Character scalars providing the x- and y-axis 
#'     titles, respectively.
#' @param colorCol Character scalar defining the column to use to color 
#'     the points (typically the adjusted p-value). 
#' @param colorThreshold Numeric scalar. If the value in \code{colorCol} 
#'     is lower than colorThreshold, the point will be indicated in red.
#' @param labelCol Character scalar defining the column to use to label 
#'     individual points in the plot.
#' @param labelValues Character vector defining the points to highlight with 
#'     labels. Must correspond to values in the \code{labelCol} column.
#' @param centerAxis Character vector, a subset of c("x", "y") or \code{NULL}.
#'     Indicates the axis that should be centered (i.e., where the axis 
#'     limits should be extended to put 0 in the center).
#' @param pointSize Character scalar, either \code{"small"} or \code{"large"},
#'     defining the style to use to plot the points. 
#' @param interactivePlot Logical scalar, whether to return an interactive 
#'     plot. 
#' 
#' @return If \code{interactivePlot} is \code{TRUE}, a \code{plotly} 
#'     object. If \code{interactivePlot} is \code{FALSE}, a \code{ggplot2}
#'     object. 
#' 
#' @author Charlotte Soneson
#' @keywords internal
#' @noRd
#' 
#' @importFrom ggplot2 ggplot theme_minimal coord_cartesian theme labs 
#'     element_text geom_point aes
#' @importFrom rlang .data
#' @importFrom ggrepel geom_text_repel
.plotScatter <- function(res, xCol, yCol, xLabel = xCol, yLabel = yCol, 
                         colorCol, colorThreshold, 
                         labelCol = NULL, labelValues = c(), centerAxis = "x", 
                         pointSize = "small", interactivePlot = FALSE) {
    .assertVector(x = res, type = "data.frame")
    .assertScalar(x = xCol, type = "character", validValues = colnames(res))
    .assertScalar(x = yCol, type = "character", validValues = colnames(res))
    .assertScalar(x = xLabel, type = "character")
    .assertScalar(x = yLabel, type = "character")
    .assertScalar(x = colorCol, type = "character", 
                  validValues = colnames(res))
    .assertScalar(x = colorThreshold, type = "numeric")
    .assertScalar(x = labelCol, type = "character", 
                  validValues = colnames(res), allowNULL = TRUE)
    .assertVector(x = labelValues, type = "character", allowNULL = TRUE)
    .assertVector(x = centerAxis, type = "character", 
                  validValues = c("x", "y"), allowNULL = TRUE)
    .assertScalar(x = pointSize, type = "character", 
                  validValues = c("small", "large"))
    .assertScalar(x = interactivePlot, type = "logical")
    
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
    gg <- ggplot2::ggplot(res, ggplot2::aes(x = .data[[xCol]], y = .data[[yCol]])) + 
        ggplot2::theme_minimal() + ggplot2::coord_cartesian(xlim = xr, ylim = yr) + 
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                       axis.title = ggplot2::element_text(size = 14),
                       title = ggplot2::element_text(size = 14)) + 
        ggplot2::labs(x = xLabel, y = yLabel)
    if (!is.null(labelCol)) {
        gg <- gg + 
            ggplot2::aes(label = .data[[labelCol]])
    }
    if (pointSize == "large") {
        gg <- gg + 
            ggplot2::geom_point(fill = "lightgrey", color = "grey", pch = 21, size = 1.5)
    } else {
        gg <- gg + 
            ggplot2::geom_point(color = "black", size = 0.25, alpha = 0.5)
    }
    if (any(res[[colorCol]] < colorThreshold)) {
        gg <- gg + 
            ggplot2::labs(subtitle = paste0("Points are indicated in red if ", 
                                            colorCol, "<", colorThreshold))
        if (pointSize == "large") {
            gg <- gg + 
                ggplot2::geom_point(
                    data = res[res[[colorCol]] < colorThreshold, , drop = FALSE],
                    fill = "red", color = "grey", pch = 21, size = 1.5
                )
        } else {
            gg <- gg + 
                ggplot2::geom_point(
                    data = res[res[[colorCol]] < colorThreshold, , drop = FALSE],
                    color = "red", size = 0.25
                )
        }
    }
    if (interactivePlot) {
        plotly::ggplotly(gg)
    } else {
        if (!is.null(labelCol) && length(labelValues) > 0) {
            .assertVector(x = labelValues, type = "character", validValues = res[[labelCol]])
            gg <- gg + 
                ggrepel::geom_text_repel(
                    data = res[res[[labelCol]] %in% labelValues, , drop = FALSE],
                    max.overlaps = Inf, size = 4, min.segment.length = 0.1)
        }
        gg
    }
}

.getColName <- function(res, validValues, aspect) {
    .assertVector(x = res, type = "data.frame")
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
#' @param meanCol,logFCCol,padjCol Character scalars indicating which 
#'     columns (from \code{res}) that will be used to represent the 
#'     mean value (x-axis), logFC (y-axis) and adjusted p-value (used 
#'     for coloring). If \code{NULL} (default), pre-specified values 
#'     will be used depending on the available columns 
#'     (\code{"logCPM"} or \code{"AveExpr"}, \code{"logFC"}, and 
#'     \code{"FDR"} or \code{"adj.P.Val"}, respectively).
#' @param padjThreshold Numeric scalar indicating the adjusted p-value 
#'     threshold to use for coloring the points. All features with 
#'     adjusted p-value below the treshold will be shown in red.
#' @param pointSize Either \code{"small"} or \code{"large"}, indicating 
#'     which of the two available plot styles that will be used.
#' @param interactivePlot Logical scalar, indicating whether an 
#'     interactive plot should be returned, in which one can hover 
#'     over the individual points and obtain further information. 
#' @param nTopToLabel Numeric scalar, indicating the number of points that 
#'     should be labeled in the plot. The points will be ranked by the 
#'     \code{padjCol} column, and the top \code{nTopToLabel} values will 
#'     be labeled by the corresponding row names. Only used if 
#'     \code{interactivePlot} is \code{FALSE}.
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
#' plotMeanDiff(res, pointSize = "large", nTopToLabel = 3)
#' 
plotMeanDiff <- function(res, meanCol = NULL, logFCCol = NULL, 
                         padjCol = NULL, padjThreshold = 0.05, pointSize = "small",
                         interactivePlot = FALSE, nTopToLabel = 0) {
    .assertScalar(nTopToLabel, type = "numeric", rngIncl = c(0, Inf))
    ## Get the columns to use
    if (is.null(meanCol)) {
        meanCol <- .getColName(res, validValues = c("logCPM", "AveExpr"), 
                               aspect = "x-axis")
    }
    if (is.null(logFCCol)) {
        logFCCol <- .getColName(res, validValues = c("logFC"), 
                                aspect = "y-axis")
    }
    if (is.null(padjCol)) {
        padjCol <- .getColName(res, validValues = c("FDR", "adj.P.Val"),
                                aspect = "highlighting")
    }
    res$feature <- rownames(res)
    if (nTopToLabel > 0) {
        labelValues <- res[order(res[[padjCol]]), "feature"][seq_len(nTopToLabel)]
    } else {
        labelValues <- c()
    }
    
    .plotScatter(res = res, xCol = meanCol, yCol = logFCCol, 
                 xLabel = meanCol, yLabel = logFCCol,
                 colorCol = padjCol, 
                 colorThreshold = padjThreshold, 
                 labelCol = "feature", labelValues = labelValues,
                 centerAxis = "y", pointSize = pointSize,
                 interactivePlot = interactivePlot)
}

#' Construct a volcano plot
#' 
#' @param res \code{data.frame} (typically output from 
#'     \code{calculateRelativeFC()}) with columns corresponding to the 
#'     log-fold change (\code{logFC}), p-value (\code{PValue} or 
#'     \code{P.Value}) and significance (\code{FDR} or \code{adj.P.Val}). 
#' @param logFCCol,pvalCol,padjCol Character scalars indicating which 
#'     columns (from \code{res}) that will be used to represent the 
#'     logFC (x-axis), p-value (y-axis) and adjusted p-value (used 
#'     for coloring). If \code{NULL} (default), pre-specified values 
#'     will be used depending on the available columns 
#'     (\code{"logFC"}, \code{"PValue"} or \code{"P.Value"}, and 
#'     \code{"FDR"} or \code{"adj.P.Val"}, respectively).
#' @param padjThreshold Numeric scalar indicating the adjusted p-value 
#'     threshold to use for coloring the points. All features with 
#'     adjusted p-value below the treshold will be shown in red.
#' @param pointSize Either \code{"small"} or \code{"large"}, indicating 
#'     which of the two available plot styles that will be used.
#' @param interactivePlot Logical scalar, indicating whether an 
#'     interactive plot should be returned, in which one can hover 
#'     over the individual points and obtain further information. 
#' @param nTopToLabel Numeric scalar, indicating the number of points that 
#'     should be labeled in the plot. The points will be ranked by the 
#'     \code{padjCol} column, and the top \code{nTopToLabel} values will 
#'     be labeled by the corresponding row names. Only used if 
#'     \code{interactivePlot} is \code{FALSE}.
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
#' plotVolcano(res, pointSize = "large", nTopToLabel = 3)
#' 
plotVolcano <- function(res, logFCCol = NULL, pvalCol = NULL, padjCol = NULL,
                        padjThreshold = 0.05, pointSize = "small",
                        interactivePlot = FALSE, nTopToLabel = 0) {
    .assertScalar(nTopToLabel, type = "numeric", rngIncl = c(0, Inf))
    ## Get the columns to use
    if (is.null(logFCCol)) {
        logFCCol <- .getColName(res, validValues = c("logFC"), 
                                aspect = "x-axis")
    }
    if (is.null(pvalCol)) {
        pvalCol <- .getColName(res, validValues = c("PValue", "P.Value"), 
                            aspect = "y-axis")
    }
    if (is.null(padjCol)) {
        padjCol <- .getColName(res, validValues = c("FDR", "adj.P.Val"),
                               aspect = "highlighting")
    }
    .assertScalar(pvalCol, type = "character", validValues = colnames(res))
    res$negLog10P <- -log10(res[[pvalCol]])
    res$feature <- rownames(res)
    if (nTopToLabel > 0) {
        labelValues <- res[order(res[[padjCol]]), "feature"][seq_len(nTopToLabel)]
    } else {
        labelValues <- c()
    }
    
    .plotScatter(res = res, xCol = logFCCol, yCol = "negLog10P", 
                 xLabel = logFCCol, yLabel = paste0("-log10(", pvalCol, ")"),
                 colorCol = padjCol, 
                 colorThreshold = padjThreshold, 
                 labelCol = "feature", labelValues = labelValues,
                 centerAxis = "x", pointSize = pointSize,
                 interactivePlot = interactivePlot)
}
