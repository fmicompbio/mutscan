#' Make pairs plot of selected assay from a SummarizedExperiment object
#'
#' Construct a pairs plot of all columns of a given assay. The lower-triangular
#' panels display the scatter plots, the upper-triangular ones print out
#' the (Pearson or Spearman) correlations, and the diagonal panels show
#' histograms of the respective columns.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @return A \code{ggplot} object.
#'
#' @param se A SummarizedExperiment object, e.g. the output of
#'     \code{summarizeExperiment}
#' @param selAssay Character scalar, the assay to use as the basis for the
#'     pairs plot.
#' @param doLog Logical scalar, whether or not to log-transform the values
#'     before plotting.
#' @param pseudocount Numeric scalar, the pseudocount to add to the values
#'     before log-transforming (if \code{doLog} is \code{TRUE}).
#' @param corMethod Either "pearson" or "spearman", the type of correlation
#'     to calculate.
#' @param histBreaks Numeric scalar, the number of breaks in the histograms
#'     to put in the diagonal panels.
#' @param pointsType Either "points", "smoothscatter", "scattermore" or 
#'     "scattermost" (the latter two require the "scattermore" package to be 
#'     installed), determining the type of plots that will be made.
#' @param corSizeMult,corSizeAdd Numeric scalars determining how the
#'     absolute correlation value is transformed into a font size. The
#'     transformation is corSizeMult * abs(corr) + corSizeAdd.
#' @param pointSize,pointAlpha Numeric scalars determining the size and
#'     opacity of points in the plot.
#' @param colorByCorrelation Logical scalar, indicating whether the correlation
#'     panels should be colored according to the correlation value.
#' @param corrColorRange Numeric vector of length 2, providing the lower and 
#'     upper limits of the color scale when coloring by correlation. Both 
#'     values should be positive; the same range is used for negative 
#'     correlations. If \code{NULL} (the default), the range is inferred from 
#'     the data.
#' @param addIdentityLine Logical scalar, indicating whether the identity line 
#'     should be added (only used if \code{pointsType = "points"}).
#'
#' @importFrom GGally eval_data_col ggpairs
#' @importFrom ggplot2 ggplot annotate theme_void ylim stat_density2d
#' @importFrom ggplot2 scale_fill_continuous geom_point theme_bw theme
#' @importFrom ggplot2 element_blank aes geom_histogram scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous geom_abline after_stat element_rect
#' @importFrom stats cor
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom grDevices hcl.colors rgb colorRamp
#' @importFrom rlang .data
#'
#' @examples
#' se <- readRDS(system.file("extdata", "GSE102901_cis_se.rds", 
#'                           package = "mutscan"))[1:200, ]
#' plotPairs(se)
#' 
plotPairs <- function(se, selAssay = "counts", doLog = TRUE, pseudocount = 1,
                      corMethod = "pearson", histBreaks = 40,
                      pointsType = "points", corSizeMult = 5,
                      corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3,
                      colorByCorrelation = TRUE, corrColorRange = NULL,
                      addIdentityLine = FALSE) {
    
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = selAssay, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = doLog, type = "logical")
    .assertScalar(x = pseudocount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = corMethod, type = "character",
                  validValues = c("pearson", "spearman"))
    .assertScalar(x = histBreaks, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = pointsType, type = "character",
                  validValues = c("smoothscatter", "points", "scattermore",
                                  "scattermost"))
    .assertScalar(x = corSizeMult, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = corSizeAdd, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pointSize, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = pointAlpha, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = colorByCorrelation, type = "logical")
    .assertVector(x = corrColorRange, type = "numeric", rngIncl = c(0, 1), 
                  len = 2, allowNULL = TRUE)
    if (!is.null(corrColorRange)) {
        stopifnot(corrColorRange[2] >= corrColorRange[1])
    }
    .assertScalar(x = addIdentityLine, type = "logical")
    if (pointsType %in% c("scattermore", "scattermost")) {
        .assertPackagesAvailable("scattermore")
    }
    
    ## --------------------------------------------------------------------- ##
    ## Define shared theme elements
    ## --------------------------------------------------------------------- ##
    ggtheme <- list(
        theme_bw(),
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    )
    
    ## --------------------------------------------------------------------- ##
    ## Correlations
    ## --------------------------------------------------------------------- ##
    ## Define function to calculate and display correlations 
    ## (for use with ggpairs)
    cor_fcn <- function(data, mapping, ...) {
        ## Get data
        xData <- eval_data_col(data, mapping$x)
        yData <- eval_data_col(data, mapping$y)
        
        ## Calculate correlation
        mainCor <- cor(xData, yData, method = corMethod, 
                       use = "pairwise.complete.obs")
        transfCor <- (abs(mainCor) - min(corRange)) / 
            (max(corRange) - min(corRange))
        ## Protect against numerical imprecision leading to values outside
        ## the allowed range
        transfCor <- min(1.0, transfCor)
        transfCor <- max(0.0, transfCor)
        
        ## Determine the color
        if (colorByCorrelation) {
            if (mainCor >= 0) {
                cols <- hcl.colors(n = 11, palette = "RdBu")[6:2]
                col <- rgb(colorRamp(
                    cols)(transfCor),
                    maxColorValue = 255)
            } else {
                # nocov start
                cols <- hcl.colors(n = 11, palette = "RdBu")[6:10]
                col <- rgb(colorRamp(
                    cols)(transfCor),
                    maxColorValue = 255)
                # nocov end
            }
        } else {
            col <- "white" #nocov
        }
        
        ## Construct plot
        ggplot(data = data, mapping = mapping) +
            annotate(x = 0.5, y = 0.5, 
                     label = round(mainCor, digits = 3),
                     geom = "text",
                     size = abs(mainCor) * corSizeMult + corSizeAdd) +
            ggtheme + ylim(c(0, 1)) +
            theme(
                panel.background = element_rect(fill = col)
            )
    }
    
    ## --------------------------------------------------------------------- ##
    ## Scatter plots
    ## --------------------------------------------------------------------- ##
    # All the point plotting functions are tested, but for some reason covr 
    # doesn't pick it up - ignore these definitions in the coverage tests
    # nocov start
    if (pointsType == "smoothscatter") {
        ## Define function to create smoothscatter-like plot 
        ## (for use with ggpairs)
        smoothscat <- function(data, mapping, ...) {
            ggplot(data = data, mapping = mapping) +
                stat_density2d(
                    aes(fill = after_stat(
                        .data$density)^0.25), geom = "tile",
                    contour = FALSE, n = 200) +
                scale_fill_continuous(low = "white", 
                                      high = "darkgreen") +
                ggtheme +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0))
        }
    } else if (pointsType == "points") {
        ## Define function to create scatter plot (for use with ggpairs)
        if (addIdentityLine) {
            plotpoints <- function(data, mapping, ...) {
                ggplot(data = data, mapping = mapping) +
                    geom_abline(slope = 1, intercept = 0, 
                                linetype = "dashed", color = "grey") + 
                    geom_point(alpha = pointAlpha, size = pointSize) +
                    ggtheme
            }
        } else {
            plotpoints <- function(data, mapping, ...) {
                ggplot(data = data, mapping = mapping) +
                    geom_point(alpha = pointAlpha, size = pointSize) +
                    ggtheme
            }
        }
    } else if (pointsType == "scattermore") {
        ## Define function to create scattermore plot (for use with ggpairs)
        if (addIdentityLine) {
            plotscattermore <- function(data, mapping, ...) {
                ggplot(data = data, mapping = mapping) +
                    geom_abline(slope = 1, intercept = 0, 
                                linetype = "dashed", color = "grey") + 
                    scattermore::geom_scattermore(alpha = pointAlpha, 
                                                  pointsize = pointSize) +
                    ggtheme
            }
        } else {
            plotscattermore <- function(data, mapping, ...) {
                ggplot(data = data, mapping = mapping) +
                    scattermore::geom_scattermore(alpha = pointAlpha, 
                                                  pointsize = pointSize) +
                    ggtheme
            }
        }
    } else if (pointsType == "scattermost") {
        ## Define function to create scattermost plot (for use with ggpairs)
        if (addIdentityLine) {
            plotscattermost <- function(data, mapping, ...) {
                ## Get data
                xData <- eval_data_col(data, mapping$x)
                yData <- eval_data_col(data, mapping$y)
                
                ggplot(data = data, mapping = mapping) +
                    geom_abline(slope = 1, intercept = 0, 
                                linetype = "dashed", color = "grey") + 
                    scattermore::geom_scattermost(xy = cbind(xData, yData), 
                                                  pointsize = pointSize) +
                    ggtheme
            }
        } else {
            plotscattermost <- function(data, mapping, ...) {
                ## Get data
                xData <- eval_data_col(data, mapping$x)
                yData <- eval_data_col(data, mapping$y)
                
                ggplot(data = data, mapping = mapping) +
                    scattermore::geom_scattermost(xy = cbind(xData, yData), 
                                                  pointsize = pointSize) +
                    ggtheme
            }
        }
    }
    # nocov end
    
    ## Define the function to use for the plots
    if (pointsType == "smoothscatter") {
        lower <- list(continuous = smoothscat)
    } else if (pointsType == "points") {
        lower <- list(continuous = plotpoints)
    } else if (pointsType == "scattermore") {
        lower <- list(continuous = plotscattermore)
    } else if (pointsType == "scattermost") {
        lower <- list(continuous = plotscattermost)
    }
    
    ## --------------------------------------------------------------------- ##
    ## Histogram
    ## --------------------------------------------------------------------- ##
    # nocov start
    diaghist <- function(data, mapping, ...) {
        gg <- ggplot(data = data, mapping = mapping) +
            geom_histogram(fill = "#F5C710", color = "grey50", 
                           bins = histBreaks) +
            ggtheme
        if (pointsType == "smoothscatter") {
            gg <- gg + scale_x_continuous(expand = c(0, 0))
        }
        gg
    }
    # nocov end
    
    ## --------------------------------------------------------------------- ##
    ## Combined plot
    ## --------------------------------------------------------------------- ##
    ## Prepare the data and plot title, depending on whether to log-transform 
    ## or not
    if (doLog) {
        mat <- log10(as.matrix(assay(se, selAssay)) + 
                         pseudocount)
        title <- paste0("log10(", selAssay,
                        ifelse(pseudocount > 0, 
                               paste0(" + ", pseudocount), ""),
                        ")")
    } else {
        mat <- as.matrix(assay(se, selAssay))
        title <- selAssay
    }
    mat <- as.data.frame(mat)
    title <- paste0(title, ", ", toupper(substr(corMethod, 1, 1)),
                    substr(corMethod, 2, nchar(corMethod)), " correlation")
    
    ## Calculate correlations and get range
    allCors <- cor(mat, method = corMethod, 
                   use = "pairwise.complete.obs")
    if (!is.null(corrColorRange)) {
        corRange <- corrColorRange
    } else {
        corRange <- range(abs(allCors[upper.tri(allCors)]))
    }
    if (length(unique(corRange)) == 1) {
        ## Having a range with width 0 leads to problems in the rescaling of 
        ## the correlations
        corRange <- c(0, 1)
    }
    
    ## Construct the pairs plot
    ggpairs(
        data = mat,
        mapping = NULL,
        columns = seq_len(ncol(mat)),
        title = title, xlab = NULL, ylab = NULL,
        upper = list(continuous = cor_fcn),
        lower = lower,
        diag = list(continuous = diaghist),
        progress = FALSE,
        axisLabels = "show"
    )
}

