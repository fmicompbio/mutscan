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
#' @param se A SummarizedExperiment object, e.g. the output of 
#'   \code{summarizeExperiment}
#' @param selAssay Character scalar, the assay to use as the basis for the 
#'   pairs plot.
#' @param doLog Logical scalar, whether or not to log-transform the values 
#'   before plotting.
#' @param pseudocount Numeric scalar, the pseudocount to add to the values 
#'   before log-transforming (if \code{doLog} is \code{TRUE}).
#' @param corMethod Either "pearson" or "spearman", the type of correlation 
#'   to calculate.
#' @param histBreaks Numeric scalar, the number of breaks in the histograms
#'   to put in the diagonal panels.
#' @param pointsType Either "points" or "smoothscatter", determining the 
#'   type of plots that will be made.
#' @param corSizeMult,corSizeAdd Numeric scalars determining how the 
#'   absolute correlation value is transformed into a font size. The 
#'   transformation is corSizeMult * abs(corr) + corSizeAdd.
#' @param pointSize,pointAlpha Numeric scalars determining the size and 
#'   opacity of points in the plot.
#'   
#' @importFrom GGally eval_data_col ggpairs wrap
#' @importFrom ggplot2 ggplot annotate theme_void ylim stat_density2d 
#'   scale_fill_continuous geom_point theme_bw theme element_blank aes
#' @importFrom stats cor
#' @importFrom SummarizedExperiment assayNames assay
#' 
plotPairs <- function(se, selAssay = "counts", doLog = TRUE, pseudocount = 1,
                      corMethod = "pearson", histBreaks = 40,
                      pointsType = "points", corSizeMult = 5, 
                      corSizeAdd = 2, pointSize = 0.1, pointAlpha = 0.3) {
  
  stopifnot(methods::is(se, "SummarizedExperiment"))
  stopifnot(is.character(selAssay) && length(selAssay) == 1 && 
              selAssay %in% SummarizedExperiment::assayNames(se))
  stopifnot(is.logical(doLog) && length(doLog) == 1)
  stopifnot(is.numeric(pseudocount) && length(pseudocount) == 1 && 
              pseudocount >= 0)
  stopifnot(is.character(corMethod) && length(corMethod) == 1 && 
              corMethod %in% c("pearson", "spearman"))
  stopifnot(is.numeric(histBreaks) && length(histBreaks) == 1 && 
              histBreaks > 0)
  stopifnot(is.character(pointsType) && length(pointsType) == 1 && 
              pointsType %in% c("smoothscatter", "points"))
  stopifnot(is.numeric(corSizeMult) && length(corSizeMult) == 1 && 
              corSizeMult > 0)
  stopifnot(is.numeric(corSizeAdd) && length(corSizeAdd) == 1 && 
              corSizeAdd >= 0)
  stopifnot(is.numeric(pointSize) && length(pointSize) == 1 && 
              pointSize > 0)
  stopifnot(is.numeric(pointAlpha) && length(pointAlpha) == 1 && 
              pointAlpha > 0)
  
  ## Define function to calculate and display correlations (for use with ggpairs)
  cor_fcn <- function(data, mapping, ...) {
    ## Get data
    xData <- GGally::eval_data_col(data, mapping$x)
    yData <- GGally::eval_data_col(data, mapping$y)

    ## Calculate correlation
    mainCor = stats::cor(xData, yData, method = corMethod)
    
    ## Construct plot
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::annotate(x = 0.5, y = 0.5, label = round(mainCor, digits = 3), 
                        geom = "text", 
                        size = abs(mainCor) * corSizeMult + corSizeAdd) +
      ggplot2::theme_void() + ggplot2::ylim(c(0, 1))
  }
  
  ## Define function to create smoothscatter-like plot (for use with ggpairs)
  smoothscat <- function(data, mapping, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::stat_density2d(ggplot2::aes(fill = ..density..^0.25), geom = "tile", 
                              contour = FALSE, n = 200) +
      ggplot2::scale_fill_continuous(low = "white", high = "dodgerblue4") + 
      ggplot2::geom_point(alpha = 0.1, shape = 20, size = pointSize, color = "grey50")
  }
  
  ## Define the function to use for the plots
  if (pointsType == "smoothscatter") {
    lower <- list(continuous = smoothscat)
  } else if (pointsType == "points") {
    lower <- list(continuous = GGally::wrap(pointsType, alpha = pointAlpha, 
                                            size = pointSize))
  } else {
    stop("Invalid 'pointsType'")
  }
  
  ## Prepare the data and plot title, depending on whether to log-transform or not
  if (doLog) {
    mat <- log10(as.matrix(SummarizedExperiment::assay(se, selAssay)) + pseudocount)
    title <- paste0("log10(", selAssay, 
                    ifelse(pseudocount > 0, paste0(" + ", pseudocount), ""),
                    ")")
  } else {
    mat <- as.matrix(assay(se, selAssay))
    title <- selAssay
  }
  mat <- as.data.frame(mat)
  title <- paste0(title, ", ", toupper(substr(corMethod, 1, 1)), 
                  substr(corMethod, 2, nchar(corMethod)), " correlation")
  
  ## Construct the pairs plot
  GGally::ggpairs(
    data = mat,
    mapping = NULL,
    columns = seq_len(ncol(mat)),
    title = title, xlab = NULL, ylab = NULL,
    upper = list(continuous = cor_fcn),
    lower = lower,
    diag = list(continuous = GGally::wrap("barDiag", fill = "cyan", 
                                          color = "grey50", bins = histBreaks)),
    progress = FALSE,
    axisLabels = "show"
  ) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

