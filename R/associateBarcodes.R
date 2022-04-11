#' Associate barcode(s) to variants
#' 
#' @param countTable A data frame containing barcode sequences, variant 
#'     sequences and corresponding read counts for each barcode-variant pair. 
#' @param barcodeCol,variantCol,countCol Character scalars indicating the 
#'     names of the columns of \code{countTable} containing barcode 
#'     sequences, variant IDs and read counts, respectively. 
#' @param minCount Numeric scalar. Pairs with fewer reads will be flagged 
#'     for exclusion.
#' @param minFracOfBarcodeCount Numeric scalar. Pairs where the read count 
#'     is lower than \code{minFracOfBarcodeCount} * total barcode count will 
#'     be flagged for exclusion. 
#' 
#' @export
#' @author Charlotte Soneson
#' 
#' @return A list with two elements: \code{counts} (a data frame with counts)
#' and \code{plots} (a list of ggplot objects).
#' 
#' @importFrom ggplot2 scale_color_manual
#' @importFrom dplyr %>% group_by mutate ungroup filter arrange desc
#'     summarize
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_histogram theme_bw labs
#'     scale_color_manual scale_x_log10 scale_y_log10 geom_point 
#'     
associateBarcodes <- function(countTable, barcodeCol, variantCol, 
                              countCol = "nbrReads", minCount = 10, 
                              minFracOfBarcodeCount = 0.5) {
    .assertScalar(x = barcodeCol, type = "character", 
                  validValues = colnames(countTable))
    .assertVector(x = countTable[[barcodeCol]], type = "character")
    .assertScalar(x = variantCol, type = "character", 
                  validValues = colnames(countTable))
    .assertVector(x = countTable[[variantCol]], type = "character")
    .assertScalar(x = countCol, type = "character", 
                  validValues = colnames(countTable))
    .assertVector(x = countTable[[countCol]], type = "numeric")
    .assertScalar(x = minCount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minFracOfBarcodeCount, type = "numeric", rngIncl = c(0, 1))
    
    ## Initialize plot list
    plots <- list()
    
    ## Add columns with total counts for barcode/variant, with and without 
    ## filtering
    countTable <- countTable %>%
        ## first aggregate all rows with the same barcode/variant pair, 
        ## in case that wasn't already done
        dplyr::group_by(.data[[barcodeCol]], .data[[variantCol]]) %>%
        dplyr::summarize("{countCol}" := sum(.data[[countCol]])) %>%
        dplyr::ungroup() %>%
        ## flag pairs with < minCount reads
        dplyr::mutate(keep = (.data[[countCol]] >= minCount)) %>%
        ## add total barcode counts
        dplyr::group_by(.data[[barcodeCol]]) %>%
        dplyr::mutate(totalBarcodeCountAll = sum(.data[[countCol]])) %>%
        dplyr::mutate(totalBarcodeCountFilt = sum(.data[[countCol]][.data$keep])) %>%
        dplyr::ungroup() %>%
        ## add total variant counts
        dplyr::group_by(.data[[variantCol]]) %>%
        dplyr::mutate(totalVariantCountAll = sum(.data[[countCol]])) %>%
        dplyr::mutate(totalVariantCountFilt = sum(.data[[countCol]][.data$keep])) %>%
        dplyr::ungroup()
    plots[["totalBcCountAll"]] <- 
        ggplot(countTable %>% dplyr::filter(!duplicated(.data[[barcodeCol]])), 
               aes(x = .data$totalBarcodeCountAll)) + 
        geom_histogram(bins = 100, fill = "lightgrey") + 
        theme_bw() + scale_x_log10() + 
        labs(title = "Total barcode count (all)", x = "Total barcode count")
    plots[["totalBcCountFiltered"]] <- 
        ggplot(countTable %>% dplyr::filter(!duplicated(.data[[barcodeCol]])), 
               aes(x = .data$totalBarcodeCountFilt)) + 
        geom_histogram(bins = 100, fill = "lightgrey") + 
        theme_bw() + scale_x_log10() + 
        labs(title = "Total barcode count (filtered)", x = "Total barcode count")
    
    ## Calculate relative frequency for each pair as the ratio between the 
    ## count for the pair and the total count for the barcode or variant, 
    ## respectively
    countTable <- countTable %>%
        dplyr::mutate(fracOfTotalBarcodeCountFilt = 
                          .data[[countCol]]/.data$totalBarcodeCountFilt,
                      fracOfTotalVariantCountFilt = 
                          .data[[countCol]]/.data$totalVariantCountFilt)
    
    ## Order by decreasing frequency and keep only the top hit for each 
    ## barcode
    countTable <- countTable %>%
        dplyr::arrange(dplyr::desc(.data[[countCol]])) %>%
        dplyr::mutate(keep = .data$keep & !duplicated(.data[[barcodeCol]]))
    
    ## Flag pairs where the count/total barcode count is too low
    countTable <- countTable %>%
        dplyr::mutate(keep = .data$keep & 
                          (.data$fracOfTotalBarcodeCountFilt >= 
                               minFracOfBarcodeCount))
    
    ## Plot count for pair vs total count for barcode
    plots[["pairVsBarcodeCount"]] <- 
        ggplot(countTable, aes(x = .data$totalBarcodeCountFilt, 
                               y = .data[[countCol]], 
                               color = .data$keep)) + 
        geom_point(size = 0.5) + theme_bw() + 
        scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red"),
                           name = "Keep") + 
        scale_x_log10() + scale_y_log10() + 
        labs(x = "Total barcode count", y = "Barcode/variant pair count")
    
    return(list(counts = countTable, plots = plots))
}



