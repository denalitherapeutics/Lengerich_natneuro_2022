#' Plot a volcano plot both static and with plotly and produce a DT::datatable
#' @param DE - data frame with differential expression results
#' (eg: toptable output, DESeq2::results ouput, etc.)
#' @param fcColumn scalar character, the column with logFoldChanges in `DE`.
#' @param pColumn scalar character, the column with pvalues in `DE`.
#' @param geneCol scalar character, the column with gene names in `DE`.
#' @param fcMin Scalar numeric, the absolute FC cutoff to be used
#' (defaults to 0.5)
#' @param pMax calar numeric, The p-value cutoff to be used (defaults to 0.05)
#' @param labels Character vector of genes to be plotted, elements must match
#' the `geneCol` of `DE`.
#' @return A ggplot object
#' @import ggplot2
#' @import ggrepel
#' @importFrom checkmate assert_choice assert_character
#' assert_number
#' @export
#' @examples
#' \dontrun{
#' volcanoPlots(DE, fcColumn, pColumn, geneCol)
#' }
simpleVolcano <- function(DE, fcColumn, pColumn, geneCol, fcMin = 0.5,
                          pMax = 0.05, labels = NULL) {
  checkmate::assert_choice(fcColumn, choices = colnames(DE))
  checkmate::assert_choice(pColumn, choices = colnames(DE))
  checkmate::assert_choice(geneCol, choices = colnames(DE))
  checkmate::assert_number(fcMin)
  checkmate::assert_number(pMax, lower = 0, upper = 1)
  checkmate::assert_character(labels, null.ok = TRUE)

  DE$diffexpressed <- "NO"
  DE$diffexpressed[DE[,fcColumn] > fcMin & DE[,pColumn] < pMax] <- "UP"
  DE$diffexpressed[DE[,fcColumn] < -fcMin & DE[,pColumn] < pMax] <- "DOWN"

  DE$logp <- -log10(DE[,pColumn])
  DE$logfc <- DE[,fcColumn]
  DE$gene <- DE[,geneCol]

  p1 <- ggplot(data = DE, aes(x=logfc, y = logp, col = diffexpressed,
                              label = gene)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
    geom_vline(xintercept = c(-fcMin, fcMin), col = "gray40",
               linetype = "dashed") +
    geom_hline(yintercept = -log10(pMax), col = "gray40", linetype = "dashed") +
    theme(legend.position = "none")

  if(!is.null(labels)){
    p1 <- p1 + ggrepel::geom_label_repel(
      label = ifelse(DE$gene %in% labels, DE$gene, NA),
      box.padding = unit(0.25, "lines"),
      hjust = 1, max.overlaps = 100)
  }

  return(p1)
}
