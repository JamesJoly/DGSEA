#' Make mountain plots from a DGSEA list
#'
#' @param DGSEA.list A list resulting from either dgsea_targeted() or dgsea_untargeted()
#' @param Gene.Set.A A character string containing the first gene set you wish to plot
#' @param Gene.Set.B A character string containing the second gene set you wish to plot
#' @param color Set to TRUE or FALSE to generate enrichment plot with color or black and white.
#'
#' @return A ggplot object demonstrating the enrichment of two gene sets
#' @export
#'
#' @examples
#' expression_data <- seq(-1,1, by = 0.077)
#' Gene <- LETTERS[seq(from = 1, to = 26)]
#' data <- as.data.frame(cbind(Gene,expression_data))
#' gene.set.A <- c("A","B","C","D","E") #positive enrichment
#' gene.set.B <- c("V","W","X","Y","Z") #negative enrichment
#' gene.set.X <- c("J","Q","O","E","F")
#' gene.set.Y <- c("D","S","K","L","R")
#' gene.set.Z <- c("G","W","P","B","T")
#' names <- c("gene.set.A","gene.set.B","gene.set.X","gene.set.Y","gene.set.Z")
#' Letter.Sets <- list(gene.set.A,gene.set.B,gene.set.X,gene.set.Y,gene.set.Z)
#' names(Letter.Sets) <- NULL
#' Gene.Sets <- list(genesets = Letter.Sets, geneset.names = names)
#'
#' dgsea.results <- dgsea_untargeted(input.df = data, gmt.list = Gene.Sets)
#'
#' make_mountain_plots(DGSEA.list = dgsea.results,
#' Gene.Set.A = "gene.set.A", Gene.Set.B = "gene.set.B")
#'
make_mountain_plots <- function(DGSEA.list, Gene.Set.A, Gene.Set.B, color = TRUE){

  if (color == TRUE){
    color.palette <- c("red","blue")
  } else if (color == FALSE){
    color.palette <- c("black","grey70")
  }

  DGSEA.list <- DGSEA.list
  plotting.info <- DGSEA.list$Mountain.Plot.Info
  ES_profiles <- plotting.info$MountainPlot
  position.of.hits <- plotting.info$Position.of.hits

  ranking.metric <- as.numeric(DGSEA.list$ranking.metric)

  ES_GeneSetA <- ES_profiles[[Gene.Set.A]]
  ES_GeneSetB <- ES_profiles[[Gene.Set.B]]

  pos.GeneSetA <- position.of.hits[[Gene.Set.A]]
  pos.GeneSetB <- position.of.hits[[Gene.Set.B]]

  myName <- paste(Gene.Set.A, " - ", Gene.Set.B, sep = "")

  DGSEA.results <- DGSEA.list$DGSEA.Results
  DGSEA.results$p_value_AB <- as.numeric(as.character(DGSEA.results$p_value_AB))
  DGSEA.results$FDR <- as.numeric(as.character(DGSEA.results$FDR))

  gsea.normalized.enrichment.score <- signif(DGSEA.results[DGSEA.results$Gene.Sets.Compared == myName,]$NES_AB, digits =3)
  gsea.p.value <- signif(DGSEA.results[DGSEA.results$Gene.Sets.Compared == myName,]$p_value_AB, digits = 3)
  gsea.q.val <- signif(DGSEA.results[DGSEA.results$Gene.Sets.Compared == myName,]$FDR, digits = 3)

  results.as.text <- paste("NES",gsea.normalized.enrichment.score,"\np-value:", gsea.p.value,
                           "\nq-value:",gsea.q.val)

  gsea.hit.indices.A <- pos.GeneSetA
  gsea.hit.indices.B <- pos.GeneSetB

  gsea.es.profile.A <- ES_GeneSetA[gsea.hit.indices.A]
  gsea.es.profile.B <- ES_GeneSetB[gsea.hit.indices.B]

  ngenes <- length(ES_GeneSetA)
  combined.ES <- c(ES_GeneSetA, ES_GeneSetB)
  max.ES <- max(combined.ES)
  min.ES <- min(combined.ES)

  mtn.plot <- NULL
  mtn.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x=c(0,pos.GeneSetA,ngenes), y = c(0,gsea.es.profile.A,0)), color = color.palette[1]) +
    ggplot2::geom_line(ggplot2::aes(x=c(0,pos.GeneSetB,ngenes), y = c(0,gsea.es.profile.B,0)), color = color.palette[2]) +
    ggplot2::labs(x = NULL, y = "ES", title = myName) + ggplot2::scale_x_continuous(expand = c(0,0)) +
    #cowplot::theme_nothing() +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(0, "pt"),
      axis.title.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0,0,0,0),"cm"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (gsea.normalized.enrichment.score > 0){
    mtn.plot <- mtn.plot + ggplot2::annotate("text",label = results.as.text, x = ngenes * 0.90, y= max.ES * 0.70, color = "black")
  } else if (gsea.normalized.enrichment.score < 0){
    mtn.plot <- mtn.plot + ggplot2::annotate("text",label = results.as.text, x = ngenes * 0.10, y= max.ES * 0.30, color = "black")
  }

  hit.markers <- ggplot2::ggplot() +
    cowplot::theme_nothing() +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA)
    ) +
    ggplot2::geom_col(ggplot2::aes(x = pos.GeneSetA, y = 1), color = color.palette[1], fill = color.palette[1]) +
    ggplot2::geom_col(ggplot2::aes(x = pos.GeneSetB, y = 1), color = color.palette[2], fill = color.palette[2]) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    aplot::xlim2(mtn.plot)

  show.rank.metric <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = 1:ngenes, y = ranking.metric), color = "grey50") +
    ggplot2::geom_area(ggplot2::aes(x = 1:ngenes, y = ranking.metric), fill = "grey50") +
    ggplot2::labs(x = "Genes", y = "Rank Metric") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0,0,0,0),"cm")
    ) +
    aplot::xlim2(mtn.plot)



  assembled <- mtn.plot + hit.markers + show.rank.metric + patchwork::plot_layout(ncol = 1, heights = c(2.5,0.5,2))
  return(assembled)
}
