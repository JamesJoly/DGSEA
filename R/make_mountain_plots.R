#' Make mountain plots from a DGSEA list
#'
#' @param DGSEA.list A list resulting from either dgsea_targeted() or dgsea_untargeted()
#' @param Gene.Set.A A character string containing the first gene set you wish to plot
#' @param Gene.Set.B A character string containing the second gene set you wish to plot
#' @param color Set to TRUE or FALSE to generate enrichment plot with color or black and white.
#'
#' @return A plot demonstrating the enrichment of two gene sets
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
  requireNamespace("DGSEA")

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

  gsea.layout <- graphics::layout(matrix(c(1, 2, 3)), heights = c(2, 0.5, 2))
  p1 <- graphics::layout.show(gsea.layout)

  # Create plots
  p1 <- graphics::par(mar = c(0, 5, 2, 2))

  combined_ES <- c(ES_GeneSetA, ES_GeneSetB)
  enrichment.score.range <- c(min(combined_ES), max(combined_ES))

  gsea.hit.indices.A <- pos.GeneSetA
  gsea.hit.indices.B <- pos.GeneSetB

  gsea.es.profile.A <- ES_GeneSetA[gsea.hit.indices.A]
  gsea.es.profile.B <- ES_GeneSetB[gsea.hit.indices.B]

  p1 <- graphics::plot(c(1, gsea.hit.indices.A, length(ES_GeneSetA)),
                       c(0, gsea.es.profile.A, 0), type = "l", col = color.palette[1], lwd = 1.5, xaxt = "n",
                       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)", main = myName,
                       ylim = enrichment.score.range,
                       #main = list(gsea.gene.set, font = 1, cex = 1),
                       panel.first = {
                         graphics::abline(h = seq(round(enrichment.score.range[1], digits = 1),
                                                  enrichment.score.range[2], 0.1),
                                          col = "gray95", lty = 2)
                         graphics::abline(h = 0, col = "gray50", lty = 2)
                       }
  )
  p1 <- p1 + graphics::lines(c(1, gsea.hit.indices.B, length(ES_GeneSetB)),
                             c(0, gsea.es.profile.B, 0), type = "l", col = color.palette[2], lwd = 1.5, xaxt = "n",
                             xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
                             ylim = enrichment.score.range,
                             #main = list(gsea.gene.set, font = 1, cex = 1),
                             panel.first = {
                               graphics::abline(h = seq(round(enrichment.score.range[1], digits = 1),
                                                        enrichment.score.range[2], 0.1),
                                                col = "gray95", lty = 2)
                               graphics::abline(h = 0, col = "gray50", lty = 2)
                             }
  )

  plot.coordinates <- graphics::par("usr")

  p1 <- p1 + {
    if(gsea.normalized.enrichment.score < 0) {
      graphics::text(length(ES_GeneSetA) * 0.05, plot.coordinates[3] * 0.95,
                     paste("p-value:", gsea.p.value, "\nq-value:",
                           gsea.q.val, "\nNES:",
                           gsea.normalized.enrichment.score), adj = c(0, 0))
      graphics::text(length(ES_GeneSetA) * 0.95, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.02),
                     paste(Gene.Set.A), col = color.palette[1], adj = c(1,1))
      graphics::text(length(ES_GeneSetA) * 0.95, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.10),
                     paste(Gene.Set.B), col = color.palette[2], adj = c(1,1))

    } else {
      graphics::text(length(ES_GeneSetA) * 0.95, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
                     paste("p-value:", gsea.p.value, "\nq-value:",
                           gsea.q.val, "\n NES:",
                           gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
      graphics::text(length(ES_GeneSetA) * 0.05, plot.coordinates[3] * 0.90,
                     paste(Gene.Set.A), col =  color.palette[1], adj = c(0,0))
      graphics::text(length(ES_GeneSetA) * 0.05, plot.coordinates[3] * 0.6,
                     paste(Gene.Set.B), col = color.palette[2], adj = c(0,0))
    }
  }

  p1 <- p1 + {
    graphics::par(mar = c(0, 5, 0, 2))
    graphics::plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
                   ylab = "", xlim = c(1, length(ES_GeneSetA)))
    graphics::abline(v = gsea.hit.indices.A, lwd = 0.75, col = color.palette[1])
    graphics::abline(v = gsea.hit.indices.B, lwd = 0.75, col = color.palette[2])
  }
  ## Need to figure out how to save Rank metric from DGSEA for the rest of this.


  metric.range <- c(min(ranking.metric),max(ranking.metric))

  graphics::par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(ranking.metric, digits = 2))
  graphics::plot(ranking.metric, type = "n", xaxs = "i",
                 xlab = "Rank in ordered gene list", xlim = c(0, length(ranking.metric)),
                 ylim = metric.range, yaxs = "i",
                 ylab = "Ranking Metric",
                 panel.first = graphics::abline(h = seq(metric.range[1] / 2,
                                                        metric.range[2] - metric.range[1] / 4,
                                                        metric.range[2] / 2), col = "gray95", lty = 2))

  graphics::barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
                    xlab = "Rank in ordered gene list", xlim = c(0, length(ranking.metric)),
                    ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
                    #ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric),
                    space = 0, add = TRUE)
  graphics::box()
  p2 <- grDevices::recordPlot()
  return(p2)
  #dev.print(pdf,file = paste(Samples[u],"DGSEA.pdf"), height = 4, width = 4)
}
