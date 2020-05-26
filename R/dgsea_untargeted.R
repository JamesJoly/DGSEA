#' This function performs untargeted DGSEA comparing all gene sets to each other
#'
#' @param input.df Data frame containing Genes in column 1 and ranking metrics for samples (e.g signal-to-noise ratio) in column 2
#' @param gmt.list gmt list containing the gene sets the user wishes to use with gene set name as column name and genes in the rows (i.e. from GSA.read.gmt())
#' @param num.permutations Number of permutations to perform, default = 1000
#' @param stat.type Character string set to "Weighted" (weight = 1) or "Classic" (score weight = 0), default is "Weighted"
#'
#' @return A list containing DGSEA results, GSEA results, and information to make mountain plots using the make_mountain_plots function
#' @export
#'
#' @examples
#' expression_data <- seq(-1,1, by = 0.077)
#' Gene <- LETTERS[seq(from = 1, to = 26)]
#' data <- as.data.frame(cbind(Gene,expression_data))
#' gene.set.A <- c("A","B","C","D","E")
#' gene.set.B <- c("V","W","X","Y","Z")
#' gene.set.X <- c("J","Q","O","E","F")
#' gene.set.Y <- c("D","S","K","L","R")
#' gene.set.Z <- c("G","W","P","B","T")
#' names <- c("gene.set.A","gene.set.B","gene.set.X","gene.set.Y","gene.set.Z")
#' Letter.Sets <- list(gene.set.A,gene.set.B,gene.set.X,gene.set.Y,gene.set.Z)
#' names(Letter.Sets) <- NULL
#' Gene.Sets <- list(genesets = Letter.Sets, geneset.names = names)
#'
#' dgsea_untargeted(input.df = data, gmt.list = Gene.Sets)
#'

dgsea_untargeted <- function(input.df, gmt.list,
                             num.permutations = 1000,
                             stat.type = "Weighted"){
  nperm = num.permutations #number of permutations
  #utils::globalVariables(c("x"), add = F)
  if (stat.type == "Classic"){
    score.weight = 0
  }
  if (stat.type == "Weighted"){
    score.weight = 1
  }

  #Read in gene expression data
  #Genes should be first column, named "Gene"
  #Samples should be columns 2:N
  data_in <- input.df

  gmt.for.reformat <- gmt.list
  Gene.Sets <- t(plyr::ldply(gmt.for.reformat$genesets, rbind)) #reformat gmt list to desired format
  colnames(Gene.Sets) <- gmt.for.reformat$geneset.names

  Gene.Sets <- as.data.frame(Gene.Sets)

  testthat::expect_is(data_in, "data.frame")
  testthat::expect_is(Gene.Sets, "data.frame")

  GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL){
    tag.indicators <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicators
    N <- length(gene.list)
    Nh <- numhits_pathway
    Nm <- N - Nh
    if (weighted.score.type == 0){
      correl.vector <- rep(1,N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicators == 1])
    norm.tag <- 1.0/sum.correl.tag
    norm.no.tag <- 1.0/Nm
    RES <- cumsum(tag.indicators * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > - min.ES) {
      #      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      #      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicators))
  } #for real ES

  GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL) {

    N <- length(gene.list)
    Nh <- numhits_pathway
    Nm <-  N - Nh

    loc.vector <- vector(length=N, mode="numeric")
    peak.res.vector <- vector(length=Nh, mode="numeric")
    valley.res.vector <- vector(length=Nh, mode="numeric")
    tag.correl.vector <- vector(length=Nh, mode="numeric")
    tag.diff.vector <- vector(length=Nh, mode="numeric")
    tag.loc.vector <- vector(length=Nh, mode="numeric")

    loc.vector[gene.list] <- seq(1, N)
    tag.loc.vector <- loc.vector[gene.set]

    tag.loc.vector <- sort(tag.loc.vector, decreasing = F)

    if (weighted.score.type == 0) {
      tag.correl.vector <- rep(1, Nh)
    } else if (weighted.score.type == 1) {
      tag.correl.vector <- correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
    } else if (weighted.score.type == 2) {
      tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
    } else {
      tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
    }

    norm.tag <- 1.0/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1.0/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)

    return(ES)

  } #for permutation ES

  Samples <- colnames(data_in)
  if (Samples[1] != "Gene"){
    stop("Please ensure that your data frame is organized with the first column to be named 'Gene'")
  }
  Samples <- Samples[-1]

  Gene.Sets.All <- colnames(Gene.Sets)

  if (length(Gene.Sets.All) > 150){
    print("Warning: Performing untargeted DGSEA with a large number of gene sets (>150) greatly increases computation time.")
  }

  annotations <- matrix(data = 0, nrow = nrow(data_in), ncol = length(Gene.Sets.All))
  colnames(annotations) <- Gene.Sets.All

  annotations <- as.data.frame(annotations)

  annotations <- cbind(data_in$Gene,annotations)
  colnames(annotations) <- c("Gene", Gene.Sets.All)
  annotations <- as.matrix(annotations)

  num.hits.pathways <- list()

  ### Annotate gene sets
  print("Annotating gene sets...")

  for (j in 1:length(Gene.Sets.All)){
    temp.pathway <- Gene.Sets[,Gene.Sets.All[j]]
    for (i in 1:nrow(annotations)){
      if (annotations[i,"Gene"] %in% temp.pathway){
        annotations[i,j+1] = "X";
      }
    }
    num.hits.pathways[[Gene.Sets.All[j]]] <- sum(annotations[,Gene.Sets.All[j]] == "X")
  }

  num.hits.pathways.df <- matrix(unlist(num.hits.pathways))
  row.names(num.hits.pathways.df) = Gene.Sets.All
  num.gene.sets.under.5 <- which(num.hits.pathways.df < 5)
  if (length(num.gene.sets.under.5) > 1){
    print("Warning: Removing gene sets with less than 5 genes observed in data set.")
    gene.sets.to.remove <- Gene.Sets.All[num.gene.sets.under.5]
    annotations[,which(colnames(annotations) %in% gene.sets.to.remove)] <- NULL
  }
  annotations <- as.data.frame(annotations)
  data_in <- merge(data_in, annotations, by = "Gene")

  data_in <- stats::na.omit(data_in)

  DGSEA.Results.All.Samples <- matrix(data = NA, nrow = 0, ncol = 13)
  colnames(DGSEA.Results.All.Samples) <- c("Gene Set A", "Gene Set B",
                                           "ES A", "NES A","p-value A", "FDR A",
                                           "ES B", "NES B","p-value B", "FDR B",
                                           "ES AB", "NES AB","p-value AB")

  GSEA.Results.All.Samples <- matrix(data = NA, nrow = 0, ncol = 7)
  colnames(GSEA.Results.All.Samples) <- c("Sample","Gene.Set","KS","KS_Normalized",
                                          "p-value","Position at Max",
                                          "FDR q-value")
  Mountain.Plot.Info.All.Samples <- list()
  rank_metric.All.Samples <- list()

  #Find out how many cores are available (if you don't already know)
  cores<-parallel::detectCores()
  #Create cluster with desired number of cores, leave one open for the machine
  #core processes
  cl <- snow::makeCluster(cores[1]-1)
  #Register cluster
  doSNOW::registerDoSNOW(cl)

  Samples <- Samples[1]

  rm(annotations)
  data_in2 <- array(data = NA)
  for (u in 1:length(Samples)){
    loop.time <- Sys.time()

    data_in2 <- cbind(subset(data_in, select = Gene.Sets.All),
                      dplyr::select(data_in, Samples[u]))  #select one Sample type and the genes and Gene.Sets.A.and.B
    data_in2[,Samples[u]] <- as.numeric(as.character(data_in2[,Samples[u]]))
    data_in2 <- data_in2[order(-data_in2[,Samples[u]]),] #sort by descending order for the rank metric
    rownames(data_in2) <- 1:nrow(data_in2) #reorder row indices for counting in for loop below

    ## Assuming first two columns in data table are Genes and Rank Metric (e.g. Foldchange, SNR)

    GSEA.Results <- matrix(data = NA, nrow = length(Gene.Sets.All), ncol = 7)
    colnames(GSEA.Results) <- c("Sample","Gene.Set","KS","KS_Normalized",
                                "p_value","Position_at_max",
                                "FDR_q_value")
    GSEA.Results <- as.data.frame(GSEA.Results)
    GSEA.Results$Gene.Set <- Gene.Sets.All
    GSEA.Results$Sample <- Samples[u]

    ions <- nrow(data_in2)

    #for plotting
    ks_results_plot <- list()
    positions.of.hits <- list()

    #ks_results_plot <- as.data.frame(ks_results_plot)
    gene.list <- 1:ions
    rank_metric <- data_in2[,Samples[u]] #Save the rank metric

    pos_gene_set <- array(data = 0, dim = nrow(data_in2), dimnames = NULL);

    ## Calculate Real KS Statistic
    for (i in 1:length(Gene.Sets.All)){
      data_in3 <- data_in2[,Gene.Sets.All[i]]
      numhits_pathway <- sum(data_in3 == "X"); #check to see if there is anything in the column (e.g. X)
      if (numhits_pathway > 1){
        pos_gene_set <- which(data_in2[,Gene.Sets.All[i]] %in% c("X"))
        KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$KS <- KS_real$ES;
        GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.All[i],]$Position_at_max <- KS_real$arg.ES;
        ks_results_plot[[Gene.Sets.All[i]]] = KS_real$RES
        positions.of.hits[[Gene.Sets.All[i]]] = pos_gene_set
      }
    }

    Mountain.Plot.Info <- list(MountainPlot = ks_results_plot, Position.of.hits = positions.of.hits)
    rm(pos_gene_set)
    rm(numhits_pathway)
    rm(data_in3)
    rm(KS_real)

    print("Calculating permutations...")

    pb <- utils::txtProgressBar(max = num.permutations, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    KSRandomArray <- matrix(data = NA, nrow = nperm, ncol = length(Gene.Sets.All))
    num.gene.sets.all <- length(Gene.Sets.All)
    `%dopar%` <- foreach::`%dopar%`
    KSRandomArray <- foreach::foreach(L = 1:nperm, .combine = "rbind",.options.snow = opts) %dopar% {
      temp.KSRandomArray <- matrix(data = NA, nrow = 1, ncol = num.gene.sets.all)
      for(i in 1:length(Gene.Sets.All)){
        numhits_pathway <- length(positions.of.hits[[Gene.Sets.All[i]]])
        pos_gene_set <- sample(1:ions,numhits_pathway)
        temp.KSRandomArray[,i] <- GSEA.EnrichmentScore2(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
      }
      temp.KSRandomArray
    }
    colnames(KSRandomArray) <- Gene.Sets.All

    rm(opts)
    rm(pb)
    KSRandomArray <- data.frame(matrix(unlist(KSRandomArray), nrow = nperm, byrow = T))
    colnames(KSRandomArray) <- Gene.Sets.All
    KSRandomArray <- stats::na.omit(KSRandomArray)

    print("Normalizing enrichment scores...")
    KSRandomArray <- as.data.frame(KSRandomArray)
    ###normalize the GSEA distribution
    KSRandomArray.Norm <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = ncol(KSRandomArray))
    colnames(KSRandomArray.Norm) <- colnames(KSRandomArray)
    avg <- 0
    KSRandomArray.temp <- 0
    for (i in 1:ncol(KSRandomArray.Norm)){
      avg <- 0
      KSRandomArray.temp <- KSRandomArray[,i]
      pos.temp <- KSRandomArray.temp[which(KSRandomArray.temp >= 0)]
      neg.temp <- KSRandomArray.temp[which(KSRandomArray.temp < 0)]

      avg.pos <- mean(pos.temp)
      avg.neg <- mean(neg.temp)

      norm.pos.temp <- pos.temp / avg.pos
      norm.neg.temp <- neg.temp / avg.neg * -1

      norm.perms <- c(norm.pos.temp,norm.neg.temp)

      KSRandomArray.Norm[,i] <- norm.perms

    }

    GSEA.NES.perms <- as.vector(KSRandomArray.Norm)
    rm(KSRandomArray.Norm)
    GSEA.NES.perms.pos <- GSEA.NES.perms[which(GSEA.NES.perms >= 0)]
    GSEA.NES.perms.neg <- GSEA.NES.perms[which(GSEA.NES.perms < 0)]
    rm(GSEA.NES.perms)
    percent.pos.GSEA <- sum(GSEA.Results$KS > 0) / length(GSEA.Results$KS)
    percent.neg.GSEA <- sum(GSEA.Results$KS < 0) / length(GSEA.Results$KS)

    # Calculate GSEA NES and p-value and FDR
    print("Calculating GSEA FDR...")
    for (i in 1:length(Gene.Sets.All)){
      temp.gene.set <- Gene.Sets.All[i]
      temp.KS <- GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS
      if (temp.KS > 0){
        pos.perms <- KSRandomArray[,temp.gene.set]
        pos.perms <- pos.perms[which(pos.perms > 0)]
        #p-val
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(pos.perms > temp.KS) / length(pos.perms),digits = 3)
        #NES
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized = signif(temp.KS / mean(pos.perms), digits = 3)
        #FDR
        percent.temp <- sum(GSEA.NES.perms.pos > GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized) / length(GSEA.NES.perms.pos)
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.pos.GSEA, digits = 3) < 1, signif(percent.temp / percent.pos.GSEA, digits = 3), 1)
      } else if (temp.KS < 0){
        neg.perms <- KSRandomArray[,temp.gene.set]
        neg.perms <- neg.perms[which(neg.perms < 0)]
        #p-val
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$p_value = signif(sum(neg.perms < temp.KS) / length(neg.perms),digits = 3)
        #NES
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized = signif(temp.KS / mean(neg.perms) * -1, digits = 3)
        #FDR
        percent.temp <- sum(GSEA.NES.perms.neg < GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$KS_Normalized) / length(GSEA.NES.perms.neg)
        GSEA.Results[GSEA.Results$Gene.Set == temp.gene.set,]$FDR_q_value = ifelse(signif(percent.temp / percent.neg.GSEA, digits = 3) < 1, signif(percent.temp / percent.neg.GSEA, digits = 3), 1)
      }
    }
    rm(GSEA.NES.perms.pos)
    rm(GSEA.NES.perms.neg)

    print("Generating DGSEA Null Distribution...")
    Permutations.Results <- list()

    combinations <- utils::combn(Gene.Sets.All, 2)

    Permutations.Results <- lapply(1:ncol(combinations), function(x){
      tempX <- combinations[1,x]
      tempY <- combinations[2,x]
      KSRandomArray[,tempX] - KSRandomArray[,tempY]
    })
    myNames <- c(paste(combinations[1,]," - ",combinations[2,], sep = ""))
    names(Permutations.Results) <- myNames
    rm(KSRandomArray)


    list.DGSEA.Results <- list()
    list.DGSEA.Results <- lapply(1:ncol(combinations), function(x){
      tempX <- combinations[1,x]
      tempY <- combinations[2,x]
      temp.Name <- paste(combinations[1,x]," - ",combinations[2,x], sep = "")

      temp.KS.A <- GSEA.Results[GSEA.Results$Gene.Set == tempX,]$KS
      temp.KS.B <- GSEA.Results[GSEA.Results$Gene.Set == tempY,]$KS

      temp.KS_Norm.A <- GSEA.Results[GSEA.Results$Gene.Set == tempX,]$KS_Normalized
      temp.KS_Norm.B <- GSEA.Results[GSEA.Results$Gene.Set == tempY,]$KS_Normalized

      temp.p.value.A <- GSEA.Results[GSEA.Results$Gene.Set == tempX,]$p_value
      temp.p.value.B <- GSEA.Results[GSEA.Results$Gene.Set == tempY,]$p_value

      temp.FDR.A <- GSEA.Results[GSEA.Results$Gene.Set == tempX,]$FDR_q_value
      temp.FDR.B <- GSEA.Results[GSEA.Results$Gene.Set == tempY,]$FDR_q_value

      permutations <- as.matrix(unlist(Permutations.Results[[temp.Name]]))
      permutations.pos <- permutations[which(permutations > 0)]
      permutations.neg <- permutations[which(permutations < 0)]

      temp.KS.AB <- temp.KS.A - temp.KS.B

      if (temp.KS.AB > 0 && length(permutations.pos) > 0){
        temp.AB.pval <- signif(sum(permutations.pos > temp.KS.AB ) / length(permutations.pos), digits = 3)
        temp.AB.KS.Norm <-  signif(temp.KS.AB / mean(permutations.pos), digits = 3)
      } else if (temp.KS.AB < 0 && length(permutations.neg) > 0){
        temp.AB.pval <- signif(sum(permutations.neg < temp.KS.AB) / length(permutations.neg), digits = 3)
        temp.AB.KS.Norm <-  signif(temp.KS.AB / mean(permutations.neg) * -1, digits = 3)
      } else {
        temp.AB.pval <- NA
        temp.AB.KS.Norm <- NA
      }
      return(list(Gene.Set.A = tempX, Gene.Set.B = tempY,
                  KS.A = temp.KS.A, KS_Normalized_A = temp.KS_Norm.A,p.value.A = temp.p.value.A, FDR.A = temp.FDR.A,
                  KS.B = temp.KS.B, KS_Normalized_B = temp.KS_Norm.B,p.value.B = temp.p.value.B, FDR.A = temp.FDR.B,
                  KS.AB = temp.KS.AB, KS.Norm.AB = temp.AB.KS.Norm, p.value.AB = temp.AB.pval))
    })
    names(list.DGSEA.Results) <- myNames

    DGSEA.Results <- data.frame(matrix(unlist(list.DGSEA.Results), nrow = ncol(combinations), byrow = T))
    #DGSEA.Results <- matrix(data = NA, nrow = 0, ncol = 13)
    colnames(DGSEA.Results) <- c("Gene Set A", "Gene Set B",
                                 "ES A","NES A","p-value A", "FDR A",
                                 "ES B","NES B","p-value B", "FDR B",
                                 "ES AB", "NES_AB","p_value_AB")

    DGSEA.Results$Sample <- Samples[u]
    DGSEA.Results <- DGSEA.Results[,c(14,1:13)]
    rm(list.DGSEA.Results)

    ### Control for false discovery rate

    list.Permutations.Normazlized <- list()

    list.Permutations.Normazlized <- lapply(1:ncol(combinations), function(x){
      temp.Name <- paste(combinations[1,x]," - ",combinations[2,x], sep = "")
      temp.permutations <- Permutations.Results[[temp.Name]]
      #Permutations.Results[[temp.Name]] <- NULL
      temp.permutations.pos <- temp.permutations[which(temp.permutations >= 0)]
      temp.permutations.neg <- temp.permutations[which(temp.permutations < 0)]
      temp.norm.perms.pos <- temp.permutations.pos / mean(temp.permutations.pos)
      temp.norm.perms.neg <- temp.permutations.neg / mean(temp.permutations.neg) * -1
      c(temp.norm.perms.pos,temp.norm.perms.neg)
    })
    Norm.Perms.Vector <- as.vector(unlist(list.Permutations.Normazlized))
    rm(Permutations.Results)

    Norm.Perms.pos <- Norm.Perms.Vector[which(Norm.Perms.Vector >= 0)]
    Norm.Perms.neg <- Norm.Perms.Vector[which(Norm.Perms.Vector < 0)]
    rm(Norm.Perms.Vector)


    DGSEA.Results$NES_AB <- as.numeric(as.character(DGSEA.Results$NES_AB))
    DGSEA.Results$p_value_AB <- as.numeric(as.character(DGSEA.Results$p_value_AB))


    real.DGSEA.NES.pos <- DGSEA.Results[which(DGSEA.Results$NES_AB >= 0),]$NES_AB
    real.DGSEA.NES.neg <- DGSEA.Results[which(DGSEA.Results$NES_AB < 0),]$NES_AB
    percent.DGSEA.NES.pos <- length(real.DGSEA.NES.pos) / length(DGSEA.Results$NES_AB)
    percent.DGSEA.NES.neg <- length(real.DGSEA.NES.neg) / length(DGSEA.Results$NES_AB)

    rm(real.DGSEA.NES.pos)
    rm(real.DGSEA.NES.neg)


    rm(list.Permutations.Normazlized)
    DGSEA.Results$Gene.Sets.Compared <- myNames


    num.pos.perms <- length(Norm.Perms.pos)
    num.neg.perms <- length(Norm.Perms.neg)

    print("Calculating DGSEA FDR...")
    filtered.combinations <- DGSEA.Results[DGSEA.Results$p_value_AB < 0.25,]$Gene.Sets.Compared
    filtered.combinations <- na.omit(filtered.combinations)
    opts <- list(progress = progress)
    pb <- utils::txtProgressBar(min = 0, max = length(filtered.combinations), style = 3)


    list.DGSEA.FDR <- list()
    list.DGSEA.FDR <- foreach::foreach(x = 1:length(filtered.combinations), .combine = "rbind",.options.snow = opts) %dopar% {
      utils::setTxtProgressBar(pb, x)
      temp.Name <- filtered.combinations[x]
      temp.NES <- DGSEA.Results[DGSEA.Results$Gene.Sets.Compared == temp.Name,]$NES_AB
      if (is.numeric(temp.NES)){
        temp.FDR <- ifelse(temp.NES > 0,
                           ifelse(((sum(Norm.Perms.pos > temp.NES) / num.pos.perms)/ percent.DGSEA.NES.pos) < 1,
                                  signif(sum(Norm.Perms.pos > temp.NES) / num.pos.perms / percent.DGSEA.NES.pos, digits = 3),1),
                           ifelse((sum(Norm.Perms.neg < temp.NES) / num.neg.perms / percent.DGSEA.NES.neg) < 1,
                                  signif(sum(Norm.Perms.neg < temp.NES) / num.neg.perms / percent.DGSEA.NES.neg, digits = 3),1))
        FDR = temp.FDR
      }
    }
    rm(opts)
    rm(pb)

    list.DGSEA.FDR <- as.data.frame(list.DGSEA.FDR)
    row.names(list.DGSEA.FDR) <- filtered.combinations
    list.DGSEA.FDR$myName <- filtered.combinations
    colnames(list.DGSEA.FDR) <- c("FDR","myName")

    #DGSEA.FDR <- as.data.frame(as.matrix(unlist(list.DGSEA.FDR)))
    DGSEA.Results$FDR <- 1
    for (i in 1:length(filtered.combinations)){
      DGSEA.Results[DGSEA.Results$Gene.Sets.Compared == filtered.combinations[i],]$FDR <-
        list.DGSEA.FDR[list.DGSEA.FDR$myName == filtered.combinations[i],]$FDR
    }


    rm(Norm.Perms.pos)
    rm(Norm.Perms.neg)

    DGSEA.Results.All.Samples <- rbind(DGSEA.Results.All.Samples,DGSEA.Results)
    GSEA.Results.All.Samples <- rbind(GSEA.Results.All.Samples,GSEA.Results)
    Mountain.Plot.Info.All.Samples <- c(Mountain.Plot.Info.All.Samples,Mountain.Plot.Info)
    rank_metric.All.Samples <- c(rank_metric.All.Samples,rank_metric)

  }

  snow::stopCluster(cl)
  rm(cl)

  return(list(DGSEA.Results = DGSEA.Results.All.Samples,
              GSEA.Results = GSEA.Results.All.Samples,
              Mountain.Plot.Info = Mountain.Plot.Info.All.Samples,
              ranking.metric = rank_metric.All.Samples))

}
