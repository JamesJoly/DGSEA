#' This function performs targeted DGSEA comparing Gene Set A to Gene Set B while using background gene sets to control for false discovery rate
#'
#' @param input.df Data frame containing Genes in column 1 and ranking metrics for samples (e.g signal-to-noise ratio) in column 2
#' @param gmt.list gmt list containing the gene sets the user wishes to use with gene set name as column name and genes in the rows (i.e. from GSA.read.gmt())
#' @param Gene.Set.A.Name Character string containing the exact name of the first gene set the user wishes to compare
#' @param Gene.Set.B.Name Character string containing the exact name of the second gene set the user wishes to compare
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
#' dgsea_targeted(input.df = data,
#' gmt.list = Gene.Sets,
#' Gene.Set.A.Name = names[1], Gene.Set.B.Name = names[2])
#'
dgsea_targeted <- function(input.df, gmt.list, Gene.Set.A.Name, Gene.Set.B.Name,
                           num.permutations = 1000, stat.type = "Weighted"){
  nperm = num.permutations #number of permutations

  if (stat.type == "Classic"){
    score.weight = 0
  }
  if (stat.type == "Weighted"){
    score.weight = 1
  }

  requireNamespace("DGSEA")

  #Genes should be first column, named "Gene"
  #Samples should be columns 2:N

  data_in <- input.df

  gmt.for.reformat <- gmt.list
  Gene.Sets <- t(plyr::ldply(gmt.for.reformat$genesets, rbind)) #reformat gmt list to desired format
  colnames(Gene.Sets) <- gmt.for.reformat$geneset.names

  testthat::expect_is(data_in, "data.frame")

  colnames(data_in)[1] <- "Gene"
  expected.number.of.genes <- length(data_in$Gene)
  actual.number.of.genes <- length(unique(data_in$Gene))
  if (actual.number.of.genes < expected.number.of.genes){
    stop("Your gene list has duplicates, please collapse your data or process in a different way to use DGSEA.")
  }


  Gene.Sets.All <- as.data.frame(Gene.Sets)

  Gene.Set.A.subset <- as.character(Gene.Sets.All[,Gene.Set.A.Name])
  Gene.Set.B.subset <- as.character(Gene.Sets.All[,Gene.Set.B.Name])

  DGSEA.Gene.Sets <- as.data.frame(cbind(Gene.Set.A.subset,Gene.Set.B.subset))
  colnames(DGSEA.Gene.Sets) <- c(Gene.Set.A.Name,Gene.Set.B.Name)
  testthat::expect_is(DGSEA.Gene.Sets, "data.frame")


  background.Gene.Sets.A.and.B <- as.data.frame(Gene.Sets.All)
  background.Gene.Sets.A.and.B[,Gene.Set.A.Name] <- NULL #remove gene sets A and B from background to not get 0s later on
  background.Gene.Sets.A.and.B[,Gene.Set.B.Name] <- NULL
  testthat::expect_is(background.Gene.Sets.A.and.B, "data.frame")




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

  Gene.Sets.A.and.B <-colnames(DGSEA.Gene.Sets)

  ### Annotate gene sets
  data_in$Gene.Set.A <- 0
  data_in$Gene.Set.B <- 0


  #Gene.Set.B (column 2 in DGSEA.Gene.Sets)
  for (i in 1:nrow(data_in)){
    if (data_in$Gene[i] %in% DGSEA.Gene.Sets[,2]){
      data_in$Gene.Set.B[i] = "X";
    }
  }
  #Gene.Set.A (column 1 in DGSEA.Gene.Sets)
  for (i in 1:nrow(data_in)){
    if (data_in$Gene[i] %in% DGSEA.Gene.Sets[,1]){
      data_in$Gene.Set.A[i] = "X";
    }
  }

  colnames(data_in) <- gsub("Gene.Set.A",Gene.Sets.A.and.B[1], colnames(data_in))
  colnames(data_in) <- gsub("Gene.Set.B",Gene.Sets.A.and.B[2], colnames(data_in))

  background.Gene.Sets.A.and.B_unique <- colnames(background.Gene.Sets.A.and.B)

  annotations <- matrix(data = 0, nrow = nrow(data_in), ncol = length(background.Gene.Sets.A.and.B_unique))
  colnames(annotations) <- background.Gene.Sets.A.and.B_unique

  annotations <- as.data.frame(annotations)

  annotations <- cbind(data_in$Gene,annotations)
  colnames(annotations) <- c("Gene", background.Gene.Sets.A.and.B_unique)
  annotations <- as.matrix(annotations)
  num.hits.pathways <- list()

  ### Annotate gene sets
  print("Annotating gene sets...")

  for (j in 1:length(background.Gene.Sets.A.and.B_unique)){
    temp.pathway <- background.Gene.Sets.A.and.B[,background.Gene.Sets.A.and.B_unique[j]]
    for (i in 1:nrow(annotations)){
      if (annotations[i,"Gene"] %in% temp.pathway){
        annotations[i,j+1] = "X";
      }
    }
    num.hits.pathways[[background.Gene.Sets.A.and.B_unique[j]]] <- sum(annotations[,background.Gene.Sets.A.and.B_unique[j]] == "X")
  }
  num.hits.pathways.df <- matrix(unlist(num.hits.pathways))
  row.names(num.hits.pathways.df) = background.Gene.Sets.A.and.B_unique
  num.gene.sets.under.5 <- which(num.hits.pathways.df < 5)
  if (length(num.gene.sets.under.5) > 1){
    print("Warning: Removing gene sets with less than 5 genes observed in data set.")
    gene.sets.to.remove <- background.Gene.Sets.A.and.B_unique[num.gene.sets.under.5]
    annotations <- annotations[,-which(colnames(annotations) %in% gene.sets.to.remove)]
  }

  gene.sets.updated <- colnames(annotations)[-1]
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

  Gene.Sets.All <- c(Gene.Sets.A.and.B,gene.sets.updated)

  rm(annotations)
  data_in2 <- array(data = NA)

  #Find out how many cores are available (if you don't already know)
  cores<-parallel::detectCores()
  #Create cluster with desired number of cores, leave one open for the machine
  #core processes
  cl <- snow::makeCluster(cores[1]-1)
  #Register cluster
  doSNOW::registerDoSNOW(cl)

  Samples <- Samples[1]
  for (u in 1:length(Samples)){
    #loop.time <- Sys.time()

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
    gene.list <- 1:ions
    rank_metric <- data_in2[,Samples[u]] #Save the rank metric

    pos_gene_set <- array(data = 0, dim = nrow(data_in2), dimnames = NULL);

    #for plotting
    ks_results_plot <- list()
    positions.of.hits <- list()

    ## Calculate Real KS Statistic
    for (i in 1:length(Gene.Sets.All)){
      data_in3 <- data_in2[,Gene.Sets.All[i]]
      numhits_pathway <- sum(data_in3 == "X"); #check to see if there is anything in the column (e.g. X)
      if (numhits_pathway > 1){
        pos_gene_set <- which(data_in2[,Gene.Sets.All[i]] %in% c("X"))

        #Calculate Real KS statistic and save the maximum (pos or neg) as well as the position of the maximum to check
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

    ## Permutations using the same random shuffled gene order for each pathway

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

    print(" ")
    print("Normalizing scores...")
    KSRandomArray <- as.data.frame(KSRandomArray)
    ###Normalize the GSEA distribution
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
    #rm(KSRandomArray.Norm)
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

    ###Make a DGSEA KSRandomArray for Glyc - All Gene.Sets
    print("Generating DGSEA Null Distribution...")

    KSRandomArray$DGSEA <- KSRandomArray[,Gene.Sets.A.and.B[1]] - KSRandomArray[,Gene.Sets.A.and.B[2]]

    KSRandomArray.Norm <- as.data.frame(KSRandomArray.Norm)
    pos.DGSEA.perms <- KSRandomArray[which(KSRandomArray$DGSEA >= 0),]$DGSEA
    neg.DGSEA.perms <- KSRandomArray[which(KSRandomArray$DGSEA < 0),]$DGSEA

    norm.pos.DGSEA.perms <- pos.DGSEA.perms / mean(pos.DGSEA.perms)
    norm.neg.DGSEA.perms <- neg.DGSEA.perms / mean(neg.DGSEA.perms) * -1

    norm.DGSEA.perms <- c(norm.pos.DGSEA.perms,norm.neg.DGSEA.perms)

    KSRandomArray.Norm$DGSEA <- norm.DGSEA.perms

    DGSEA_RandomArray_A <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = length(background.Gene.Sets.A.and.B_unique))
    DGSEA_RandomArray_A <- as.data.frame(DGSEA_RandomArray_A)
    DGSEA_RandomArray_B <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = length(background.Gene.Sets.A.and.B_unique))
    DGSEA_RandomArray_B <- as.data.frame(DGSEA_RandomArray_B)

    DGSEA_Real_A <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(background.Gene.Sets.A.and.B_unique)))

    DGSEA_Real_B <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(background.Gene.Sets.A.and.B_unique)))

    for(i in 1:length(background.Gene.Sets.A.and.B_unique)){
      identifier <- paste("DGSEA A", background.Gene.Sets.A.and.B_unique[i])
      colnames(DGSEA_RandomArray_A)[i] <- identifier
      DGSEA_RandomArray_A[,i] <- KSRandomArray[,Gene.Sets.A.and.B[1]] - KSRandomArray[,background.Gene.Sets.A.and.B_unique[i]]
      colnames(DGSEA_Real_A)[i] <- identifier
      DGSEA_Real_A[,i] <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[1],]$KS -GSEA.Results[GSEA.Results$Gene.Set == background.Gene.Sets.A.and.B_unique[i],]$KS

      identifier <- paste("DGSEA B", background.Gene.Sets.A.and.B_unique[i])
      colnames(DGSEA_RandomArray_B)[i] <- identifier
      DGSEA_RandomArray_B[,i] <- KSRandomArray[,Gene.Sets.A.and.B[2]] - KSRandomArray[,background.Gene.Sets.A.and.B_unique[i]]
      colnames(DGSEA_Real_B)[i] <- identifier
      DGSEA_Real_B[,i] <-  GSEA.Results[GSEA.Results$Gene.Set == background.Gene.Sets.A.and.B_unique[i],]$KS - GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[2],]$KS
    }

    DGSEA_Real <- cbind(DGSEA_Real_A,DGSEA_Real_B)
    DGSEA_RandomArray <- cbind(DGSEA_RandomArray_A, DGSEA_RandomArray_B)

    #normalize the DGSEA distribution
    DGSEA_Background <- colnames(DGSEA_RandomArray)
    DGSEA_Background.norm <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = length(DGSEA_Background))
    colnames(DGSEA_Background.norm) <- DGSEA_Background
    Real.DGSEA_Background.NES <- matrix(data = NA, nrow = 1, ncol = length(DGSEA_Background))
    colnames(Real.DGSEA_Background.NES) <- DGSEA_Background
    for (i in 1:length(DGSEA_Background)){
      avg <- 0

      DGSEA_RandomArray.temp <- DGSEA_RandomArray[,DGSEA_Background[i]]
      pos.temp <- DGSEA_RandomArray.temp[which(DGSEA_RandomArray.temp >= 0)]
      neg.temp <- DGSEA_RandomArray.temp[which(DGSEA_RandomArray.temp < 0)]

      avg.pos <- mean(pos.temp)
      avg.neg <- mean(neg.temp)

      norm.pos.temp <- pos.temp / avg.pos
      norm.neg.temp <- neg.temp / avg.neg * -1

      norm.deeps <- c(norm.pos.temp,norm.neg.temp)
      if (is.finite(avg.pos)){
        DGSEA_Background.norm[,DGSEA_Background[i]] <- norm.deeps
      }
      real_dES <- DGSEA_Real[1,DGSEA_Background[i]]
      if (real_dES >= 0){
        Real.DGSEA_Background.NES[1,DGSEA_Background[i]] <- real_dES / avg.pos
      } else if (real_dES < 0){
        Real.DGSEA_Background.NES[1,DGSEA_Background[i]] <- real_dES / avg.neg * -1
      }
    }

    #Remove NA
    DGSEA_Background.norm <- DGSEA_Background.norm[,!is.na(colSums(DGSEA_Background.norm))]
    Real.DGSEA_Background.NES <- Real.DGSEA_Background.NES[,!is.na(colSums(Real.DGSEA_Background.NES))]
    #Now we have a distribution of DGSEA NES across all Gene.Sets.A.and.B and a distribution of Real DGSEA NES
    #We can use this for FDR calcs.

    DGSEA.Results <- matrix(data = NA, nrow = 1, ncol = 16)
    colnames(DGSEA.Results) <- c("Sample","Gene Set A", "Gene Set B", "Gene.Sets.Compared",
                                 "ES A","NES A","p-value A", "FDR A",
                                 "ES B","NES B","p-value B", "FDR B",
                                 "ES AB", "NES_AB","p_value_AB", "FDR")
    DGSEA.Results <- as.data.frame(DGSEA.Results)
    DGSEA.Results$Sample <- Samples[u]
    DGSEA.Results$`Gene Set A` <- Gene.Sets.A.and.B[1]
    DGSEA.Results$`Gene Set B` <- Gene.Sets.A.and.B[2]
    DGSEA.Results$Gene.Sets.Compared <- paste(Gene.Sets.A.and.B[1]," - ", Gene.Sets.A.and.B[2], sep = "")

    DGSEA.Results$`ES A` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[1],]$KS
    DGSEA.Results$`NES A` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[1],]$KS_Normalized
    DGSEA.Results$`p-value A` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[1],]$p_value
    DGSEA.Results$`FDR A` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[1],]$FDR_q_value

    DGSEA.Results$`ES B` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[2],]$KS
    DGSEA.Results$`NES B` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[2],]$KS_Normalized
    DGSEA.Results$`p-value B` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[2],]$p_value
    DGSEA.Results$`FDR B` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[2],]$FDR_q_value

    DGSEA.Results$`ES AB` <- GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[1],]$KS -
      GSEA.Results[GSEA.Results$Gene.Set == Gene.Sets.A.and.B[2],]$KS

    DGSEA.perms <- KSRandomArray$DGSEA

    DGSEA.perms.pos <- DGSEA.perms[which(DGSEA.perms > 0)]
    DGSEA.perms.neg <- DGSEA.perms[which(DGSEA.perms <= 0)]

    DGSEA.Norm.perms.pos <- DGSEA.perms.pos / mean(DGSEA.perms.pos)
    DGSEA.Norm.perms.neg <- DGSEA.perms.neg / mean(DGSEA.perms.neg) * -1

    percent.pos.DGSEA <- sum(DGSEA_Real > 0) / length(DGSEA_Real)
    percent.neg.DGSEA <- sum(DGSEA_Real < 0) / length(DGSEA_Real)

    if (DGSEA.Results$`ES AB` > 0){
      DGSEA.Results$NES_AB <- DGSEA.Results$`ES AB` / mean(DGSEA.perms.pos)
      DGSEA.Results$p_value_AB <- sum(DGSEA.perms.pos > DGSEA.Results$`ES AB`) / length(DGSEA.perms.pos)
      DGSEA.Results$FDR <- ifelse(sum(DGSEA.Norm.perms.pos > DGSEA.Results$NES_AB) / (length(DGSEA.Norm.perms.pos) * percent.pos.DGSEA) < 1,
                                  sum(DGSEA.Norm.perms.pos > DGSEA.Results$NES_AB) / (length(DGSEA.Norm.perms.pos) * percent.pos.DGSEA), 1)
    } else if(DGSEA.Results$`ES AB` < 0){
      DGSEA.Results$NES_AB <- DGSEA.Results$`ES AB` / mean(DGSEA.perms.neg) * -1
      DGSEA.Results$p_value_AB <- sum(DGSEA.perms.neg < DGSEA.Results$`ES AB`) / length(DGSEA.perms.neg)
      DGSEA.Results$FDR <- ifelse(sum(DGSEA.Norm.perms.neg < DGSEA.Results$NES_AB) / (length(DGSEA.Norm.perms.neg) * percent.neg.DGSEA) < 1,
                                  sum(DGSEA.Norm.perms.neg < DGSEA.Results$NES_AB) / (length(DGSEA.Norm.perms.neg) * percent.neg.DGSEA), 1)
    }

    DGSEA.Results.All.Samples <- rbind(DGSEA.Results.All.Samples,DGSEA.Results)
    GSEA.Results.All.Samples <- rbind(GSEA.Results.All.Samples,GSEA.Results)
    Mountain.Plot.Info.All.Samples <- c(Mountain.Plot.Info.All.Samples,Mountain.Plot.Info)
    rank_metric.All.Samples <- c(rank_metric.All.Samples, rank_metric)

  }

  snow::stopCluster(cl)
  rm(cl)

  return(list(DGSEA.Results = DGSEA.Results.All.Samples, GSEA.Results = GSEA.Results.All.Samples,
              Mountain.Plot.Info = Mountain.Plot.Info.All.Samples,
              ranking.metric = rank_metric.All.Samples))

}

