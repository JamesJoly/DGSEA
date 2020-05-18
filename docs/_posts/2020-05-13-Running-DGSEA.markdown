---
layout: post
title:  "Running DGSEA"
date:   2020-05-13 14:48:58 -0700
categories: DGSEA
---

Once you've installed DGSEA into your R library, you can run either a targeted analysis or untargeted analysis. Targeted analyses are generally used when a hypothesis is known *a priori*, whereas untargeted analyses are for hypothesis generation.

Regardless, the inputs to the DGSEA functions are similar to that of the GSEA PreRanked function in the Broad Institute's GSEA Java applet. An example for how to format the data is shown below. Similar code was used to generate Figure 2 in our paper (Joly et al. 2020). 

{% highlight R %}
# Generate 16,0000 random gene names
gene.names <- as.character(c(1:16000))
gene.names <- paste("Gene",gene.names, sep = "")
{% endhighlight %}

{% highlight R %}
# Generate gene sets, here we choose gene set of sizes 10, 25, 50, 100
Gene.Set.A.100 <- gene.names[1:100]
Gene.Set.B.100 <- gene.names[101:200]

Gene.Set.A.10 <- gene.names[1:10]
Gene.Set.B.10 <- gene.names[101:110]

Gene.Set.A.50 <- gene.names[1:50]
Gene.Set.B.50 <- gene.names[101:150]

Gene.Set.A.25 <- gene.names[1:25]
Gene.Set.B.25 <- gene.names[101:125]
{% endhighlight %}

{% highlight R %}
# Assemble these to have the output that a .gmt file would have
names <- c("Gene.Set.A.100","Gene.Set.B.100","Gene.Set.A.10", "Gene.Set.B.10",
           "Gene.Set.A.50", "Gene.Set.B.50", "Gene.Set.A.25", "Gene.Set.B.25")
sets <- list(Gene.Set.A.100,Gene.Set.B.100,Gene.Set.A.10, Gene.Set.B.10,
             Gene.Set.A.50, Gene.Set.B.50, Gene.Set.A.25, Gene.Set.B.25)
Gene.Sets.gmt.list <- list(genesets = sets, geneset.names = names) #Same format as GSA::GSA.read.gmt()
{% endhighlight %}

{% highlight R %}
# Generate data with standard normal distribution
expression.data <- rnorm(length(gene.names), mean = 0, sd = 1)
{% endhighlight %}

{% highlight R %}
# Replace Gene Set A and Gene Set B expression values with +X and -X
# Here you can play with X to see how strong of an enrichment is necessary
X = 1
positions.Gene.Set.A.10 <- which(gene.names %in% Gene.Set.A.10)
positions.Gene.Set.B.10 <- which(gene.names %in% Gene.Set.B.10)
expression.data[positions.Gene.Set.A.10] <- rnorm(length(positions.Gene.Set.A.10), mean = X, sd = 1)
expression.data[positions.Gene.Set.B.10] <- rnorm(length(positions.Gene.Set.B.10), mean = -X, sd = 1)
{% endhighlight %}

{% highlight R %}
# Assemble GSEA Pre-ranked format
library(dplyr)
Pre.ranked.form <- cbind(gene.names, expression.data) %>% as.data.frame()
colnames(Pre.ranked.form) <-c("Gene", "Expression.data")
{% endhighlight %}

{% highlight R %}
#Run dgsea untargeted
library(DGSEA)

untargeted <- dgsea_untargeted(Pre.ranked.form, gmt.list = Gene.Sets.gmt.list)
DGSEA.results.untargeted <- untargeted$DGSEA.Results
# Gene Sets: Gene.Set.A.10 and Gene.Set.B.10 should have strong positive and negative enrichments
{% endhighlight %}

{% highlight R %}
# Run targeted dgsea

library(DGSEA)

geneset.names <- Gene.Sets.gmt.list$geneset.names
targeted <- dgsea_targeted(Pre.ranked.form, gmt.list = Gene.Sets.gmt.list,
                           Gene.Set.A.Name = geneset.names[1], Gene.Set.B.Name = geneset.names[2])
DGSEA.results.targeted <- targeted$DGSEA.Results
{% endhighlight %}

[jekyll-docs]: https://jekyllrb.com/docs/home
[DGSEA-gh]:   https://github.com/JamesJoly/DGSEA

