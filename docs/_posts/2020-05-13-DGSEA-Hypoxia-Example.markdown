---
layout: post
title:  "How to use DGSEA with your expression data"
date:   2020-05-13 14:48:58 -0700
categories: DGSEA
---

Now that we've installed and validated DGSEA, let's test the algorithm with real data and generate enrichment plots. 

In our manuscript, we used cellular response to hypoxia as a positive control to detect a shift towards glycolysis and away from oxidative phosphorylation. We used RNA sequencing data from 31 breast cancer cell lines subjected to 1% and 20% oxygen from:

Ye,I.C. et al. (2018) Molecular Portrait of Hypoxia in Breast Cancer: A Prognostic Signature and Novel HIF-Regulated Genes. Mol. Cancer Res. MCR, 16, 1889–1901.

I have truncated this data to only contain two of the 31 cell lines (MCF10A and MCF12A). 

{% highlight R %}
# Import data
library(RCurl)
data_in <- getURL("https://raw.githubusercontent.com/JamesJoly/DGSEA/master/MCF10A_MCF12A_Log2Ratio.csv")
data_in <- read.csv(text = data_in)
colnames(data_in)[1] <- "Gene"
{% endhighlight %}

{% highlight R %}
# Reformat to the GSEA PreRanked format
MCF10A <- data_in[,1:2]
MCF12A <- data_in[,-2]
{% endhighlight %}

{% highlight R %}
# Import gene sets
library(GSA)
kegg.pathways <- GSA.read.gmt("https://raw.githubusercontent.com/JamesJoly/DGSEA/master/KEGG_metabolic_pathways.gmt")
{% endhighlight %}

{% highlight R %}
# Perform targeted analysis examining Core Glycolysis - OxPhos
library(DGSEA)
set.names <- kegg.pathways$geneset.names

targeted.MCF10A <- dgsea_targeted(MCF10A, kegg.pathways,
                                  set.names[1], set.names[2])

targeted.MCF12A <- dgsea_targeted(MCF12A, kegg.pathways,
                                  set.names[1], set.names[2])
{% endhighlight %}

{% highlight R %}
# Generate Mountain Plots
MCF10A.mtn.plot <- make_mountain_plots(targeted.MCF10A, set.names[1], set.names[2])
MCF10A.mtn.plot
{% endhighlight %}

{:refdef: style="text-align: center;"}
![MCF10A_mtn_plot](https://raw.githubusercontent.com/JamesJoly/DGSEA/master/docs/assets/images/MCF10A_hypoxia_png.png){:height="70%" width="70%"}
{: refdef}


{% highlight R %}
# Generate Mountain Plots
MCF12A.mtn.plot <- make_mountain_plots(targeted.MCF12A, set.names[1], set.names[2])
MCF12A.mtn.plot
{% endhighlight %}

{:refdef: style="text-align: center;"}
![MCF12A_mtn_plot](https://raw.githubusercontent.com/JamesJoly/DGSEA/master/docs/assets/images/MCF12A_hypoxia_png.png){:height="70%" width="70%"}
{: refdef}

If you find DGSEA helpful, please cite us at: 

Joly, J. H., Lowry, W. E., and Graham, N. A. (2020) Differential Gene Set Enrichment Analysis: A statistical approach to quantify the relative enrichment of two gene sets. Bioinformatics. [10.1093/bioinformatics/btaa658][bioinformatics_link]


[DGSEA-gh]:   https://github.com/JamesJoly/DGSEA
[bioinformatics_link]:   https://doi.org/10.1093/bioinformatics/btaa658

