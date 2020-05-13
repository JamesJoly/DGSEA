---
layout: post
title:  "Installing DGSEA"
date:   2020-05-13 14:48:58 -0700
categories: DGSEA
---

DGSEA can be found on GitHub at https://github.com/JamesJoly/DGSEA. We are in the process of putting DGSEA on Bioconductor, but for now it can be directly installed using devtools in R:

{% highlight ruby %}
if (!require(devtools)) {
  install.packages('devtools')
}    
devtools::install_github('JamesJoly/DGSEA')
{% endhighlight %}

Check out the Running DGSEA page to run your first analysis. File all bugs/feature requests at [DGSEA’s GitHub repo][DGSEA-gh].

[jekyll-docs]: https://jekyllrb.com/docs/home
[DGSEA-gh]:   https://github.com/JamesJoly/DGSEA

