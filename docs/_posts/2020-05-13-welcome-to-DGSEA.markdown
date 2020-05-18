---
layout: post
title:  "Installing DGSEA"
date:   2020-05-13 14:48:58 -0700
categories: DGSEA
---

DGSEA can be found on [GitHub]. We are in the process of putting DGSEA on Bioconductor, but for now it can be directly installed using devtools in R:

{% highlight R %}
if (!require(devtools)) {
  install.packages('devtools')
}    
devtools::install_github('JamesJoly/DGSEA')
{% endhighlight %}

Check out the [Running DGSEA][https://jamesjoly.github.io/DGSEA/dgsea/2020/05/13/Running-DGSEA.html] page to see how to format your data to run your first analysis. Please file all bugs/feature requests at [DGSEA’s GitHub repo][DGSEA-gh].

[jekyll-docs]: https://jekyllrb.com/docs/home
[DGSEA-gh]:   https://github.com/JamesJoly/DGSEA
[GitHub]: https://github.com/JamesJoly/DGSEA
