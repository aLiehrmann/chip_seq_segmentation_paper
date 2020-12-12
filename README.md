# Increased peak detection accuracy in over-dispersed ChIP-seq data with supervised segmentation models

**Motivation:** Histone modification constitutes a basic mechanism for the genetic regulation of gene expression. In early 2000s, a powerful technique has emerged that couples chromatin immunoprecipitation with high-throughput sequencing (ChIP-seq). This technique provides a direct survey of the DNA regions associated to these modifications. In order to realize the full potential of this technique, increasingly sophisticated statistical algorithms have been developed or adapted to analyze the massive amount of data it generates. Many of these algorithms were built around natural assumptions such as the Poisson one to model the noise in the count data. In this work we start from these natural assumptions and show that it is possible to improve upon them.

**Results:** The results of our comparisons on seven reference datasets of histone modifications (H3K36me3 \& H3K4me3) suggest that natural assumptions are not always realistic under application conditions. We show that the unconstrained multiple changepoint detection model, with alternative noise assumptions and a suitable setup, reduces the over-dispersion exhibited by count data and turns out to detect peaks more accurately than algorithms which rely on these natural assumptions.

**Keywords:** ChIP-seq; histone modifications; over-dispersion; peak calling; multiple changepoint detection; likelihood inference; supervised learning

*Data are available [here](https://www.dropbox.com/s/ipwe1ecb0mg12u8/data.tar.gz?dl=0).*