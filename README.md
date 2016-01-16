# fastPerm: Fast approximation of small p-values in permutation tests by partitioning the permutation space

This is an R package for quickly approximating small permutation p-values (e.g. p < 10^-6) for two-sample functions of the mean. It does this by partitioning the permutation space in such a way that there is a trend in the p-values across the partitions. It then calculates p-values in partitions that are cheap to evaluate, predicts p-values in partitions that are expensive to evaluate, and takes a weighted sum to get an overall p-value. A paper describing this method in more detail is currently under review.

fastPerm currently supports ratios and differences of the means. Functions of the median are experimental in this version of fastPerm, and may not be reliable.

## Installation

```{r}
# install.packages("devtools")
# library("devtools")
install_github("bdsegal/fastPerm")
```

## Examples

```{r}
library(fastPerm)

x <- rexp(100, 4)
y <- rexp(100, 2)
mStopRatioMean(x, y)
fastPerm(x, y)

x <- rnorm(110, 0, 1)
y <- rnorm(110, 1, 1)
mStopDiffMean(x, y)
fastPerm(x, y, testStat = diffMean)
```
