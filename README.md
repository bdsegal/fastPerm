# fastPerm

This is an R package for quickly approximating small p-values (e.g. p < 10<sup>-6</sup>) in two-sample permutation tests, in which the test statistic is the difference or ratio of the means. It does this by partitioning the permutations in such a way that there is a trend in the p-values across the partitions. It then calculates p-values in partitions that are cheap to evaluate, predicts p-values in partitions that are expensive to evaluate, and takes a weighted sum to get an overall p-value. A paper describing this method in detail is currently under review.

fastPerm currently supports ratios and differences of the means (`testStat = ratioMean` and `testStat = diffMean`, respectively).

This package also implements stochastic approximation Monte Carlo (Yu et al., 2011) in the `SAMC` function.

## Installation

```{r}
# install.packages("devtools")
library("devtools")
install_github("bdsegal/fastPerm")
```

## Examples

```{r}
library(fastPerm)

# Ratio of means
x <- rexp(100, 5)
y <- rexp(100, 2)
# Expected stopping partition, based on asymptotic approximation
mStopRatioMean(x, y)
# Fast permutation test, defaults to 'testStat = ratioMean'
fastPerm(x, y)

# Difference of means
x <- rnorm(110, 0, 1)
y <- rnorm(95, 1, 1)
# Expected stopping partition, based on asymptotic approximation
mStopDiffMean(x, y)
# Fast permutation test
fastPerm(x, y, testStat = diffMean)
```

## References

Kai Yu, Faming Liang, Julia Ciampa, and Nilanjan Chatterjee. Efficient p-value evaluation for resampling-based tests. Biostatistics, pages 1-11, 2011.

Brian Segal, Thomas Braun, Michael Elliott, and Hui Jiang. Fast approximation of small p-values in permutation tests by partitioning the permutations. Submitted. Working draft available online at [https://arxiv.org/abs/1605.03992](https://arxiv.org/abs/1605.03992).