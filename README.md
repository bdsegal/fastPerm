# fastPerm

This is an R package for quickly approximating small p-values (e.g. p < 10<sup>-6</sup>) for the difference and ratio of means in two-sample tests as described by Segal et al. (2018).

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

Segal, B. D., Braun, T., Elliott, M. R. and Jiang, H. (2018). Fast approximation of small p-values in permutation tests by partitioning the permutations. Biometrics, 74(1), 196-206. [doi:10.1111/biom.12731](http://dx.doi.org/10.1111/biom.12731).

Yu, K., Liang, F., Ciampa, J., and Chatterjee, N. (2011). Efficient p-value evaluation for
resampling-based tests. Biostatistics 12, 582–593.
