# hhg-test
Heller-Heller-Gorfine multivariate test of association

## Getting Started

`inversions.c` must be compiled from C source code into a [MEX file](https://www.mathworks.com/help/matlab/ref/mex.html).

## Prerequisites

This code was tested in MATLAB R2017b.
A [MATLAB-compatible C compiler](https://www.mathworks.com/support/compilers.html) is required, as is the [MATLAB Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/index.html).

## Usage

`HHGPermutationTest` takes five input arguments:
1. `X`: The first input data matrix. Rows are assumed to represent samples, and columns are assumed to represent dimensions.
2. `Y`: The second input data matrix. `Y` must contain the same number of samples as `X`.
3. `nperm`: The number of permutations to perform. The default value is 100.
4. `maxN`: The maximum sample size for which the full distance matrices will be held in memory. If the sample size exceeds `maxN`, the distance matrices will be computed incrementally and stored on disk.

Once computation is completed, `HHGPermutationTest` returns up to three output arguments:
1. `p`: The p-value of the permutation test. If `p < 1/nperm`, a value of 0 is returned.
2. `t`: The value of the HHG test statistic.
3. `pstat`: The values of the HHG test statistic for each permutation.

## References

1. Heller, R., Heller, Y., & Gorfine, M. (2012). "A consistent multivariate test of association based on ranks of distances." _Biometrika_, 100(2), 503-510. [(Link)](https://academic.oup.com/biomet/article/100/2/503/202568)

## License

This project is licensed under the [MIT License](LICENSE.txt).
