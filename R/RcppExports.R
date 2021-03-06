# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

SplineBlendAll <- function(z, knots, splineOrder, nBins) {
    .Call('_BSplineMI_SplineBlendAll', PACKAGE = 'BSplineMI', z, knots, splineOrder, nBins)
}

hist2d_C <- function(weights, i1, i2) {
    .Call('_BSplineMI_hist2d_C', PACKAGE = 'BSplineMI', weights, i1, i2)
}

entropy2d_C <- function(weights, i1, i2) {
    .Call('_BSplineMI_entropy2d_C', PACKAGE = 'BSplineMI', weights, i1, i2)
}

entropyHist2d_C <- function(weights, i1, i2, nBins, nSamples) {
    .Call('_BSplineMI_entropyHist2d_C', PACKAGE = 'BSplineMI', weights, i1, i2, nBins, nSamples)
}

calcMIfromWeights <- function(entropy, weights, threads) {
    .Call('_BSplineMI_calcMIfromWeights', PACKAGE = 'BSplineMI', entropy, weights, threads)
}

