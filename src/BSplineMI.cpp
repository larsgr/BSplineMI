#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;


double SplineBlend( int curBin, int splineOrder, NumericVector knots, double v, int nBins){
  double value = 0.0;
  double d1,d2;

  if (splineOrder == 1) {

    if (((knots[curBin] <= v) && (v < knots[curBin+1])) ||
        ( abs(v - knots[curBin+1]) < 1e-10 && (curBin+1 == nBins) ) ){
      value = 1.0;
    } else {
      value = 0.0;
    }

  } else {

    d1 = knots[curBin+splineOrder-1] - knots[curBin];
    d2 = knots[curBin+splineOrder] - knots[curBin+1];

    if ( (d1 == 0) && (d2 == 0) ){
      value = 0;
    } else if (d1 == 0) {
      value = (knots[curBin+splineOrder] - v) / d2 * SplineBlend(curBin+1,splineOrder-1,knots,v,nBins);
    } else if (d2 == 0) {
      value = (v - knots[curBin]) / d1 * SplineBlend(curBin,splineOrder-1,knots,v,nBins);
    } else {
      value = (v - knots[curBin]) / d1 * SplineBlend(curBin,splineOrder-1,knots,v,nBins) +
        (knots[curBin+splineOrder] - v) / d2 * SplineBlend(curBin+1,splineOrder-1,knots,v,nBins);
    }
  }
  if (value < 0)  {
    //  this does not happen...
    value = 0; // rounding sometimes makes this < 0, e.g. -0.000000001
  }
  return(value);
}


// [[Rcpp::export]]
NumericVector SplineBlendAll( NumericVector z, NumericVector knots, int splineOrder, int nBins){
  int len = z.size();
  NumericVector retVal(nBins*len);

  for( int i = 0; i < len; i++){
    for(int curBin = 0; curBin < nBins; curBin++ ){
      retVal[curBin+i*nBins] = SplineBlend( curBin, splineOrder, knots, z[i], nBins);
    }
  }

  return(retVal);
}


// [[Rcpp::export]]
NumericVector hist2d_C(const NumericVector weights,
                       const int i1,
                       const int i2){
  // weights is a 3d array [bins x samples x genes ]
  NumericVector dim = weights.attr("dim");
  int nBins = dim[0];
  int nSamples = dim[1];

  int w1 = nBins*nSamples*i1; // add nBins*nSamples*i1 to get weights for gene i1
  int w2 = nBins*nSamples*i2; // add nBins*nSamples*i2 to get weights for gene i2

  int i,j,k;

  NumericMatrix retVal = NumericMatrix( nBins, nBins);

  for(i = 0; i < nSamples; i++){
    for(j = 0; j < nBins; j++){
      for(k = 0; k < nBins; k++){
        retVal(j,k) += weights[j + i*nBins + w1] * weights[k + i*nBins + w2];
      }
    }
  }

  // returns joint probability density matrix [bins x bins]:
  return( retVal / nSamples );
}



// [[Rcpp::export]]
double entropy2d_C(const NumericVector weights,
                   const int i1,
                   const int i2){
  NumericVector pMat = hist2d_C(weights, i1, i2);

  // calculate entropy
  // -sum(p * log2(p), na.rm = T)
  double H = 0.0; //returned entropy
  double p;
  int psize = pMat.size();
  for(int i = 0; i < psize; i++){
    p = pMat[i];
    // need to skip the 0's as it would generate NaN
    if( p > 0 ){
      H -= p*log2(p);
    }
  }
  return(H);
}

// [[Rcpp::export]]
double entropyHist2d_C(const NumericVector& weights,
                       const int i1,
                       const int i2,
                       const int nBins,
                       const int nSamples){

  // allocate joint probability matrix
  int psize = nBins*nBins;
  double *pMat = (double *)calloc(psize, sizeof(double));

  // weights is a 3d array [bins x samples x genes ]
  int w1 = nBins*nSamples*i1; // add nBins*nSamples*i1 to get weights for gene i1
  int w2 = nBins*nSamples*i2; // add nBins*nSamples*i2 to get weights for gene i2

  int i,j,k;

  for(i = 0; i < nSamples; i++){
    for(j = 0; j < nBins; j++){
      for(k = 0; k < nBins; k++){
        pMat[j + nBins*k] += weights[j + i*nBins + w1] * weights[k + i*nBins + w2];
      }
    }
  }

  // calculate entropy
  // -sum(p * log2(p), na.rm = T)
  double H = 0.0; //returned entropy

  for( i = 0; i < psize; i++){
    // divide by number of samples to get probability
    double p = pMat[i]/nSamples;
    // need to skip the 0's as it would generate NaN
    if( p > 0 ){
      H -= p*log2(p);
    }
  }

  free(pMat);

  return(H);
}

// [[Rcpp::export]]
NumericMatrix calcMIfromWeights(const NumericVector entropy,
                                const NumericVector weights,
                                const int threads){

  // get number of genes
  int nGenes = entropy.size();

  // get number of bins and samples
  NumericVector dim = weights.attr("dim");
  int nBins = dim[0];
  int nSamples = dim[1];

  // allocate MI matrix
  NumericMatrix mi = NumericMatrix(nGenes, nGenes);

#ifdef _OPENMP
  if ( threads > 1 )
    omp_set_num_threads( threads );
#else
  if( threads > 1 )
    Rf_warning("Parallel execution with openMP not supported. Running in single thread. Install openMP and recompile package for parallel execution.");
#endif

  #pragma omp parallel for shared(mi) schedule(dynamic)
  for( int i = 1; i < nGenes; i++){
    for( int j = 0; j < i; j++){
      // double miSingle = testdummy(weights, nBins, nSamples);
      double miSingle = entropy[j] + entropy[i] -
                        entropyHist2d_C(weights,i,j,nBins,nSamples);

      // In case the expression for one gene has no entropy the mi should be 0,
      // but precision error will produce not quite 0. Fix:
      if( miSingle < 0.000001 )
        miSingle = 0;

      mi[i + j*nGenes] = miSingle;
      mi[j + i*nGenes] = miSingle;
    }
  }

  return(mi);
}
