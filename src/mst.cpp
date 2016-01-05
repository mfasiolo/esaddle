
#include "esaddle.h"

/*
 * Fast computation of Euclidean Minimum Spanning Tree using MLPACK library
 *
 * INPUT
 * X_: d by n matrix, where each column represents a different point in d dimensional space
*/

SEXP mst(SEXP X_)
{
  using namespace Rcpp;
    
  try{
    
    const arma::mat X = as<arma::mat>(X_);
    
    mlpack::emst::DualTreeBoruvka<> dtb(X);
    
    arma::mat out;
    
    dtb.ComputeMST(out);
    
    return wrap(out);
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return Rcpp::wrap(NA_REAL);
  
}