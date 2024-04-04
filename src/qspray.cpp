/* Based on original code by Robin Hankin */

#include "qspray.h"


// -------------------------------------------------------------------------- //
qspray prepare(const Rcpp::List& Powers, const Rcpp::StringVector& coeffs) {
  qspray S;

  for(signed int i = 0; i < Powers.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers(i);
    gmpq coeff(Rcpp::as<std::string>(coeffs(i)));
    if(coeff != 0) {
      powers pows(Exponents.begin(), Exponents.end());
      simplifyPowers(pows);
      S[pows] += coeff;
    }
  }
  // Now remove zero entries:
  qspray::iterator it = S.begin();
  while(it != S.end()) {
    if(it->second == 0) {
      it = S.erase(it);  //  in C++11, erase() returns *next* iterator
    } else {
      ++it;  // else just increment the iterator
    }
  }
  
  return S;
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_maker(const Rcpp::List& Powers,
                        const Rcpp::StringVector& coeffs) {
  return retval(prepare(Powers, coeffs));
}


// -------------------------------------------------------------------------- //
qspray makeQspray(const Rcpp::List& Powers, const Rcpp::StringVector& coeffs) {
  qspray S;
  for(int i = 0; i < Powers.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers(i);
    gmpq coeff(Rcpp::as<std::string>(coeffs(i)));
    powers pows(Exponents.begin(), Exponents.end());
    S[pows] = coeff;
  }
  return S;
}


// -------------------------------------------------------------------------- //
Rcpp::List retval(const qspray& S) {  // used to return a list to R

  // In this function, returning a zero-row matrix results in a
  // segfault ('memory not mapped').  So we check for 'S' being zero
  // size and, if so, return a special Nil value.  This corresponds to
  // an empty spray object.

  if(S.size() == 0) {
    return Rcpp::List::create(Rcpp::Named("powers") = R_NilValue,
                              Rcpp::Named("coeffs") = R_NilValue);
  } else {
    Rcpp::List Powers(S.size());
    powers pows;
    unsigned int row = 0, col = 0;
    Rcpp::StringVector Coeffs(S.size());
    unsigned int i = 0;
    for(auto it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      Rcpp::IntegerVector Exponents(pows.size());
      col = 0;
      for(auto ci = pows.begin(); ci != pows.end(); ++ci) {
        Exponents(col++) = *ci;
      }
      Powers(row++) = Exponents;
      Coeffs(i++) = q2str(it->second);
    }
    return Rcpp::List::create(Rcpp::Named("powers") = Powers,
                              Rcpp::Named("coeffs") = Coeffs);
  }
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_deriv(
    const Rcpp::List& Powers, const Rcpp::StringVector& coeffs,   
    const Rcpp::IntegerVector& n
){
  Qspray<gmpq> Q(makeQspray(Powers, coeffs));
  std::vector<unsigned int> orders(n.begin(), n.end());
  Qspray<gmpq> Qprime = Q.deriv(orders);
  return retval(Qprime.get());
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_add(const Rcpp::List& Powers1,
                      const Rcpp::StringVector& coeffs1,
                      const Rcpp::List& Powers2,
                      const Rcpp::StringVector& coeffs2) {
  Qspray<gmpq> Q1(makeQspray(Powers1, coeffs1));
  Qspray<gmpq> Q2(makeQspray(Powers2, coeffs2));  
  Qspray<gmpq> Q = Q1 + Q2;
  return retval(Q.get());
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_subtract(const Rcpp::List& Powers1,
                           const Rcpp::StringVector& coeffs1,
                           const Rcpp::List& Powers2,
                           const Rcpp::StringVector& coeffs2) {
  Qspray<gmpq> Q1(makeQspray(Powers1, coeffs1));
  Qspray<gmpq> Q2(makeQspray(Powers2, coeffs2));  
  Qspray<gmpq> Q = Q1 - Q2;
  return retval(Q.get());
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_mult(const Rcpp::List& Powers1,
                       const Rcpp::StringVector& coeffs1,
                       const Rcpp::List& Powers2,
                       const Rcpp::StringVector& coeffs2) {
  Qspray<gmpq> Q1(makeQspray(Powers1, coeffs1));
  Qspray<gmpq> Q2(makeQspray(Powers2, coeffs2));  
  Qspray<gmpq> Q = Q1 * Q2;
  return retval(Q.get());
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
bool qspray_equality(const Rcpp::List& Powers1,
                     const Rcpp::StringVector& coeffs1,
                     const Rcpp::List& Powers2,
                     const Rcpp::StringVector& coeffs2) {
  Qspray<gmpq> Q1(makeQspray(Powers1, coeffs1));
  Qspray<gmpq> Q2(makeQspray(Powers2, coeffs2));  
  return Q1 == Q2;
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_power(const Rcpp::List& Powers,
                        const Rcpp::StringVector& coeffs,
                        unsigned int n) {
   Qspray<gmpq> Q(makeQspray(Powers, coeffs));
   return retval(Q.power(n).get());
}
