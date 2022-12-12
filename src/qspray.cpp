#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <Rcpp.h>
#include "gmp.h"

typedef std::vector<signed int> powers;

class PowersHasher {
 public:
  size_t operator()(const powers& exponents) const {
    // thanks to Steffan Hooper for advice
    std::size_t seed = 0;
    for(auto& i : exponents) {
      seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

typedef std::unordered_map<powers, CGAL::Gmpq, PowersHasher> qspray;

std::string q2str(CGAL::Gmpq r) {
  CGAL::Gmpz numer = r.numerator();
  CGAL::Gmpz denom = r.denominator();
  size_t n = mpz_sizeinbase(numer.mpz(), 10) + 2;
  size_t d = mpz_sizeinbase(denom.mpz(), 10) + 2;
  char* cnumer = new char[n];
  char* cdenom = new char[d];
  cnumer = mpz_get_str(cnumer, 10, numer.mpz());
  cdenom = mpz_get_str(cdenom, 10, denom.mpz());
  std::string snumer = cnumer;
  std::string sdenom = cdenom;
  delete[] cnumer;
  delete[] cdenom;
  return snumer + "/" + sdenom;
}

powers simplifyPowers(powers pows) {
  int n = pows.size();
  if(n == 0) {
    return pows;
  }
  n--;
  powers::iterator it = pows.end();
  bool zero = pows[n--] == 0;
  while(zero && n >= 0) {
    it--;
    zero = pows[n--] == 0;
  }
  pows.erase(it, pows.end());
  return pows;
}

qspray prepare(const Rcpp::List Powers, const Rcpp::StringVector coeffs) {
  qspray S;
  // powers pows;
  powers spows;
  signed int i, j;
  qspray::iterator it;

  for(i = 0; i < Powers.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers(i);
    CGAL::Gmpq coeff(Rcpp::as<std::string>(coeffs(i)));
    if(coeff != 0) {
      // pows.clear();
      powers pows(Exponents.begin(), Exponents.end());
      // for(j = 0; j < Exponents.size(); j++) {
      //   pows.push_back(Exponents(j));
      // }
      spows = simplifyPowers(pows);
      S[spows] += coeff;
    }
  }  // i loop closes

  // Now remove zero entries:
  it = S.begin();
  while(it != S.end()) {
    if(it->second == 0) {
      it = S.erase(it);  //  in C++11, erase() returns *next* iterator
    } else {
      ++it;  // else just increment the iterator
    }
  }
  return (S);
}

Rcpp::List makeindex(
    const qspray S) {  // takes a spray, returns the matrix of indices
  const unsigned int ncol = S.begin()->first.size();
  Rcpp::List out(S.size());  // index
  powers pows;
  unsigned int row = 0, col = 0;

  for(auto it = S.begin(); it != S.end(); ++it) {
    pows = it->first;
    Rcpp::IntegerVector Exponents(pows.size());
    col = 0;
    for(auto ci = pows.begin(); ci != pows.end(); ++ci) {
      Exponents(col++) = *ci;
    }
    out(row++) = Exponents;
  }
  return (out);
}

Rcpp::StringVector makevalue(const qspray S) {  // takes a spray, returns data
  Rcpp::StringVector out(S.size());             // data
  unsigned int i = 0;
  qspray::const_iterator it;  // it iterates through a sparse array

  for(it = S.begin(); it != S.end(); ++it) {
    out(i++) = q2str(it->second);  // initialize-and-fill is more efficient than
                                   // out.push_back(it->second)
  }
  return (out);
}

Rcpp::List retval(const qspray& S) {  // used to return a list to R

  // In this function, returning a zero-row matrix results in a
  // segfault ('memory not mapped').  So we check for 'S' being zero
  // size and, if so, return a special Nil value.  This corresponds to
  // an empty spray object.

  if(S.size() == 0) {
    return Rcpp::List::create(Rcpp::Named("index") = R_NilValue,
                              Rcpp::Named("value") = R_NilValue);
  } else {
    return Rcpp::List::create(Rcpp::Named("index") = makeindex(S),
                              Rcpp::Named("value") = makevalue(S));
  }
}

// [[Rcpp::export]]
Rcpp::List qspray_maker(const Rcpp::List& Powers,
                        const Rcpp::StringVector& coeffs) {
  return retval(prepare(Powers, coeffs));
}

// [[Rcpp::export]]
Rcpp::List qspray_add(const Rcpp::List& Powers1,
                      const Rcpp::StringVector& coeffs1,
                      const Rcpp::List& Powers2,
                      const Rcpp::StringVector& coeffs2) {
  qspray S1, S2;              // basic idea is S1 = S1+S2;
  qspray::const_iterator it;  // it iterates through a sparse array
  powers pows;

  S1 = prepare(Powers1, coeffs1);
  S2 = prepare(Powers2, coeffs2);

  // the "meat" of this function, namely S1=S1+S2 (S1 += S2):
  for(it = S2.begin(); it != S2.end(); ++it) {
    pows = it->first;
    S1[pows] +=
        S2[pows];  // S1 to be returned.  NB: S1 may have increased in size.
    if(S1[pows] == 0) {
      S1.erase(pows);
    }
  }

  return retval(S1);
}

void growPowers(powers& pows, signed int m, signed int n) {
  for(signed int i = m; i < n; i++) {
    pows.push_back(0);
  }
}

void harmonize(powers& pows1, powers& pows2) {
  signed int n1 = pows1.size();
  signed int n2 = pows2.size();
  if(n1 < n2) {
    growPowers(pows1, n1, n2);
  } else {
    growPowers(pows2, n2, n1);
  }
}

qspray prod(const qspray S1, const qspray S2) {
  qspray Sout;
  qspray::const_iterator it1, it2;
  powers powssum;
  unsigned int i;

  // the "meat" of this function:  Sout=S1*S2
  for(it1 = S1.begin(); it1 != S1.end(); ++it1) {
    const CGAL::Gmpq r1 = it1->second;
    if(r1 != 0) {
      powers pows1 = it1->first;
      for(it2 = S2.begin(); it2 != S2.end(); ++it2) {
        const CGAL::Gmpq r2 = it2->second;
        if(r2 != 0) {
          powers pows2 = it2->first;
          harmonize(pows1, pows2);
          powssum.clear();
          for(i = 0; i < pows1.size(); i++) {
            powssum.push_back(pows1[i] + pows2[i]);  // meat 1: powers add
          }
          Sout[powssum] += r1 * r2;
        }  // meat 2: coefficients multiply
      }
    }
  }
  //    for(spray::iterator it=S3.begin(); it != S3.end(); ++it){
  //    if(it->second ==0){S3.erase(it);}
  // }

  return Sout;
}

// [[Rcpp::export]]
Rcpp::List qspray_mult(const Rcpp::List& Powers1,
                       const Rcpp::StringVector& coeffs1,
                       const Rcpp::List& Powers2,
                       const Rcpp::StringVector& coeffs2) {
  return retval(prod(prepare(Powers1, coeffs1), prepare(Powers2, coeffs2)));
}

// [[Rcpp::export]]
void test() {
  powers pows = {1};
  growPowers(pows, pows.size(), 4);
  Rcpp::Rcout << pows.size();
  powers spows = simplifyPowers(pows);
  Rcpp::Rcout << spows.size();
}