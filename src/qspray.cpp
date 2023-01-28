/* Based on original code by Robin Hankin */

#include <Rcpp.h>
#include <boost/multiprecision/gmp.hpp>
#include <complex.h>
typedef std::vector<signed int>                             powers;
typedef boost::multiprecision::mpq_rational                 gmpq;
typedef boost::multiprecision::mpz_int                      gmpi;
typedef std::complex<gmpq>                                  qcplx;

// -------------------------------------------------------------------------- //
qcplx qxmult(qcplx z1, qcplx z2) {
  gmpq r1 = z1.real();
  gmpq i1 = z1.imag();
  gmpq r2 = z2.real();
  gmpq i2 = z2.imag();
  qcplx result(r1*r2 - i1*i2, r1*i2+r2*i1);
  return result;
}


// -------------------------------------------------------------------------- //
qcplx qxpow(qcplx z, unsigned k) {
  qcplx result(gmpq("1"), gmpq("0"));
  while(k) {
    if(k & 1) {
      result = qxmult(result, z);
    }
    k >>= 1;
    z = qxmult(z, z);
  }
  return result;
}


// -------------------------------------------------------------------------- //
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

typedef std::unordered_map<powers, gmpq, PowersHasher> qspray;


// -------------------------------------------------------------------------- //
std::string q2str(gmpq r) {
  const gmpi numer = boost::multiprecision::numerator(r);
  const gmpi denom = boost::multiprecision::denominator(r);
  mpz_t p;
  mpz_init(p);
  mpz_set(p, numer.backend().data());
  mpz_t q;
  mpz_init(q);
  mpz_set(q, denom.backend().data());
  const size_t n = mpz_sizeinbase(p, 10) + 2;
  const size_t d = mpz_sizeinbase(q, 10) + 2;
  char* cnumer = new char[n];
  char* cdenom = new char[d];
  cnumer = mpz_get_str(cnumer, 10, p);
  cdenom = mpz_get_str(cdenom, 10, q);
  std::string snumer = cnumer;
  std::string sdenom = cdenom;
  delete[] cnumer;
  delete[] cdenom;
  mpz_clear(p);
  mpz_clear(q);
  return snumer + "/" + sdenom;
}


// -------------------------------------------------------------------------- //
void simplifyPowers(powers& pows) {
  int n = pows.size();
  if(n == 0) {
    return;
  }
  n--;
  powers::iterator it = pows.end();
  bool zero = pows[n] == 0;
  while(zero && n >= 0) {
    it--;
    n--;
    zero = pows[n] == 0;
  }
  if(n == -1) {
    pows = {};
  } else {
    pows.erase(it, pows.end());
  }
}


// -------------------------------------------------------------------------- //
qspray prepare(const Rcpp::List Powers, const Rcpp::StringVector coeffs) {
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
qspray makeQspray(const Rcpp::List Powers, const Rcpp::StringVector coeffs) {
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
// [[Rcpp::export]]
Rcpp::StringVector evalQxspray(const Rcpp::List Powers,
                               const Rcpp::StringVector coeffs,
                               const Rcpp::StringVector v_re,
                               const Rcpp::StringVector v_im) {
  qspray S = makeQspray(Powers, coeffs);

  std::vector<qcplx> v;
  int n = v_re.size();
  v.reserve(n);
  for(int i = 0; i < n; i++) {
    qcplx vi(gmpq(Rcpp::as<std::string>(v_re(i))),
             gmpq(Rcpp::as<std::string>(v_im(i))));
    v.emplace_back(vi);
  }

  powers pows;
  gmpq coef;

  qcplx result(gmpq("0"), gmpq("0"));

  for(auto it = S.begin(); it != S.end(); ++it) {
    pows = it->first;
    coef = it->second;
    qcplx term(coef, gmpq("0"));
    int i = 0;
    for(auto ci = pows.begin(); ci != pows.end(); ++ci) {
      term = qxmult(term, qxpow(v[i++], *ci));
    }
    result += term;
  }

  return Rcpp::StringVector::create(q2str(result.real()), q2str(result.imag()));
}


// -------------------------------------------------------------------------- //
Rcpp::List makepowers(const qspray S) {  // returns the list of powers
  Rcpp::List Powers(S.size());
  powers pows;
  unsigned int row = 0, col = 0;

  for(auto it = S.begin(); it != S.end(); ++it) {
    pows = it->first;
    Rcpp::IntegerVector Exponents(pows.size());
    col = 0;
    for(auto ci = pows.begin(); ci != pows.end(); ++ci) {
      Exponents(col++) = *ci;
    }
    Powers(row++) = Exponents;
  }
  return Powers;
}


// -------------------------------------------------------------------------- //
Rcpp::StringVector makecoeffs(const qspray S) {  // returns coefficients
  Rcpp::StringVector Coeffs(S.size());
  unsigned int i = 0;
  qspray::const_iterator it;
  for(it = S.begin(); it != S.end(); ++it) {
    Coeffs(i++) = q2str(it->second);
  }
  return Coeffs;
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
    return Rcpp::List::create(Rcpp::Named("powers") = makepowers(S),
                              Rcpp::Named("coeffs") = makecoeffs(S));
  }
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_maker(const Rcpp::List& Powers,
                        const Rcpp::StringVector& coeffs) {
  return retval(prepare(Powers, coeffs));
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_add(const Rcpp::List& Powers1,
                      const Rcpp::StringVector& coeffs1,
                      const Rcpp::List& Powers2,
                      const Rcpp::StringVector& coeffs2) {
  qspray::const_iterator it;
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);
  powers pows;
  
  for(it = S2.begin(); it != S2.end(); ++it) {
    pows = it->first;
    S1[pows] += S2[pows];
    if(S1[pows] == 0) {
      S1.erase(pows);
    }
  }

  return retval(S1);
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_subtract(const Rcpp::List& Powers1,
                           const Rcpp::StringVector& coeffs1,
                           const Rcpp::List& Powers2,
                           const Rcpp::StringVector& coeffs2) {
  qspray::const_iterator it;
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);
  powers pows;
  
  for(it = S2.begin(); it != S2.end(); ++it) {
    pows = it->first;
    S1[pows] -= S2[pows];
    if(S1[pows] == 0) {
      S1.erase(pows);
    }
  }

  return retval(S1);
}


// -------------------------------------------------------------------------- //
powers growPowers(powers pows, signed int m, signed int n) {
  powers gpows;
  gpows.reserve(n);
  for(signed int i = 0; i < m; i++) {
    gpows.emplace_back(pows[i]);
  }
  for(signed int i = m; i < n; i++) {
    gpows.emplace_back(0);
  }
  return gpows;
}


// -------------------------------------------------------------------------- //
qspray prod(const qspray S1, const qspray S2) {
  qspray Sout;
  qspray::const_iterator it1, it2;
  powers powssum;
  signed int i;

  for(it1 = S1.begin(); it1 != S1.end(); ++it1) {
    const gmpq r1 = it1->second;
    if(r1 != 0) {
      powers pows1 = it1->first;
      signed int n1 = pows1.size();
      for(it2 = S2.begin(); it2 != S2.end(); ++it2) {
        const gmpq r2 = it2->second;
        if(r2 != 0) {
          powers pows2 = it2->first;
          signed int n2 = pows2.size();
          powssum.clear();
          if(n1 < n2) {
            powers gpows = growPowers(pows1, n1, n2);
            powssum.reserve(n2);
            for(i = 0; i < n2; i++) {
              powssum.emplace_back(gpows[i] + pows2[i]);
            }
          } else if(n1 > n2) {
            powers gpows = growPowers(pows2, n2, n1);
            powssum.reserve(n1);
            for(i = 0; i < n1; i++) {
              powssum.emplace_back(pows1[i] + gpows[i]);
            }
          } else {
            powssum.reserve(n1);
            for(i = 0; i < n1; i++) {
              powssum.emplace_back(pows1[i] + pows2[i]);
            }
          }
          Sout[powssum] += r1 * r2;
        }
      }
    }
  }

  return Sout;
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_mult(const Rcpp::List& Powers1,
                       const Rcpp::StringVector& coeffs1,
                       const Rcpp::List& Powers2,
                       const Rcpp::StringVector& coeffs2) {
  return retval(
      prod(makeQspray(Powers1, coeffs1), makeQspray(Powers2, coeffs2)));
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
bool qspray_equality(const Rcpp::List& Powers1,
                     const Rcpp::StringVector& coeffs1,
                     const Rcpp::List& Powers2,
                     const Rcpp::StringVector& coeffs2) {
  powers pows;
  qspray::const_iterator it;
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);

  if(S1.size() != S2.size()) {
    return false;
  }

  for(it = S1.begin(); it != S1.end(); ++it) {
    pows = it->first;
    if(S1[pows] != S2[pows]) {
      return false;
    } else {
      S2.erase(pows);
    }
  }
  // at this point, S1[v] == S2[v] for every index 'v' of S1;  S1\subseteq S2.
  // We need to check that every element of S2 has been accounted for:
  if(S2.empty()) {
    return true;
  } else {
    return false;
  }
}


// -------------------------------------------------------------------------- //
qspray unit() {
  qspray out;
  powers pows(0);
  gmpq one(1);
  out[pows] = one;
  return out;
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_power(const Rcpp::List& Powers,
                        const Rcpp::StringVector& coeffs,
                        unsigned int n) {
  qspray out;
  if(n >= 1) {
    const qspray S = makeQspray(Powers, coeffs);
    if(n == 1) {
      out = S;
    } else {
      out = unit();
      for(; n > 0; n--) {
        out = prod(S, out);
      }
    }
  } else {
    out = unit();
  }

  return retval(out);
}
