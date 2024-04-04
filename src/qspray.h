#ifndef __HEADER__
#define __HEADER__

#include <Rcpp.h>
#include <boost/multiprecision/gmp.hpp>
#include <complex.h>
typedef std::vector<signed int>                             powers;
typedef boost::multiprecision::mpq_rational                 gmpq;
typedef boost::multiprecision::mpz_int                      gmpi;
typedef std::complex<gmpq>                                  qcplx;

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

std::string q2str(gmpq);
void simplifyPowers(powers&);
powers growPowers(powers, signed int, signed int);
qspray makeQspray(const Rcpp::List&, const Rcpp::StringVector&); 
Rcpp::List retval(const qspray&);


// -------------------------------------------------------------------------- //
template<typename T>
class Qspray {

  std::unordered_map<powers,T,PowersHasher> S;

public:
  // constructors -----
  Qspray()
    : S()
      {}

  Qspray(const std::unordered_map<powers,T,PowersHasher>& S_) 
    : S(S_) 
      {}

  Qspray(int k)
    : S(scalarQspray(T(k)))
      {}
  
  // methods -----
  std::unordered_map<powers,T,PowersHasher> get() {
    return S;
  } 

  bool empty() {
    return S.empty();
  }

  bool isNull() {
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    T zero;
    bool result = true;
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      if(it->second != zero) {
        result = false;
        break;
      }
    } 
    return result;
  }

  bool operator==(const Qspray<T>& Q) {
    typename std::unordered_map<powers,T,PowersHasher> SS(S);
    typename std::unordered_map<powers,T,PowersHasher> S2(Q.S);
    if(S.size() != S2.size()) {
      return false;
    }
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      if(SS[pows] != S2[pows]) {
        return false;
      } else {
        S2.erase(pows);
      }
    }
    // at this point, S2[k] == S[k] for every index 'k' of S;
    // it remains to check that every element of Q has been accounted for:
    if(S2.empty()) {
      return true;
    } else {
      return false;
    }
  }

  Qspray<T> operator-() {
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;  
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      S[pows] = -it->second;
    }
    return Qspray<T>(S);
  }

  Qspray<T> operator+=(const Qspray<T>& Q) {
    typename std::unordered_map<powers,T,PowersHasher> S2 = Q.S;
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    const T zero; // il s'agira de définir T() pour les ratio of polynomials :-)
    for(it = S2.begin(); it != S2.end(); ++it) {
      pows = it->first;
      S[pows] += it->second;
      if(S[pows] == zero) {
        S.erase(pows);
      }
    }
    return Qspray<T>(S);
  }

  Qspray<T> operator+(const Qspray<T>& Q2) {
    Qspray<T> Q(S);
    Q += Q2;
    return Q;
  }

  Qspray<T> operator-=(const Qspray<T>& Q) {
    typename std::unordered_map<powers,T,PowersHasher> S2 = Q.S;
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    const T zero; // il s'agira de définir T() pour les ratio of polynomials :-)
    for(it = S2.begin(); it != S2.end(); ++it) {
      pows = it->first;
      S[pows] -= it->second;
      if(S[pows] == zero) {
        S.erase(pows);
      }
    }
    return Qspray<T>(S);
  }

  Qspray<T> operator-(const Qspray<T>& Q2) {
    Qspray<T> Q(S);
    Q -= Q2;
    return Q;
  }

  Qspray<T> operator*=(const Qspray<T>& Q) {
    typename std::unordered_map<powers,T,PowersHasher> S2 = Q.S;
    typename std::unordered_map<powers,T,PowersHasher> Sout;
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it1, it2;
    const T zero;
    powers powssum;
    signed int i;
    for(it1 = S.begin(); it1 != S.end(); ++it1) {
      const T r1 = it1->second;
      if(r1 != zero) {
        powers pows1 = it1->first;
        signed int n1 = pows1.size();
        for(it2 = S2.begin(); it2 != S2.end(); ++it2) {
          const gmpq r2 = it2->second;
          if(r2 != zero) {
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
    S = Sout;
    return Qspray<T>(Sout);
  }

  Qspray<T> operator*(const Qspray<T>& Q2) {
    Qspray<T> Q(S);
    Q *= Q2;
    return Q;
  }

  Qspray<T> power(unsigned int n) {
    typename std::unordered_map<powers,T,PowersHasher> u;
    powers pows(0);
    T      one(1);
    u[pows] = one;
    Qspray<T> Result(u);
    Qspray<T> Q(S);
    while(n) {
      if(n & 1) {
        Result *= Q;
      }
      n >>= 1;
      Q *= Q;
    }
    return Result;
  }
  
  Qspray<T> deriv(std::vector<unsigned int> n) {
    typename std::unordered_map<powers,T,PowersHasher> Sprime;
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers v;
    signed int j, J, nj, expnt;
    signed int N = n.size();
    T zero;

    for(it = S.begin(); it != S.end(); ++it) {
      std::vector<signed int> exponents = it->first;
      J = exponents.size();
      if(J < N) {
        continue;
      }
      T coeff = S[it->first];
      for(j = 0; j < N; j++) {
        nj = n[j];
        while((nj > 0) && (coeff != zero)) { // while loop because it might not run at all
          coeff *= exponents[j]; // multiply coeff first, then decrement exponent 
          exponents[j]--;
          nj--;
        }
      }
      if(coeff != zero) {
        v.clear();
        for(j = 0; j < J; j++) {
          expnt = exponents[j];
          v.push_back(expnt);
        }
        simplifyPowers(v);
        Sprime[v] += coeff;  // increment because v is not row-unique any more
      }
    }  // loop closes
    
    return Qspray<T>(Sprime);
  }

};


// -------------------------------------------------------------------------- //
template <typename T>
Qspray<T> scalarQspray(T x) {
  typename std::unordered_map<powers,T,PowersHasher> singleton;
  powers pows(0);
  singleton[pows] = x;
  return Qspray<T>(singleton);
}

template <typename T>
int numberOfVariables(const Qspray<T>& Q) {
  std::unordered_map<powers,T,PowersHasher> S = Q.get();
  int d = 0;
  for(const auto& term : S) {
    int n = term.first.size();
    if(n > d) {
      d = n;
    }
  }
  return d;
}

Qspray<gmpq> gcdQsprays(const Qspray<gmpq>& Q1, const Qspray<gmpq>& Q2) {
  return scalarQspray<gmpq>(1);
}

int lexLeadingIndex(std::vector<powers> expnts) {
  const int n = expnts.size();
  int i = 0;
  while(i < n-1) {
    powers vi = expnts[i];
    for(int j = i + 1; j < n; j++) {
      powers vj = expnts[j];
      bool vjmax = std::lexicographical_compare(
        std::begin(vi), std::end(vi), std::begin(vj), std::end(vj)
      );
      if(vjmax) {
        i = j - 1;
        break;
      } else if(j == n-1) {
        return i;
      }
    }
    i++;
  }
  return i;
}

Rcpp::List leadingTerm(Qspray<gmpq>& Q, int d) {
  qspray S = Q.get();
  std::vector<powers> pows;
  std::vector<gmpq>   coeffs;
  pows.reserve(S.size());
  coeffs.reserve(S.size());
  for(const auto& term : S) {
    pows.emplace_back(term.first);
    coeffs.emplace_back(term.second);
  }
  int index = lexLeadingIndex(pows);
  powers leadingPows = pows[index];
  int npows = leadingPows.size();
  if(npows < d) {
    leadingPows = growPowers(leadingPows, npows, d);
  }
  std::string leadingCoeff = q2str(coeffs[index]);
  Rcpp::IntegerVector powsRcpp(leadingPows.begin(), leadingPows.end());
  return Rcpp::List::create(
    Rcpp::Named("powers") = powsRcpp,
    Rcpp::Named("coeff")  = leadingCoeff
  );
}

bool divides(Rcpp::List f, Rcpp::List g) {
  Rcpp::IntegerVector pows_f = f["powers"];
  Rcpp::IntegerVector pows_g = g["powers"];
  int n = pows_f.size();
  int i = 0;
  bool out = true;
  while(out && i < n) {
    out = out && (pows_f(i) <= pows_g(i));
    i++;
  }
  return out;
}

Qspray<gmpq> quotient(Rcpp::List f, Rcpp::List g) {
  Rcpp::IntegerVector pows_f = f["powers"];
  std::string coeff_f        = f["coeff"];
  Rcpp::IntegerVector pows_g = g["powers"];
  std::string coeff_g        = g["coeff"];
  gmpq qcoeff_f(coeff_f);
  gmpq qcoeff_g(coeff_g);
  Rcpp::IntegerVector powsRcpp = pows_f - pows_g;
  gmpq qcoeff = qcoeff_f / qcoeff_g;
  qspray S;
  powers pows(powsRcpp.begin(), powsRcpp.end());
  simplifyPowers(pows);
  S[pows] = qcoeff;
  return Qspray<gmpq>(S);
}

std::pair<Qspray<gmpq>,Qspray<gmpq>> qsprayDivision(
  Qspray<gmpq>& p, Qspray<gmpq>& g, int d
) {
  Rcpp::List LTg = leadingTerm(g, d);
  Qspray<gmpq> q;
  Qspray<gmpq> r;
  bool divoccured;
  while(!p.empty()) {
    divoccured = false;
    Rcpp::List LTp = leadingTerm(p, d);
    if(divides(LTg, LTp)) {
      Qspray<gmpq> Qtnt = quotient(LTp, LTg);
      p -= Qtnt * g;
      q += Qtnt;
      divoccured = true;
    } 
    if(!divoccured) {
      Rcpp::IntegerVector powsRcpp = LTp["powers"];
      std::string coeff            = LTp["coeff"];
      gmpq   coef(coeff);
      powers pows(powsRcpp.begin(), powsRcpp.end());
      simplifyPowers(pows);
      qspray LTpspray;
      LTpspray[pows] = coef;
      Qspray<gmpq> ltp(LTpspray);
      r += ltp;
      p -= ltp;
    }
  }
  return std::pair<Qspray<gmpq>,Qspray<gmpq>>(q, r);
}

Qspray<gmpq> QuotientQsprays(Qspray<gmpq>& A, Qspray<gmpq>& B, int d) {
  return qsprayDivision(A, B, d).first;
}

// ---------------------------------------------------------------------------//
template<typename T>
class RatioOfQsprays {

  Qspray<T> numerator;
  Qspray<T> denominator;
  int       dimension;

public:
  // constructors -----
  RatioOfQsprays()
    : numerator(scalarQspray<T>(T(0))), 
      denominator(scalarQspray<T>(T(1))),
      dimension(0)
      {}

  RatioOfQsprays(const Qspray<T>& numerator_, const Qspray<T>& denominator_) 
    : numerator(numerator_), 
      denominator(denominator_),
      dimension(std::max<int>(numberOfVariables(numerator_), numberOfVariables(denominator_)))
      {}

  RatioOfQsprays(int k)
    : numerator(scalarQspray<T>(T(k))), 
      denominator(scalarQspray<T>(T(1))),
      dimension(0)
      {}
  
  // methods -----
  std::pair<Qspray<T>,Qspray<T>> get() {
    return std::make_pair<Qspray<T>,Qspray<T>>(numerator, denominator);
  } 

  void simplify() {
  	Qspray<T> G = gcdQsprays(numerator, denominator);
  	numerator   = QuotientQsprays(numerator, G);
  	denominator = QuotientQsprays(denominator, G);
  }

  RatioOfQsprays<T> operator+=(const RatioOfQsprays<T>& ROQ2) {
  	numerator = numerator * ROQ2.denominator + ROQ2.numerator * denominator;
  	denominator = denominator * ROQ2.denominator;
  	RatioOfQsprays R(numerator, denominator);
  	R.simplify();
  	return R;
  }

};


#endif