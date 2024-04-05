#ifndef __HEADER__
#define __HEADER__

#include <Rcpp.h>
#include <boost/multiprecision/gmp.hpp>
#include <complex.h>
typedef std::vector<signed int>                             powers;
typedef boost::multiprecision::mpq_rational                 gmpq;
typedef boost::multiprecision::mpz_int                      gmpi;
typedef std::complex<gmpq>                                  qcplx;

std::string q2str(gmpq);
void simplifyPowers(powers&);
powers growPowers(powers, signed int, signed int);

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
  
  // methods ----------
  std::unordered_map<powers,T,PowersHasher> get() {
    return S;
  } 

  bool empty() {
    return S.empty();
  }

  bool isNull() {
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    T zero(0);
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

  int numberOfVariables() {
    int d = 0;
    for(const auto& term : S) {
      int n = term.first.size();
      if(n > d) {
        d = n;
      }
    }
    return d;
  }

  bool isConstant() {
    int nterms = S.size();
    bool result = false;
    if(nterms == 0) {
      result = true;
    } else if(nterms == 1) {
      powers pows(0);
      if(auto search = S.find(pows); search != S.end()) {
        result = true;
      }
    }
    return result;
  }

  T constantTerm() {
    powers pows(0);
    return S[pows];
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
    const T zero(0);
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
    const T zero(0);
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
    const T zero(0);
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
    T zero(0);
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

Qspray<gmpq> makeQspray(const Rcpp::List&, const Rcpp::StringVector&); 
Rcpp::List returnQspray(Qspray<gmpq>);


// -------------------------------------------------------------------------- //
template <typename T>
Qspray<T> scalarQspray(T x) {
  typename std::unordered_map<powers,T,PowersHasher> singleton;
  powers pows(0);
  singleton[pows] = x;
  return Qspray<T>(singleton);
}

static Qspray<gmpq> gcdQsprays(const Qspray<gmpq>& Q1, const Qspray<gmpq>& Q2) {
  return scalarQspray<gmpq>(1);
}

static int lexLeadingIndex(std::vector<powers> expnts) {
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

static Rcpp::List leadingTerm(Qspray<gmpq>& Q, int d) {
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

static bool divides(Rcpp::List f, Rcpp::List g) {
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

static Qspray<gmpq> quotient(Rcpp::List f, Rcpp::List g) {
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

static std::pair<Qspray<gmpq>,Qspray<gmpq>> qsprayDivision(
  Qspray<gmpq>& p, Qspray<gmpq>& g
) {
  int d = std::max<int>(p.numberOfVariables(), g.numberOfVariables());
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

static Qspray<gmpq> QuotientQsprays(Qspray<gmpq>& A, Qspray<gmpq>& B) {
  return qsprayDivision(A, B).first;
}

// ---------------------------------------------------------------------------//
template<typename T>
class RatioOfQsprays {

  Qspray<T> numerator;
  Qspray<T> denominator;
  int       dimension;

public:
  // constructors ---------------
  RatioOfQsprays()
    : numerator(scalarQspray<T>(T(0))), 
      denominator(scalarQspray<T>(T(1))),
      dimension(0)
      {}

  RatioOfQsprays(Qspray<T> numerator_, Qspray<T> denominator_) 
    : numerator(numerator_), 
      denominator(denominator_),
      dimension(
        std::max<int>(
          numerator_.numberOfVariables(), denominator_.numberOfVariables()
        )
      )
      {}

  RatioOfQsprays(int k)
    : numerator(scalarQspray<T>(T(k))), 
      denominator(scalarQspray<T>(T(1))),
      dimension(0)
      {}
  
  // methods --------------------
  Qspray<T> getNumerator() {
    return numerator;
  }

  Qspray<T> getDenominator() {
    return denominator;
  }

  void simplify() {
  	Qspray<T> G = gcdQsprays(numerator, denominator);
  	numerator   = QuotientQsprays(numerator, G);
  	denominator = QuotientQsprays(denominator, G);
    if(denominator.isConstant()) {
      Qspray<T> d = scalarQspray<T>(T(1) / denominator.constantTerm());
      numerator   *= d;
      denominator *= d;
    }
  }

  RatioOfQsprays<T> operator+=(const RatioOfQsprays<T>& ROQ2) {
  	numerator   = numerator * ROQ2.denominator + ROQ2.numerator * denominator;
  	denominator = denominator * ROQ2.denominator;
  	RatioOfQsprays ROQ(numerator, denominator);
  	ROQ.simplify();
  	return ROQ;
  }

  RatioOfQsprays<T> operator+(const RatioOfQsprays<T>& ROQ2) {
    RatioOfQsprays<T> ROQ(numerator, denominator);
    ROQ += ROQ2;
    return ROQ;
  }

  RatioOfQsprays<T> operator-=(const RatioOfQsprays<T>& ROQ2) {
    numerator   = numerator * ROQ2.denominator - ROQ2.numerator * denominator;
    denominator = denominator * ROQ2.denominator;
    RatioOfQsprays ROQ(numerator, denominator);
    ROQ.simplify();
    return ROQ;
  }

  RatioOfQsprays<T> operator-(const RatioOfQsprays<T>& ROQ2) {
    RatioOfQsprays<T> ROQ(numerator, denominator);
    ROQ -= ROQ2;
    return ROQ;
  }

  RatioOfQsprays<T> operator*=(const RatioOfQsprays<T>& ROQ2) {
    numerator   = numerator * ROQ2.numerator;
    denominator = denominator * ROQ2.denominator;
    RatioOfQsprays ROQ(numerator, denominator);
    ROQ.simplify();
    return ROQ;
  }

  RatioOfQsprays<T> operator*(const RatioOfQsprays<T>& ROQ2) {
    RatioOfQsprays<T> ROQ(numerator, denominator);
    ROQ *= ROQ2;
    return ROQ;
  }

  RatioOfQsprays<T> operator/=(const RatioOfQsprays<T>& ROQ2) {
    numerator   = numerator * ROQ2.denominator;
    denominator = denominator * ROQ2.numerator;
    RatioOfQsprays ROQ(numerator, denominator);
    ROQ.simplify();
    return ROQ;
  }

  RatioOfQsprays<T> operator/(const RatioOfQsprays<T>& ROQ2) {
    RatioOfQsprays<T> ROQ(numerator, denominator);
    ROQ /= ROQ2;
    return ROQ;
  }

};

// -------------------------------------------------------------------------- //
static Rcpp::List returnRatioOfQsprays(RatioOfQsprays<gmpq> ROQ) {
  return Rcpp::List::create(
    Rcpp::Named("numerator")   = returnQspray(ROQ.getNumerator()),
    Rcpp::Named("denominator") = returnQspray(ROQ.getDenominator())
  );
}


// -------------------------------------------------------------------------- //
static RatioOfQsprays<gmpq> makeRatioOfQsprays(
  const Rcpp::List& Numerator, 
  const Rcpp::List& Denominator
) {
  Rcpp::List Powers1 = Numerator["powers"];
  Rcpp::List Powers2 = Denominator["powers"];
  Rcpp::StringVector coeffs1 = Numerator["coeffs"];
  Rcpp::StringVector coeffs2 = Denominator["coeffs"];
  qspray S1;
  for(int i = 0; i < Powers1.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers1(i);
    gmpq coeff(Rcpp::as<std::string>(coeffs1(i)));
    powers pows(Exponents.begin(), Exponents.end());
    S1[pows] = coeff;
  }
  qspray S2;
  for(int i = 0; i < Powers2.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers2(i);
    gmpq coeff(Rcpp::as<std::string>(coeffs2(i)));
    powers pows(Exponents.begin(), Exponents.end());
    S2[pows] = coeff;
  }
  Qspray<gmpq> Q1(S1);
  Qspray<gmpq> Q2(S2);
  RatioOfQsprays<gmpq> ROQ(Q1, Q2);
  return ROQ;
}

// -------------------------------------------------------------------------- //
typedef Qspray<RatioOfQsprays<gmpq>>                                   SymbolicQspray;
typedef std::unordered_map<powers, RatioOfQsprays<gmpq>, PowersHasher> symbolicQspray;


static Rcpp::List returnSymbolicQspray(SymbolicQspray SQ) { // used to return a list to R
  symbolicQspray S = SQ.get();
  if(S.size() == 0) {
    return Rcpp::List::create(Rcpp::Named("powers") = R_NilValue,
                              Rcpp::Named("coeffs") = R_NilValue);
  } else {
    Rcpp::List Powers(S.size());
    powers pows;
    unsigned int row = 0, col = 0;
    Rcpp::List Coeffs(S.size());
    unsigned int i = 0;
    for(auto it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      Rcpp::IntegerVector Exponents(pows.size());
      col = 0;
      for(auto ci = pows.begin(); ci != pows.end(); ++ci) {
        Exponents(col++) = *ci;
      }
      Powers(row++) = Exponents;
      Coeffs(i++) = returnRatioOfQsprays(it->second);
    }
    return Rcpp::List::create(Rcpp::Named("powers") = Powers,
                              Rcpp::Named("coeffs") = Coeffs);
  }
}

// -------------------------------------------------------------------------- //
static SymbolicQspray makeSymbolicQspray(
  const Rcpp::List& Powers, const Rcpp::List& Coeffs
) {
  symbolicQspray S;
  for(int i = 0; i < Powers.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers(i);
    powers pows(Exponents.begin(), Exponents.end());
    Rcpp::List coeff = Coeffs(i);
    Rcpp::List numerator   = coeff["numerator"];
    Rcpp::List denominator = coeff["denominator"];
    S[pows] = makeRatioOfQsprays(numerator, denominator);
  }
  return SymbolicQspray(S);
}



#endif