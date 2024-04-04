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


// ---------------------------------------------------------------------------//
template<typename T>
class RatioOfQsprays {

  Qspray<T> numerator;
  Qspray<T> denominator; 

public:
  // constructors -----
  RatioOfQsprays()
    : numerator(scalarQspray<T>(0)), denominator(scalarQspray<T>(1))
      {}

  RatioOfQsprays(const Qspray<T>& numerator_, const Qspray<T>& denominator_) 
    : numerator(numerator_), denominator(denominator_)
      {}

  RatioOfQsprays(int k)
    : numerator(scalarQspray<T>(k)), denominator(scalarQspray<T>(1))
      {}
  
  // methods -----
  std::pair<Qspray<T>,Qspray<T>> get() {
    return std::make_pair<Qspray<T>,Qspray<T>>(numerator, denominator);
  } 

};


#endif