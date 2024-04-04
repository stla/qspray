/* Based on original code by Robin Hankin */

#include "qspray.h"

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
// Note: this class is not used currently
template<typename T>
class Qspray {
  const std::unordered_map<powers,T,PowersHasher> S;
public:
  Qspray()
    : S()
      {}

  Qspray(const std::unordered_map<powers,T,PowersHasher> &S_) 
    : S(S_) 
      {}
    
  bool isNull() {
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    bool result = true;
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      if(it->second != -it->second) {
        result = false;
        break;
      }
    } // insuffisant !!!!!!!!!!!
    return result;
  }

  bool operator==(const Qspray<T>& S2) {
    if(S.size() != S2.size()) {
      return false;
    }
    Qspray<T> Q(S2);
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      if(S[pows] != Q[pows]) {
        return false;
      } else {
        Q.erase(pows);
      }
    }
    // at this point, Q[k] == S[k] for every index 'k' of S;
    // it remains to check that every element of Q has been accounted for:
    if(Q.empty()) {
      return true;
    } else {
      return false;
    }
  }

  Qspray<T> operator-() {
    Qspray<T> Q(S);
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;  
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      Q[pows] = -it->second;
    }
    return Q;
  }

  Qspray<T> operator+(const Qspray<T>& S2) {
    Qspray<T> Q(S2);
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it;
    powers pows;
    const T zero(); // il s'agira de définir T() pour les ratio of polynomials :-)
    for(it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      T a = it->second + Q[pows]; // quid si pas de clé pows? https://stackoverflow.com/q/78272300/1100107
      if(a == zero) {
        Q.erase(pows); // il a peut-être pas de clé pows!
      } else {
        Q[pows] = a; 
      }
    }
    return Q;
  }

  Qspray<T> operator-(const Qspray<T>& S2) {
    return S + (-S2);
  }

  Qspray<T> operator*(const Qspray<T> S2) {
    Qspray<T> Sout;
    typename std::unordered_map<powers,T,PowersHasher>::const_iterator it1, it2;
    const T zero = T();
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
    return Sout;    
  }
  

};


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
Rcpp::List qspray_deriv(
    const Rcpp::List& Powers, const Rcpp::StringVector& coeffs,   
    const Rcpp::IntegerVector& n
){
  qspray S;
  powers v;
  signed int i, j, J, nj, expnt;
  signed int N = n.size();
  signed int nterms = coeffs.size();
  std::vector<std::vector<signed int>> Powers_out(nterms);
  std::vector<signed int> sizes(nterms);

  for(i = 0 ; i < nterms; i++) {
    Rcpp::IntegerVector Exponents = Powers(i);
    J = Exponents.size();
    sizes[i] = J;
    Powers_out[i].reserve(J);
    for(j = 0 ; j < J; j++){
      Powers_out[i].emplace_back(Exponents(j));
    }
  }
  
  for(i = 0; i < nterms; i++) {
    std::vector<signed int> exponents = Powers_out[i];
    J = sizes[i];
    if(J < N) {
      continue;
    }
    gmpq coeff(Rcpp::as<std::string>(coeffs(i)));
    for(j = 0; j < N; j++) {
      nj = n(j);
      while((nj > 0) && (coeff != 0)) { // while loop because it might not run at all
        coeff *= exponents[j]; // multiply coeff first, then decrement exponent 
        exponents[j]--;
        nj--;
      }
    }
    if(coeff != 0) {
      v.clear();
      for(j = 0; j < J; j++) {
        expnt = exponents[j];
        v.push_back(expnt);
      }
      simplifyPowers(v);
      S[v] += coeff;  // increment because v is not row-unique any more
    }
  }  // i loop closes
  
  return retval(S);
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_maker(const Rcpp::List& Powers,
                        const Rcpp::StringVector& coeffs) {
  return retval(prepare(Powers, coeffs));
}


// -------------------------------------------------------------------------- //
qspray add(const qspray& S1, const qspray& S2) {
  qspray S = S1;
  qspray::const_iterator it;
  powers pows;  
  for(it = S2.begin(); it != S2.end(); ++it) {
    pows = it->first;
    S[pows] += it->second;
    if(S[pows] == 0) {
      S.erase(pows);
    }
  }
  return S;
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_add(const Rcpp::List& Powers1,
                      const Rcpp::StringVector& coeffs1,
                      const Rcpp::List& Powers2,
                      const Rcpp::StringVector& coeffs2) {
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);
  qspray S = add(S1, S2);  
  return retval(S);
}


// -------------------------------------------------------------------------- //
qspray subtract(const qspray& S1, const qspray& S2) {
  qspray S = S1;
  qspray::const_iterator it;
  powers pows;  
  for(it = S2.begin(); it != S2.end(); ++it) {
    pows = it->first;
    const gmpq coeff = it->second;
    S[pows] -= coeff; 
    if(S[pows] == 0) {
      S.erase(pows);
    }
  }
  return(S);
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
  qspray S = subtract(S1, S2);
  return retval(S);
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


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
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

Rcpp::List leadingTerm(const qspray& S, int d) {
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

qspray quotient(Rcpp::List f, Rcpp::List g) {
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
  return S;
}

// [[Rcpp::export]]
Rcpp::List qsprayDivisionRcpp(
  Rcpp::List Powers1, Rcpp::StringVector coeffs1,
  Rcpp::List Powers2, Rcpp::StringVector coeffs2,
  int d
) {
  qspray p = makeQspray(Powers1, coeffs1);
  qspray g = makeQspray(Powers2, coeffs2);
  Rcpp::List LTg = leadingTerm(g, d);
  qspray q;
  qspray r;
  bool divoccured;
  while(!p.empty()) {
    divoccured = false;
    Rcpp::List LTp = leadingTerm(p, d);
    if(divides(LTg, LTp)) {
      qspray qtnt = quotient(LTp, LTg);
      p = subtract(p, prod(qtnt, g));
      q = add(q, qtnt);
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
      r = add(r, LTpspray);
      p = subtract(p, LTpspray);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("Q") = retval(q),
    Rcpp::Named("R") = retval(r)
  );
}

// [[Rcpp::export]]
Rcpp::List BBdivisionRcpp(
  Rcpp::List Powers, Rcpp::StringVector coeffs,
  Rcpp::List gs, Rcpp::List LTgs, int d
) {
  int ngs = gs.size();
  qspray p = makeQspray(Powers, coeffs);
  qspray r;
  int i;
  bool divoccured;
  while(!p.empty()) {
    i = 0;
    divoccured = false;
    Rcpp::List LTp = leadingTerm(p, d);
    while(i < ngs && !divoccured) {
      Rcpp::List LTg = LTgs(i);
      if(divides(LTg, LTp)) {
        Rcpp::List g  = gs(i);
        qspray gspray = makeQspray(g["powers"], g["coeffs"]);
        qspray qtnt = quotient(LTp, LTg);
        p = subtract(p, prod(qtnt, gspray));
        divoccured = true;
      } else {
        i++;
      }
    }
    if(!divoccured) {
      Rcpp::IntegerVector powsRcpp = LTp["powers"];
      std::string coeff = LTp["coeff"];
      gmpq coef(coeff);
      powers pows(powsRcpp.begin(), powsRcpp.end());
      simplifyPowers(pows);
      qspray LTpspray;
      LTpspray[pows] = coef;
      r = add(r, LTpspray);
      p = subtract(p, LTpspray);
    }
  }
  return retval(r);
}
