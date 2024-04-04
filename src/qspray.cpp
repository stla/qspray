/* Based on original code by Robin Hankin */

#include "qspray.h"

// -------------------------------------------------------------------------- //
qcplx qxmult(qcplx z1, qcplx z2) {
  gmpq r1 = z1.real();
  gmpq i1 = z1.imag();
  gmpq r2 = z2.real();
  gmpq i2 = z2.imag();
  qcplx result(r1*r2 - i1*i2, r1*i2 + r2*i1);
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
template <typename T>
Qspray<T> scalarQspray(T x) {
  typename std::unordered_map<powers,T,PowersHasher> singleton;
  powers pows(0);
  singleton[pows] = x;
  return Qspray<T>(singleton);
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
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);
  Qspray<gmpq> Q1(S1);
  Qspray<gmpq> Q2(S2);  
  Qspray<gmpq> Q = Q1 + Q2;
  return retval(Q.get());
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_subtract(const Rcpp::List& Powers1,
                           const Rcpp::StringVector& coeffs1,
                           const Rcpp::List& Powers2,
                           const Rcpp::StringVector& coeffs2) {
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);
  Qspray<gmpq> Q1(S1);
  Qspray<gmpq> Q2(S2);  
  Qspray<gmpq> Q = Q1 - Q2;
  return retval(Q.get());
}


// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List qspray_mult(const Rcpp::List& Powers1,
                       const Rcpp::StringVector& coeffs1,
                       const Rcpp::List& Powers2,
                       const Rcpp::StringVector& coeffs2) {
  // return retval(
  //     prod(makeQspray(Powers1, coeffs1), makeQspray(Powers2, coeffs2)));
  qspray S1 = makeQspray(Powers1, coeffs1);
  qspray S2 = makeQspray(Powers2, coeffs2);
  Qspray<gmpq> Q1(S1);
  Qspray<gmpq> Q2(S2);  
  Qspray<gmpq> Q = Q1 * Q2;
  return retval(Q.get());
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
  Qspray<gmpq> Q1(S1);
  Qspray<gmpq> Q2(S2);  
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

Rcpp::List leadingTerm(Qspray<gmpq>& Q, int d) {
  return leadingTerm(Q.get(), d);
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
  Qspray<gmpq> p(makeQspray(Powers1, coeffs1));
  Qspray<gmpq> g(makeQspray(Powers2, coeffs2));
  Rcpp::List LTg = leadingTerm(g, d);
  Qspray<gmpq> q;
  Qspray<gmpq> r;
  bool divoccured;
  while(!p.empty()) {
    divoccured = false;
    Rcpp::List LTp = leadingTerm(p, d);
    if(divides(LTg, LTp)) {
      qspray qtnt = quotient(LTp, LTg);
      Qspray<gmpq> Qtnt(qtnt);
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
  return Rcpp::List::create(
    Rcpp::Named("Q") = retval(q.get()),
    Rcpp::Named("R") = retval(r.get())
  );
}

// [[Rcpp::export]]
Rcpp::List BBdivisionRcpp(
  Rcpp::List Powers, Rcpp::StringVector coeffs,
  Rcpp::List gs, Rcpp::List LTgs, int d
) {
  int ngs = gs.size();
  Qspray<gmpq> p(makeQspray(Powers, coeffs));
  Qspray<gmpq> r;
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
        Qspray<gmpq> gspray(makeQspray(g["powers"], g["coeffs"]));
        Qspray<gmpq> qtnt(quotient(LTp, LTg));
        p -= qtnt * gspray;
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
      Qspray<gmpq> ltp(LTpspray);
      r += ltp;
      p -= ltp;
    }
  }
  return retval(r.get());
}
