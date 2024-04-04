#include "qspray.h"

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
        Qspray<gmpq> qtnt = quotient(LTp, LTg);
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
