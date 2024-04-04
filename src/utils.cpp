#include "qspray.h"

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
  while(zero && n > 0) {
    it--;
    n--;
    zero = pows[n] == 0;
  }
  if(zero) {
    pows = {};
  } else {
    pows.erase(it, pows.end());
  }
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

