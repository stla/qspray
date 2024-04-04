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

#endif