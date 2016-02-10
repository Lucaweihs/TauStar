/***
 * Copyright (C) 2016 Luca Weihs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "AsymCdfIntegrandEvaluator.h"
#include "RcppArmadillo.h"
using namespace Rcpp;

/***
 * Computes the infinite summation
 *
 * \sum_{i=0}^{\infty} 1 / (offset + i)^exponent
 *
 * up to some specified max error. If the function cannot guarentee the desired
 * error then a warning is printed.
 */
double hurwitzZeta(double exponent, double offset, double maxError) {
  if (exponent <= 1.0) {
    Rcpp::stop("Exponent in hurwitzZeta must be > 1");
  }
  double maxSumLength = pow(2, 27);
  double sumLengthAsDouble = fmax(
    ceil(pow(maxError, -1.0 / exponent)) - floor(offset) - 1, 15.0
  );
  int sumLength
    = (sumLengthAsDouble > maxSumLength) ?maxSumLength : sumLengthAsDouble;
  if (sumLengthAsDouble > maxSumLength) {
    Rprintf("WARNING: computation of hurwitz zeta may not provide"
              "desired level of accuracy.\n");
  }
  long double sum = 0;
  for (int i = 0; i <= sumLength; i++) {
    sum += 1.0 / pow(offset + i, exponent);
  }
  long double tailInt = (exponent - 1) * 1.0 / pow(sumLength + offset + 1,
                         exponent - 1);
  return sum + tailInt;
}

/***
 * A function used in calculating the characteristic function for the
 * asymptotic distribution t* statistic in the continuous case. In particular
 * this computes a simple sum of complex values along a grid.
 */
std::complex<double> gridSum(std::complex<double> v, int sideLen) {
  std::complex<double> sum = 0;
  for(int i = 1; i <= sideLen; i++) {
    for(int j = 1; j <= sideLen; j++) {
      sum += -0.5 * log(1.0 + v / pow(1.0 * i * j, 2.0));
    }
  }
  return(sum);
}

/***
 * Computes one of the square roots of the infinite product
 *
 * \prod_{j=1}^{\infty} (1 + v / (i^2 * j^2))^(-1)
 *
 * which square root is returned is undefined.
 */
std::complex<double> sinhProd(std::complex<double> v, int i) {
  std::complex<double> a = M_PI * sqrt(v) / (1.0 * i);
  return(sqrt(a / sinh(a)));
}

/***
 * Computes
 *
 * (-1)^k \zeta(2k, h)^2 / (2k)
 *
 * where \zeta is the Hurwitz-Zeta function. This computation is only up to some
 * max specified error.
 */
double aCoef(int k, int h, double maxError) {
  double maxHurZeta = 1.0 / ((2 * k - 1) * pow(h - 1, 2 * k - 1));
  double maxErrorForHurZeta = -maxHurZeta +
    (1 / 2.0) * sqrt(4 * maxHurZeta * maxHurZeta + 8 * k * maxError);
  int sign = (k % 2 == 0) ? 1 : -1;
  return(sign * std::pow(hurwitzZeta(2 * k, h,
                                     maxErrorForHurZeta), 2.0) / (2 * k));
}

/***
 * A function used in calculating the characteristic function for the
 * asymptotic distribution t* statistic in the continuous case. In particular
 * this is the infinite tail summation of that computation. As usualy this is
 * computed only up to some specified error.
 */
std::complex<double> tailSum(std::complex<double> v, int h, double maxError) {
  std::complex<double> sum = 0;
  std::complex<double> vProd = 1;
  double absV = abs(v);

  // Half of the max error comes from truncating the summation
  double factor = absV / pow(h, 4.0);
  int sumLength;
  if (factor >= 1) {
    Rprintf("WARNING: h chosen for tailSum is too small and may not result in"
              "inaccuracies. Choose h so that |v|/h^4 < 1 (best if < 1/2).");
    sumLength = 100;
  } else {
    sumLength = fmax(
      ceil((-log(maxError / 2.0) +
        4 * log(h) +
        2 * log(M_PI * (6 * h - 5) / pow(6 * (1 - 2 * h), 2)))
             / (-log(factor))) + 2,
             10);
  }

  // The second half of the error comes from approximating the a coefficients
  // in the summation. we want 1/4 the error from the first term, 1/8 for the
  // second, 1/16 from the third, ... This follows the empirical observation
  // that attaining small error rates from the first term takes many many terms.
  double maxErrorForA = maxError / 2;
  for(int i = 1; i <= sumLength; i++) {
    vProd *= v;
    maxErrorForA /= 2.0 * absV;
    sum += aCoef(i, h, maxErrorForA) * vProd;
  }
  return sum;
}

/***
 * The characteristic function of the asymptotic distribution of t* in the
 * continuous case computed up to some specified error.
 */
std::complex<double> asymContCharFunction(double t, double maxError) {
  if(t == 0) {
    return 1;
  }
  std::complex<double> v(0, 36 * (-2.0 * t) / pow(M_PI, 4.0));
  int h = ceil(pow(2.0 * abs(v), 1.0 / 4.0)) + 2;
  std::complex<double> sum = -gridSum(v, h - 1);
  for (int i = 1; i <= h - 1; i++) {
    sum += 2.0 * log(sinhProd(v, i));
  }
  sum += tailSum(v, h, maxError / abs(exp(sum)));
  return(exp(sum));
}


std::complex<double> AsymCdfIntegrandEvaluator::integrand(double x, double t, double maxError) {
  std::complex<double> val;
  std::complex<double> i(0, 1);
  if (t == 0) {
    val = x;
  } else {
    val = asymContCharFunction(t, maxError * t / 2.0) *
      (1.0 - exp(-i * t * x)) / (i * t);
  }
  return val / (2.0 * M_PI);
}
