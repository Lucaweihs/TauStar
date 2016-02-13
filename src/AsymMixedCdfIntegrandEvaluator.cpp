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

#include "AsymMixedCdfIntegrandEvaluator.h"
using namespace Rcpp;

typedef AsymMixedCdfIntegrandEvaluator AMCIE;

AMCIE::AsymMixedCdfIntegrandEvaluator(arma::vec eigP): eigenP(eigP) {}

int piRemSign(double x) {
  if (x == 0) {
    return 0;
  } else if (x > 0) {
    return (fmod(x, 2 * M_PI) <= M_PI) ? 1 : -1;
  } else {
    return (fmod(x, 2 * M_PI) >= -M_PI) ? 1 : -1;
  }
}

int getSinhSign(double rate) {
  int j = 0;
  double sum = 0;
  double remainder = 0.5 * rate * M_PI * M_PI / 6.0;
  while (fabs(remainder) >= M_PI ||
         (piRemSign(sum) != piRemSign(sum + remainder))) {
    j++;
    double v = rate / ((1.0 * j) * j);
    sum += 0.5 * asin(v / sqrt(1 + v * v));
    remainder -= 0.5 * v;
    if (j % 10000 == 0) {
      break;
    }
  }
  return piRemSign(sum);
}

std::complex<double> AMCIE::integrand(double x, double t, double maxError) {
  if (t == 0) {
    return x / (2 * M_PI);
  }
  std::complex<double> val;
  std::complex<double> I(0, 1);

  std::complex<double> sum = 0;
  std::complex<double> v(0, 12.0 * (-2.0 * t) / (M_PI * M_PI));
  double precision = pow(10, -15);
  for(int i = 0; i < eigenP.size(); i++) {
    if (fabs(eigenP[i]) > precision) {
      int sign = getSinhSign((v * eigenP[i]).imag());
      std::complex<double> sinhProdVal = sinhProd(v * eigenP[i], 1);
      if (sinhProdVal.imag() * sign <= 0) {
        sinhProdVal *= -1;
      }
      sum += log(sinhProdVal);
    }
  }
  return 1 / (2 * M_PI) * exp(sum) * (1.0 - exp(-I * t * x)) / (I * t);
}
