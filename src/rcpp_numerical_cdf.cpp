/***
 * Copyright (C) 2015 Luca Weihs
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

/***
 * A collection of functions that implement the main functionality of the
 * algorithms described in
 */

#include<RcppArmadillo.h>
//#include<stdio.h>
//#include<ctype.h>
#include<math.h>
#include<queue>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
std::complex<double> sinhProd(std::complex<double> v, int i) {
  std::complex<double> a = M_PI * pow(v, 0.5) / (1.0 * i);
  return(pow(sinh(a) / a, -.5));
}

//[[Rcpp::export]]
double hurwitzZeta(double exponent, double offset, int sumLength) {
  return(
    arma::accu(
      1.0 / arma::pow(offset + arma::linspace<arma::vec>(0, sumLength, sumLength + 1),
                    exponent))
  );
}

//[[Rcpp::export]]
double aCoef(int k, int h) {
  int sign = (k % 2 == 0) ? 1 : -1;
  return(sign * std::pow(hurwitzZeta(2 * k, h, 10000), 2.0) / (2 * k));
}

//[[Rcpp::export]]
std::complex<double> tailSum(std::complex<double> v, int h, int sumLength) {
  std::complex<double> sum = 0;
  for(int i = 1; i <= sumLength; i++) {
    sum += aCoef(i,h) * pow(v, i);
  }
  return sum;
}

//[[Rcpp::export]]
std::complex<double> gridSum(std::complex<double> v, int sideLen) {
  std::complex<double> sum = 0;
  for(int i = 1; i <= sideLen; i++) {
    for(int j = 1; j <= sideLen; j++) {
      sum += -0.5 * log(1.0 + v / pow(1.0 * i * j, 2.0));
    }
  }
  return(sum);
}

//[[Rcpp::export]]
std::complex<double> asymCharFunction(double t, double maxError) {
  if(t == 0) {
    return 1;
  }
  std::complex<double> v(0, 36 * (-2.0 * t) / pow(M_PI, 4.0));
  int h = ceil(pow(2.0 * abs(v), 1/4)) + 2;
  int l = ceil(
    -log2(maxError) +
      4 * log2(h) +
      2 * log2(M_PI * (6 * h - 5) / pow(6 * (1 - 2 * h), 2))
    ) + 2;
  std::complex<double> sum = -gridSum(v, h - 1);
  for (int i = 1; i <= h - 1; i++) {
    sum += 2.0 * log(sinhProd(v, i));
  }
  sum += tailSum(v, h, l);
  return(exp(sum));
}

class RiemannInt {
public:
  double lower;
  double upper;
  double funcError;
  std::complex<double> lowerEval;
  std::complex<double> upperEval;

  RiemannInt(double l, double u,
             std::complex<double> eval, bool isLowerEval,
             double errFromFunc) {
    lower = l;
    upper = u;
    if (isLowerEval) {
      lowerEval = eval;
      upperEval = asymCharFunction(u, errFromFunc / 2.0);
    } else {
      lowerEval = asymCharFunction(l, errFromFunc / 2.0);
      upperEval = eval;
    }
    funcError = errFromFunc;
  }

  RiemannInt(double l, double u, double errFromFunc) {
    lower = l;
    upper = u;
    lowerEval = asymCharFunction(l, errFromFunc / 2.0);
    upperEval = asymCharFunction(u, errFromFunc / 2.0);
    funcError = errFromFunc;
  }

  double width() {
    return upper - lower;
  }

  RiemannInt splitLower() {
    RiemannInt ri(lower, (upper - lower) / 2 + lower, lowerEval, true, funcError);
    return ri;
  }

  RiemannInt splitUpper() {
    RiemannInt ri((upper - lower) / 2 + lower, upper, upperEval, false, funcError);
    return ri;
  }
};

std::complex<double> asymCdfIntegrand(double x, double t, double maxError) {
  std::complex<double> val;
  std::complex<double> i(0, 1);
  if (t == 0) {
    val = x;
  } else {
    val = asymCharFunction(t, maxError) * (1.0 - exp(-i * t * x)) / (i * t);
  }
  return val / (2.0 * M_PI);
}

//[[Rcpp::export]]
double asymCdf(double x, double maxError) {
  double T = 215.0; // TODO: hardcoded for now
  double integrandError = pow(10.0, -16);

  int numInts = 200;
  double intWidth = T / numInts;
  std::vector<double> oldVec(numInts);
  double intVal = 0;
  for (int j = 0; j < numInts; j++) {
    if (j == 0) {
      oldVec[j] = asymCdfIntegrand(x, j * intWidth, integrandError).real();
    } else {
      oldVec[j] = 2.0 * asymCdfIntegrand(x, j * intWidth, integrandError).real();
    }
    intVal += oldVec[j];
  }
  intVal *= intWidth;

  double oldIntVal = std::numeric_limits<double>::infinity();

  int k = 0;
  while(fabs(intVal / oldIntVal - 1) > maxError * 0.1 && k++ < 1){
    oldIntVal = intVal;
    intVal = 0;
    numInts *= 2;
    intWidth /= 2;
    std::vector<double> newVec(numInts);
    for(int j = 0; j < numInts; j++) {
      if (j % 2 == 0) {
        newVec[j] = oldVec[j / 2];
      } else {
        newVec[j] = 2.0 * asymCdfIntegrand(x, j * intWidth, integrandError).real();
      }
      intVal += newVec[j];
    }
    intVal *= intWidth;
    oldVec = newVec;
  }

  return intVal;
}




