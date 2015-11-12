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
#include<algorithm>
#include<math.h>
#include<queue>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

std::complex<double> sinhProd(std::complex<double> v, int i) {
  std::complex<double> a = M_PI * pow(v, 0.5) / (1.0 * i);
  return(pow(sinh(a) / a, -.5));
}

double hurwitzZeta(double exponent, double offset, double maxError) {
  double maxSumLength = pow(2, 27);
  double sumLengthAsDouble = fmax(
    ceil(pow(maxError, -1.0 / exponent)) - floor(offset) - 1, 15.0
  );
  int sumLength = (sumLengthAsDouble > maxSumLength) ? maxSumLength : sumLengthAsDouble;
  if (sumLengthAsDouble > maxSumLength) {
    Rprintf("WARNING: computation of hurwitz zeta may not provide desired level of accuracy.\n");
  }
  long double sum = 0;
  for (int i = 0; i <= sumLength; i++) {
    sum += 1.0 / pow(offset + i, exponent);
  }
  long double tailInt = (exponent - 1) * 1.0 / pow(sumLength + offset + 1, exponent - 1);
  return sum + tailInt;
}

double aCoef(int k, int h, double maxError) {
  double maxHurZeta = 1.0 / ((2 * k - 1) * pow(h - 1, 2 * k - 1));
  double maxErrorForHurZeta = -maxHurZeta +
    (1 / 2.0) * sqrt(4 * maxHurZeta * maxHurZeta + 8 * k * maxError);
  int sign = (k % 2 == 0) ? 1 : -1;
  return(sign * std::pow(hurwitzZeta(2 * k, h, maxErrorForHurZeta), 2.0) / (2 * k));
}

std::complex<double> tailSum(std::complex<double> v, int h, double maxError) {
  std::complex<double> sum = 0;
  std::complex<double> vProd = 1;
  double absV = abs(v);

  // Half of the max error comes from truncating the summation
  double factor = absV / pow(h, 4.0);
  int sumLength;
  if (factor >= 1) {
    Rprintf("WARNING: h chosen for tailSum is too small and may not result in inaccuracies. Choose h so that |v|/h^4 < 1 (best if < 1/2).");
    sumLength = 100;
  } else {
    sumLength = fmax(
      ceil((-log(maxError / 2.0) +
        4 * log(h) +
        2 * log(M_PI * (6 * h - 5) / pow(6 * (1 - 2 * h), 2))) / (-log(factor))) + 2,
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

std::complex<double> gridSum(std::complex<double> v, int sideLen) {
  std::complex<double> sum = 0;
  for(int i = 1; i <= sideLen; i++) {
    for(int j = 1; j <= sideLen; j++) {
      sum += -0.5 * log(1.0 + v / pow(1.0 * i * j, 2.0));
    }
  }
  return(sum);
}

std::complex<double> asymCharFunction(double t, double maxError) {
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

std::complex<double> asymCdfIntegrand(double x, double t, double maxError) {
  std::complex<double> val;
  std::complex<double> i(0, 1);
  if (t == 0) {
    val = x;
  } else {
    val = asymCharFunction(t, maxError * t / 2.0) * (1.0 - exp(-i * t * x)) / (i * t);
  }
  return val / (2.0 * M_PI);
}

std::complex<double> asymPdfIntegrand(double x, double t, double maxError) {
  std::complex<double> val;
  std::complex<double> i(0, 1);
  if (t == 0) {
    val = x;
  } else {
    val = asymCharFunction(t, maxError * t / 2.0) * exp(-i * t * x);
  }
  return val / (2.0 * M_PI);
}

//[[Rcpp::export]]
double HoeffIndCdfRCPP(double x, double maxError) {
  if (x <= 0) {
    return 0;
  }
  double T = 215.0; // TODO: hardcoded for now
  double integrandError = maxError * .0001;

  int numInts = 500;
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
  //while (fabs(intVal / oldIntVal - 1) > maxError) {
  while (fabs(intVal - oldIntVal) > maxError) {
    //printf("k=%d\n", k);
    //printf("intVal=%f\n", oldIntVal);
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
    if (k++ == 5) {
      Rprintf("WARNING: convergence criteria not met for asymptotic cdf computation.\n");
      break;
    }
  }

  return intVal;
}

//[[Rcpp::export]]
double HoeffIndPdfRCPP(double x, double maxError) {
  if (x <= 0) {
    return 0;
  }
  double T = 100.0; // TODO: hardcoded for now
  double integrandError = maxError * .0001;

  int numInts = 1000;
  double intWidth = T / numInts;
  std::vector<double> oldVec(numInts);
  double intVal = 0;
  for (int j = 0; j < numInts; j++) {
    if (j == 0) {
      oldVec[j] = asymPdfIntegrand(x, j * intWidth, integrandError).real();
    } else {
      oldVec[j] = 2.0 * asymPdfIntegrand(x, j * intWidth, integrandError).real();
    }
    intVal += oldVec[j];
  }
  intVal *= intWidth;

  double oldIntVal = std::numeric_limits<double>::infinity();

  int k = 0;
  while (fabs(intVal - oldIntVal) > maxError) {
    oldIntVal = intVal;
    intVal = 0;
    numInts *= 2;
    intWidth /= 2;
    std::vector<double> newVec(numInts);
    for(int j = 0; j < numInts; j++) {
      if (j % 2 == 0) {
        newVec[j] = oldVec[j / 2];
      } else {
        newVec[j] = 2.0 * asymPdfIntegrand(x, j * intWidth, integrandError).real();
      }
      intVal += newVec[j];
    }
    intVal *= intWidth;
    oldVec = newVec;
    if (k++ == 5) {
      Rprintf("WARNING: convergence criteria not met for asymptotic pdf computation.\n");
      break;
    }
  }
  return intVal;
}

// [[Rcpp::export]]
std::vector<std::complex<double> > eigenForDiscreteProbs(arma::vec p) {
  arma::vec cdf(p.size());
  arma::vec q(p.size());
  arma::vec q1(p.size());
  cdf[0] = p[0];
  q[0] = p[0] * cdf[0];
  for (int i = 1; i < p.size(); i++) {
    cdf[i] = cdf[i - 1] + p[i];
    q[i] = q[i-1] + p[i] * cdf[i];
  }
  q1[0] = p[0] * (1 - cdf[0]);
  for (int i = 1; i < p.size(); i++) {
    q1[i] = q1[i-1] + p[i] * (1 - cdf[i]);
  }

  arma::mat symMat(p.size(), p.size());
  for(int i = 0; i < p.size(); i++) {
    for(int j = i; j < p.size(); j++) {
      if (i == 0) {
        symMat(i,j) = -(1 - cdf[j]) * (1 - cdf[j]);
      } else {
        symMat(i,j) = -(1 - cdf[j]) * (1 - cdf[j]) - cdf[i-1] * cdf[i-1];
      }

      if (i != j) {
        symMat(i,j) += cdf[i] * (1 - cdf[i]) + (q1[j-1] - q1[i]);
        symMat(j,i) = symMat(i,j);
        symMat(j,i) *= p[i];
      }
      symMat(i,j) *= p[j];
    }
  }
  arma::Col<std::complex<double> > eigenVals = arma::eig_gen(symMat);
  std::vector<std::complex<double> > eigenValsVec(eigenVals.size());
  for (int i = 0; i < eigenVals.size(); i++) {
    eigenValsVec[i] = eigenVals[i];
  }
  return eigenValsVec;
}



