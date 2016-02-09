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
  return(sign * std::pow(hurwitzZeta(2 * k, h, maxErrorForHurZeta), 2.0) / (2 * k));
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
 * The characteristic function of the asymptotic distribution of t* in the
 * continuous case computed up to some specified error.
 */
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

/***
 * A class whose only purpose is to be subclassed. Provides an interface
 * to get an integrand (taking two parameters) up to some maximum error.
 */
class IntegrandEvaluator {
public:
  virtual std::complex<double>
    integrand(double x, double t, double maxError) = 0;
};

/***
 * An IntegrandEvaluator used when numerically inverting the characteristic
 * function of the asymptotic distribution of t* to get the CDF of that
 * distribution.
 */
class AsymCdfIntegrandEvaluator : public IntegrandEvaluator {
public:
  std::complex<double> integrand(double x, double t, double maxError) {
    std::complex<double> val;
    std::complex<double> i(0, 1);
    if (t == 0) {
      val = x;
    } else {
      val = asymCharFunction(t, maxError * t / 2.0) * (1.0 - exp(-i * t * x)) / (i * t);
    }
    return val / (2.0 * M_PI);
  }
};

/***
 * An IntegrandEvaluator used when numerically inverting the characteristic
 * function of the asymptotic distribution of t* to get the PDF of that
 * distribution.
 */
class AsymPdfIntegrandEvaluator : public IntegrandEvaluator {
public:
  std::complex<double> integrand(double x, double t, double maxError) {
    std::complex<double> val;
    std::complex<double> i(0, 1);
    if (t == 0) {
      val = x;
    } else {
      val = asymCharFunction(t, maxError * t / 2.0) * exp(-i * t * x);
    }
    return val / (2.0 * M_PI);
  }
};

/***
 * Takes an integrand evaluator which corresponds to the integrand of a
 * numerical characteristic function inversion and then performs the inversion
 * using that integrand evaluator. The integration is done on an interval of
 * size [-T,T] and is
 */
double numericalCfInversion(
    IntegrandEvaluator* intEval,
    double x,
    double T,
    double convCrit,
    int maxIter
) {
  if (x <= 0) {
    return 0;
  }
  double integrandError = convCrit * .0001;

  int numInts = 500;
  double intWidth = T / numInts;
  std::vector<double> oldVec(numInts);
  double intVal = 0;
  for (int j = 0; j < numInts; j++) {
    if (j == 0) {
      oldVec[j] = intEval->integrand(x, j * intWidth, integrandError).real();
    } else {
      oldVec[j] = 2.0 * intEval->integrand(x, j * intWidth, integrandError).real();
    }
    intVal += oldVec[j];
  }
  intVal *= intWidth;

  double oldIntVal = std::numeric_limits<double>::infinity();

  int k = 0;
  while (fabs(intVal - oldIntVal) > convCrit) {
    oldIntVal = intVal;
    intVal = 0;
    numInts *= 2;
    intWidth /= 2;
    std::vector<double> newVec(numInts);
    for(int j = 0; j < numInts; j++) {
      if (j % 2 == 0) {
        newVec[j] = oldVec[j / 2];
      } else {
        newVec[j] = 2.0 * intEval->integrand(x, j * intWidth, integrandError).real();
      }
      intVal += newVec[j];
    }
    intVal *= intWidth;
    oldVec = newVec;
    //Rprintf("Error: %f Int val: %f\n", fabs(intVal - oldIntVal), intVal);
    if (k++ == maxIter) {
      Rprintf("WARNING: convergence criteria not met for characteristic function inversion.\n");
      break;
    }
  }

  return intVal;
}

//[[Rcpp::export]]
arma::vec HoeffIndCdfRCPP(arma::vec x, double maxError) {
  AsymCdfIntegrandEvaluator acie;
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmin(fmax(numericalCfInversion((IntegrandEvaluator*) &acie, x[i], 215.0, maxError, 7), 0), 1); // TODO: 215.0 hardcoded for now
  }
  return pdfVals;
}

//[[Rcpp::export]]
arma::vec HoeffIndPdfRCPP(arma::vec x, double maxError) {
  AsymPdfIntegrandEvaluator apie;
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion((IntegrandEvaluator*) &apie, x[i], 100.0, maxError, 7), 0); // TODO: 100.0 hardcoded for now
  }
  return pdfVals;
}

// [[Rcpp::export]]
arma::vec eigenForDiscreteProbs(arma::vec p) {
  arma::vec cdf(p.size());
  arma::vec q1(p.size());
  cdf[0] = p[0];
  for (int i = 1; i < p.size(); i++) {
    cdf[i] = cdf[i - 1] + p[i];
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
        symMat(j,i) *= sqrt(p[i] * p[j]);
      }
      symMat(i,j) *= sqrt(p[i] * p[j]);
    }
  }
  return arma::eig_sym(symMat);
}

class AsymDiscreteCdfIntegrandEvaluator : public IntegrandEvaluator {
protected:
  arma::vec eigenP;
  arma::vec eigenQ;
public:
  AsymDiscreteCdfIntegrandEvaluator(arma::vec eigP, arma::vec eigQ) {
    eigenP = eigP;
    eigenQ = eigQ;
  }

  std::complex<double> integrand(double x, double t, double maxError) {
    if (t == 0) {
      return x / (2 * M_PI);
    }
    std::complex<double> val;
    std::complex<double> I(0, 1);

    std::complex<double> sum = 0;
    for(int i = 0; i < eigenP.size(); i++) {
      for(int j = 0; j < eigenQ.size(); j++) {
        sum += -0.5 * log(1.0 - 8.0 * I * t * eigenP[i] * eigenQ[j]);
      }
    }
    return 1 / (2 * M_PI) * exp(sum) * (1.0 - exp(-I * t * x)) / (I * t);
  }
};

class AsymDiscretePdfIntegrandEvaluator : public IntegrandEvaluator {
protected:
  arma::vec eigenP;
  arma::vec eigenQ;
public:
  AsymDiscretePdfIntegrandEvaluator(arma::vec eigP, arma::vec eigQ) {
    eigenP = eigP;
    eigenQ = eigQ;
  }

  std::complex<double> integrand(double x, double t, double maxError) {
    if (t == 0) {
      return x / (2 * M_PI);
    }
    std::complex<double> val;
    std::complex<double> I(0, 1);

    std::complex<double> sum = 0;
    for(int i = 0; i < eigenP.size(); i++) {
      for(int j = 0; j < eigenQ.size(); j++) {
        sum += -0.5 * log(1.0 - 8.0 * I * t * eigenP[i] * eigenQ[j]);
      }
    }
    return 1 / (2 * M_PI) * exp(sum) * exp(-I * t * x);
  }
};

// [[Rcpp::export]]
arma::vec HoeffIndDiscreteCdfRCPP(arma::vec x, arma::vec eigenP, arma::vec eigenQ,
                                  double maxError) {
  AsymDiscreteCdfIntegrandEvaluator adcie(eigenP, eigenQ);
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = fmin(fmax(numericalCfInversion((IntegrandEvaluator*) &adcie, x[i], 215.0, maxError, 7), 0), 1);
  }
  return cdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndDiscretePdfRCPP(arma::vec x, arma::vec eigenP, arma::vec eigenQ,
                                  double maxError) {
  AsymDiscretePdfIntegrandEvaluator adpie(eigenP, eigenQ);
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion((IntegrandEvaluator*) &adpie, x[i], 100.0, maxError, 7), 0);
  }
  return pdfVals;
}

class AsymMixedCdfIntegrandEvaluator : public IntegrandEvaluator {
protected:
  arma::vec eigenP;

public:
  AsymMixedCdfIntegrandEvaluator(arma::vec eigP) {
    eigenP = eigP;
  }
  static int piRemSign(double x) {
    if (x == 0) {
      return 0;
    } else if (x > 0) {
      return (fmod(x, 2 * M_PI) <= M_PI) ? 1 : -1;
    } else {
      return (fmod(x, 2 * M_PI) >= -M_PI) ? 1 : -1;
    }
  }

  static int getSinhSign(double rate) {
    int j = 0;
    double sum = 0;
    double remainder = 0.5 * rate * M_PI * M_PI / 6.0;
    while (fabs(remainder) >= M_PI || (piRemSign(sum) != piRemSign(sum + remainder))) {
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

  std::complex<double> integrand(double x, double t, double maxError) {
    if (t == 0) {
      return x / (2 * M_PI);
    }
    std::complex<double> val;
    std::complex<double> I(0, 1);

    std::complex<double> sum = 0;
    std::complex<double> v(0, 12.0 * (-2.0 * t) / (M_PI * M_PI));
    double precision = pow(10, -15);
    for(int i = 0; i < eigenP.size() && fabs(eigenP[i]) > precision; i++) {
      int sign = getSinhSign((v * (-eigenP[i])).imag());
      std::complex<double> sinhProdVal = sinhProd(v * (-eigenP[i]), 1);
      if (sinhProdVal.imag() * sign <= 0) {
        sinhProdVal *= -1;
      }
      sum += log(sinhProdVal);
    }
    return 1 / (2 * M_PI) * exp(sum) * (1.0 - exp(-I * t * x)) / (I * t);
  }
};

class AsymMixedPdfIntegrandEvaluator : public IntegrandEvaluator {
protected:
  arma::vec eigenP;

public:
  AsymMixedPdfIntegrandEvaluator(arma::vec eigP) {
    eigenP = eigP;
  }

  std::complex<double> integrand(double x, double t, double maxError) {
    if (t == 0) {
      return x / (2 * M_PI);
    }
    std::complex<double> val;
    std::complex<double> I(0, 1);

    std::complex<double> sum = 0;
    std::complex<double> v(0, 12.0 * (-2.0 * t) / (M_PI * M_PI));
    double precision = pow(10, -15);
    for(int i = 0; i < eigenP.size() && fabs(eigenP[i]) > precision; i++) {
      int sign = AsymMixedCdfIntegrandEvaluator::getSinhSign((v * (-eigenP[i])).imag());
      std::complex<double> sinhProdVal = sinhProd(v * (-eigenP[i]), 1);
      if (sinhProdVal.imag() * sign <= 0) {
        sinhProdVal *= -1;
      }
      sum += log(sinhProdVal);
    }
    return 1 / (2 * M_PI) * exp(sum) * exp(-I * t * x);
  }
};

// [[Rcpp::export]]
arma::vec HoeffIndMixedCdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedCdfIntegrandEvaluator amcie(eigenP);
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = fmin(fmax(numericalCfInversion((IntegrandEvaluator*) &amcie, x[i], 215.0, maxError, 10), 0), 1);
  }
  return cdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndMixedPdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedPdfIntegrandEvaluator ampie(eigenP);
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion((IntegrandEvaluator*) &ampie, x[i], 215.0, maxError, 10), 0);
  }
  return pdfVals;
}



