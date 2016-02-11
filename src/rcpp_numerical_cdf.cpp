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
* A collection of functions that allow for the computation of the asymptotic
* distributions described in Nandy, Weihs, Drton (2016).
*/

#include<RcppArmadillo.h>
#include<algorithm>
#include<math.h>
#include<queue>
#include "AsymCdfIntegrandEvaluator.h"
#include "AsymPdfIntegrandEvaluator.h"
#include "AsymDiscreteCdfIntegrandEvaluator.h"
#include "AsymDiscretePdfIntegrandEvaluator.h"
#include "AsymMixedCdfIntegrandEvaluator.h"
#include "AsymMixedPdfIntegrandEvaluator.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void doubleWidth(arma::vec& positions, arma::vec& values,
                 IntegrandEvaluator& intEval, double x,
                 double integrandError) {
  int oldSize = values.size();
  double oldMaxPosition = positions(positions.size() - 1);
  double newMaxPosition = 2 * oldMaxPosition;
  double intervalSize = (newMaxPosition - oldMaxPosition) / oldSize;

  positions.resize(oldSize * 2);
  values.resize(oldSize * 2);
  for (int i = oldSize; i < 2 * oldSize; i++) {
    positions(i) = positions(i - 1) + intervalSize;
    values(i) = 2.0 * intEval.integrand(x, positions(i), integrandError).real();
  }
}

void bisect(arma::vec& positions, arma::vec& values,
            IntegrandEvaluator& intEval, double x,
            double integrandError) {
  int oldSize = values.size();

  positions.resize(oldSize * 2 - 1);
  values.resize(oldSize * 2 - 1);
  for (int i = oldSize - 1; i > 0; i--) {
    positions(2 * i) = positions(i);
    values(2 * i) = values(i);
  }
  for (int i = 1; i < positions.size(); i += 2) {
    positions(i) = (positions(i + 1) + positions(i - 1)) / 2.0;
    values(i) = 2.0 * intEval.integrand(x, positions(i), integrandError).real();
  }
}

double riemannIntegrate(const arma::vec& positions, const arma::vec& values) {
  if (positions(0) != 0 && positions.size() >= 2) {
    stop("riemannIntegrate expects the first position to be 0 and"
           " there must be at least 2 positions.");
  }
  double sum = 0;

  // Special case handling of 0th position
  double bisectVal = positions(1) / 2;
  sum += values(0) * (2 * bisectVal);

  // Looping through position 1,...,n - 2
  for (int i = 1; i <= positions.size() - 2; i++) {
    double lastBisectVal = bisectVal;
    bisectVal = (positions(i + 1) + positions(i)) / 2;
    double intervalWidth = (bisectVal - lastBisectVal);
    sum += values(i) * intervalWidth;
  }

  // Special case handling of last position
  int lastInd = values.size() - 1;
  sum += values(lastInd) * (2 * (positions(lastInd) - bisectVal));

  return sum;
}

/***
* Takes an integrand evaluator which corresponds to the integrand of a
* numerical characteristic function inversion and then performs the inversion
* using that integrand evaluator. The integration is done on an interval of
* size [-T,T] and is
*/
double numericalCfInversion(IntegrandEvaluator& intEval, double x, double T,
                            double convCrit, int maxIter) {
  // We implicitly assume that we work only with CDFs that are supported
  // strictly on the non-negative axis. As such we return 0 when x <= 0.
  if (x <= 0) {
    return 0;
  }
  double integrandError = convCrit * .0001;

  int numInts = 10;
  double intWidth = T / numInts;
  arma::vec positions(numInts);
  arma::vec values(numInts);

  for (int j = 0; j < numInts; j++) {
    positions(j) = j * intWidth;
    if (j == 0) {
      values(j) = intEval.integrand(x, positions(j), integrandError).real();
    } else {
      values(j) = 2.0*intEval.integrand(x, positions(j), integrandError).real();
    }
  }

  double oldIntVal = riemannIntegrate(positions, values);

  bisect(positions, values, intEval, x, integrandError);
  double intVal = riemannIntegrate(positions, values);
  double bisectChange = fabs(oldIntVal - intVal) + convCrit + 1;
  oldIntVal = intVal;

  doubleWidth(positions, values, intEval, x, integrandError);
  double widthChange = fabs(oldIntVal - intVal) + convCrit + 1;

  int k = 0;
  while (k < 5 || (fmax(bisectChange, widthChange) >= convCrit && k < maxIter)) {
    oldIntVal = intVal;
    if (bisectChange > widthChange) {
      bisect(positions, values, intEval, x, integrandError);
      intVal = riemannIntegrate(positions, values);
      bisectChange = fabs(oldIntVal - intVal);
    } else {
      doubleWidth(positions, values, intEval, x, integrandError);
      intVal = riemannIntegrate(positions, values);
      widthChange = fabs(oldIntVal - intVal);
    }
    k++;
  }

  if (k == maxIter) {
    Rprintf("WARNING: max iterations reached in"
              " characteristic function inversion.\n");
  }

  return intVal;
}

double boundInZeroOne(double x) {
  return fmin(fmax(x, 0), 1);
}

//[[Rcpp::export]]
arma::vec HoeffIndCdfRCPP(arma::vec x, double maxError) {
  AsymCdfIntegrandEvaluator acie;
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = boundInZeroOne(
      numericalCfInversion(acie, x[i], 50.0, maxError, 12)); // TODO: 50.0
                                                       // hardcoded for now
  }
  return cdfVals;
}

//[[Rcpp::export]]
arma::vec HoeffIndPdfRCPP(arma::vec x, double maxError) {
  AsymPdfIntegrandEvaluator apie;
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(
      numericalCfInversion(apie, x[i], 50.0, maxError, 12),
      0); // TODO: 100.0 hardcoded for now
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

// [[Rcpp::export]]
arma::vec HoeffIndDiscreteCdfRCPP(arma::vec x, arma::vec eigenP,
                                  arma::vec eigenQ, double maxError) {
  AsymDiscreteCdfIntegrandEvaluator adcie(eigenP, eigenQ);
  double k = 2 * (eigenQ.size() - 1) * (eigenP.size() - 1) + 1;
  double intervalWidth = std::min(std::max(pow((maxError / (2 * (k - 2))),
                                               1 / (1 - k / 2)) / 4, 100.0),
                                               4000000.0);
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = boundInZeroOne(numericalCfInversion(adcie, x[i], intervalWidth,
                                                     maxError, 17));
  }
  return cdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndDiscretePdfRCPP(arma::vec x, arma::vec eigenP,
                                  arma::vec eigenQ, double maxError) {
  AsymDiscretePdfIntegrandEvaluator adpie(eigenP, eigenQ);
  double k = 2 * (eigenQ.size() - 1) * (eigenP.size() - 1) + 1;
  double intervalWidth = std::min(std::max(pow((maxError / (2 * (k - 2))),
                                               1 / (1 - k / 2)) / 4, 100.0),
                                               4000000.0);
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion(adpie, x[i], 400.0, maxError, 17),
                      0);
  }
  return pdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndMixedCdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedCdfIntegrandEvaluator amcie(eigenP);
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = boundInZeroOne(
      numericalCfInversion(amcie, x[i], 20.0, maxError, 12));
  }
  return cdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndMixedPdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedPdfIntegrandEvaluator ampie(eigenP);
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion(ampie, x[i], 20.0, maxError, 12), 0);
  }
  return pdfVals;
}

