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
      oldVec[j] =
        2.0 * intEval->integrand(x, j * intWidth, integrandError).real();
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
        newVec[j] =
          2.0 * intEval->integrand(x, j * intWidth, integrandError).real();
      }
      intVal += newVec[j];
    }
    intVal *= intWidth;
    oldVec = newVec;
    if (k++ == maxIter) {
      Rprintf("WARNING: convergence criteria not met for"
                " characteristic function inversion.\n");
      break;
    }
  }

  return intVal;
}

double boundInZeroOne(double x) {
  return fmin(fmax(x, 0), 1);
}

//[[Rcpp::export]]
arma::vec HoeffIndCdfRCPP(arma::vec x, double maxError) {
  AsymCdfIntegrandEvaluator acie;
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = boundInZeroOne(
      numericalCfInversion((IntegrandEvaluator*) &acie,
                           x[i], 215.0, maxError, 7)); // TODO: 215.0
                                                       // hardcoded for now
  }
  return pdfVals;
}

//[[Rcpp::export]]
arma::vec HoeffIndPdfRCPP(arma::vec x, double maxError) {
  AsymPdfIntegrandEvaluator apie;
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(
      numericalCfInversion((IntegrandEvaluator*) &apie, x[i], 100.0,
                           maxError, 7),
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
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = boundInZeroOne(
      numericalCfInversion((IntegrandEvaluator*) &adcie, x[i], 215.0, maxError,
                           7));
  }
  return cdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndDiscretePdfRCPP(arma::vec x, arma::vec eigenP,
                                  arma::vec eigenQ, double maxError) {
  AsymDiscretePdfIntegrandEvaluator adpie(eigenP, eigenQ);
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion((IntegrandEvaluator*) &adpie, x[i],
                                           100.0, maxError, 7), 0);
  }
  return pdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndMixedCdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedCdfIntegrandEvaluator amcie(eigenP);
  arma::vec cdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    cdfVals[i] = boundInZeroOne(
      numericalCfInversion((IntegrandEvaluator*) &amcie, x[i], 215.0, maxError,
                           10));
  }
  return cdfVals;
}

// [[Rcpp::export]]
arma::vec HoeffIndMixedPdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedPdfIntegrandEvaluator ampie(eigenP);
  arma::vec pdfVals(x.size());
  for (int i = 0; i < x.size(); i++) {
    pdfVals[i] = fmax(numericalCfInversion((IntegrandEvaluator*) &ampie, x[i],
                                           215.0, maxError, 10), 0);
  }
  return pdfVals;
}

