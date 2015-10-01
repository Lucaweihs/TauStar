#include <Rcpp.h>
#include"red_black_tree.h"
#include<stdio.h>
#include<ctype.h>
#include<math.h>

using namespace Rcpp;

void DoubDest(void* a) {
	free((double*)a);
}
int DoubComp(const void* a,const void* b) {
	if( *(double*)a > *(double*)b) return(1);
	if( *(double*)a < *(double*)b) return(-1);
	return(0);
}
void DoubPrint(const void* a) {
	Rprintf("%f",*(double*)a);
}
void InfoPrint(void* a) { }
void InfoDest(void *a) { }
void bubbleSort(long double * vec, int n) {
	int i,j,k;
	long double tmp;
	for(i=n-1; i > 0; i--) {
		for(j = 0; j < i; j++) {
			if(vec[j+1] < vec[j]) {
				tmp = vec[j+1];
				vec[j+1] = vec[j];
				vec[j] = tmp;
			}
		}
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector TStarFastTiesRCPP(NumericVector xNumeric, NumericVector yNumeric) {
	int n = xNumeric.size();
	int i,j,k;
	rb_red_blk_tree* tree;
	rb_red_blk_tree* revTree;
	long double concord = 0;
	long double discord = 0;
	double minY;
	double maxY;
	int numYMin;
	int numYMax;
	int numGreaterMin;
	int numLessMax;
	int top, middle, bottom;
	long double c = 0;
	long double d = 0;
	int *savedYsInds;
	int numSavedYs = 0;
	int totalYInTree = 0;
	double lastX = 0;
	NumericVector tauStar;

	double *x;
	double *y;
	x = (double*) malloc(sizeof(double)*n);
	y = (double*) malloc(sizeof(double)*n);

	// Copy input data to something understood by the rbtree functions
	for(i=0; i<n; i++) {
		x[i] = xNumeric[i];
		y[i] = yNumeric[i];
	}

	savedYsInds = (int*) malloc(sizeof(int)*n);

	tree = RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);
	revTree = RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);

	for(i=0; i<n-1; i++) {
		if(lastX == x[i] && i != 0) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(k=0; k<numSavedYs; k++) {
				totalYInTree += 1;
				RBTreeInsert(tree,&(y[savedYsInds[k]]),0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		for(j=i+1; j<n; j++) {

			minY = fmin(y[i], y[j]);
			maxY = fmax(y[i], y[j]);

			numLessMax = RBNumLessThan(tree, &maxY);
			numGreaterMin = RBNumGreaterThan(tree, &minY);
			top = RBNumGreaterThan(tree, &maxY);
			bottom = RBNumLessThan(tree, &minY);

			if(minY != maxY) {
				middle = numGreaterMin + numLessMax - totalYInTree;
			} else {
				middle = 0;
			}

			concord += ((bottom*(bottom-1))/2) + ((top*(top-1))/2);

			if(minY != maxY) {
				discord += ((middle*(middle-1))/2) + top*middle + top*bottom + middle*bottom;
				numYMin = totalYInTree - numGreaterMin - bottom;
				numYMax = totalYInTree - top - numLessMax;
				discord += numYMin*numGreaterMin + numYMax*(numLessMax-numYMin);
			}
		}
	}

	// Now run everything in reverse to get rid of
	// quadruples incorrectly counted as discordant
	numSavedYs = 0;
	totalYInTree = 0;
	lastX = 0;

	for(i=n-1; i>0; i--) {
		if(lastX == x[i] && i != n-1) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(k=0; k<numSavedYs; k++) {
				totalYInTree += 1;
				RBTreeInsert(revTree,&(y[savedYsInds[k]]),0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		for(j=i-1; j>=0; j--) {

			minY = fmin(y[i], y[j]);
			maxY = fmax(y[i], y[j]);

			if(minY == maxY) {
				bottom = RBNumLessThan(revTree, &maxY);
				top = RBNumGreaterThan(revTree, &maxY);

				discord -= top*bottom;
			}
		}
	}

	RBTreeDestroy(tree);
	RBTreeDestroy(revTree);
	free(savedYsInds);
	free(x);
	free(y);

	c = 16*concord - 8*discord;
	d = (c < 0) ? -1 : 1;
	tauStar = NumericVector::create(d*expl(logl(d*c) - (logl(n) + logl(n-1) + logl(n-2) + logl(n-3))));
	return(tauStar);
}

// [[Rcpp::export]]
Rcpp::NumericVector VTStarFastTiesRCPP(NumericVector xNumeric, NumericVector yNumeric) {
	int n = xNumeric.size();
	int i,j,k;
	rb_red_blk_tree* tree;
	rb_red_blk_tree* revTree;
	long double concord = 0;
	long double discord = 0;
	double minY;
	double maxY;
	int numYMin;
	int numYMax;
	int numGreaterMin;
	int numLessMax;
	int top, middle, bottom;
	long double c = 0;
	long double d = 0;
	int *savedYsInds;
	int numSavedYs = 0;
	int totalYInTree = 0;
	double lastX = 0;
	NumericVector tauStar;

	double *x;
	double *y;
	x = (double*) malloc(sizeof(double)*n);
	y = (double*) malloc(sizeof(double)*n);

	// Copy input data to something understood by the rbtree functions
	for(i=0; i<n; i++) {
		x[i] = xNumeric[i];
		y[i] = yNumeric[i];
	}

	savedYsInds = (int*) malloc(sizeof(int)*n);

	tree = RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);
	revTree = RBTreeCreate(DoubComp, DoubDest, InfoDest, DoubPrint, InfoPrint);

	for(i=0; i<n; i++) {
		if(lastX == x[i] && i != 0) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(k=0; k<numSavedYs; k++) {
				totalYInTree += 1;
				RBTreeInsert(tree,&(y[savedYsInds[k]]),0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		top = RBNumGreaterThan(tree, &(y[i]));
		bottom = RBNumLessThan(tree, &(y[i]));
		concord += .5*((bottom*(bottom-1))/2) + .5*((top*(top-1))/2) + .25*(top+bottom);
		for(j=i+1; j<n; j++) {

			minY = fmin(y[i], y[j]);
			maxY = fmax(y[i], y[j]);

			numLessMax = RBNumLessThan(tree, &maxY);
			numGreaterMin = RBNumGreaterThan(tree, &minY);
			top = RBNumGreaterThan(tree, &maxY);
			bottom = RBNumLessThan(tree, &minY);

			if(minY != maxY) {
				middle = numGreaterMin + numLessMax - totalYInTree;
			} else {
				middle = 0;
			}

			concord += ((bottom*(bottom-1))/2) + ((top*(top-1))/2) + .5*(top+bottom);

			if(minY != maxY) {
				discord += ((middle*(middle-1))/2) + top*middle + top*bottom + middle*bottom;
				numYMin = totalYInTree - numGreaterMin - bottom;
				numYMax = totalYInTree - top - numLessMax;
				discord += numYMin*numGreaterMin + numYMax*(numLessMax-numYMin);
			}
		}
	}

	// Now run everything in reverse to get rid of
	// quadruples incorrectly counted as discordant
	numSavedYs = 0;
	totalYInTree = 0;
	lastX = 0;

	for(i=n-1; i>0; i--) {
		if(lastX == x[i] && i != n-1) {
			savedYsInds[numSavedYs] = i;
			numSavedYs += 1;
		} else {
			for(k=0; k<numSavedYs; k++) {
				totalYInTree += 1;
				RBTreeInsert(revTree,&(y[savedYsInds[k]]),0);
			}
			lastX = x[i];
			numSavedYs = 1;
			savedYsInds[0] = i;
		}
		for(j=i-1; j>=0; j--) {

			minY = fmin(y[i], y[j]);
			maxY = fmax(y[i], y[j]);

			if(minY == maxY) {
				bottom = RBNumLessThan(revTree, &maxY);
				top = RBNumGreaterThan(revTree, &maxY);

				discord -= top*bottom;
			}
		}
	}

	RBTreeDestroy(tree);
	RBTreeDestroy(revTree);
	free(savedYsInds);
	free(x);
	free(y);

	c = 16*concord - 8*discord;
	d = (c < 0) ? -1 : 1;
	tauStar = NumericVector::create(d*expl(logl(d*c) - 4*logl(n)));
	return(tauStar);
}