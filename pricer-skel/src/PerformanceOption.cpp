/* 
 * File:   PerformanceOption.cpp
 * Author: samotd
 * 
 * Created on September 16, 2016, 11:06 AM
 */

#include "PerformanceOption.hpp"

PerformanceOption::PerformanceOption() {
}

PerformanceOption::PerformanceOption(const PerformanceOption& orig) {
}

PerformanceOption(double monT, int monNbTimeStep, int maSize) {
   T_ = monT;
   nbTimeSteps_ = monNbTimeStep;
   size_ = maSize;
}

PerformanceOption::~PerformanceOption() {
}

virtual double BasketOption::payoff(const PnlMat* path) {

    double sommeGenerale = 0;
    double sommeDenominateur = 0;
    double sommeNumerateur = 0;
    PnlMat* poids = pnl_mat_create(1, path->n);
    int N = (path->n) - 1;
    for (int j = 0; j < path->n; j++) {
        poids[0, j] = 1 / (path->n);
    }

    for (int i = 1; i < path->m; i++) {

        for (int j = 0; j < path->n; j++) {
            sommeNumerateur = poids[j] * path[i, j];
        }
        for (int j = 0; j < path->n; j++) {
            sommeNumerateur = poids[j] * path[i - 1, j];
        }
        sommeGenerale = sommeGenerale + (sommeNumerateur / sommeDenominateur) - 1;

        sommeNumerateur = 0;
        sommeDenominateur = 0;
    }

    if (sommeGenerale < 0) {
        sommeGenerale = 0;
    }

    sommeGenerale = sommeGenerale + 1;

    return sommeGenerale;
}