/* 
 * File:   PerformanceOption.cpp
 * Author: samotd
 * 
 * Created on September 16, 2016, 11:06 AM
 */

#include "PerformanceOption.hpp"
#include <iostream>

PerformanceOption::PerformanceOption() {
}

PerformanceOption::PerformanceOption(const PerformanceOption& orig) {
}

PerformanceOption::PerformanceOption(double monT, int monNbTimeStep, int maSize) {
   T_ = monT;
   nbTimeSteps_ = monNbTimeStep;
   size_ = maSize;
}

PerformanceOption::~PerformanceOption() {
}

 double PerformanceOption::payoff(const PnlMat* path) {

    double sommeGenerale = 0.0;
    double sommeDenominateur = 0.0;
    double sommeNumerateur = 0.0;
    double payoffcoeff = 1.0 / (double)(path->n);
    double p = 0;
    PnlMat* poids = pnl_mat_create(1, path->n);
    for (int j = 0; j < path->n; j++) {
          pnl_mat_set(poids,0,j,payoffcoeff);
    }

    for (int i = 1; i < path->m; i++) {

        for (int j = 0; j < path->n; j++) {
            sommeNumerateur += pnl_mat_get(poids,0,j) * pnl_mat_get(path,i,j);
        
            sommeDenominateur += pnl_mat_get(poids,0,j) * pnl_mat_get(path,i-1,j);
        }

        p = (sommeNumerateur / sommeDenominateur) - 1.0;
        
        if (p < 0) {
            sommeGenerale += 0;
        }
        else {
            sommeGenerale +=  p;
        }

        sommeNumerateur = 0.0;
        sommeDenominateur = 0.0;
    }

    sommeGenerale += 1;
    pnl_mat_free(&poids);

    return sommeGenerale;
}