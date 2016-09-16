/* 
 * File:   AsianOption.cpp
 * Author: samotd
 * 
 * Created on September 16, 2016, 10:40 AM
 */

#include "AsianOption.hpp"

AsianOption::AsianOption() {
}

AsianOption::AsianOption(const AsianOption& orig) {
}

AsianOption::AsianOption(double monT, int monNbTimeStep, int maSize, double monStrike){
   T_ = monT;
   nbTimeSteps_ = monNbTimeStep;
   size_ = maSize;
   strike = monStrike;
}

AsianOption::~AsianOption() {
}

virtual double BasketOption::payoff(const PnlMat* path) {
    double sommeGenerale = 0;
    double sommeInterne = 0;
    PnlMat* poids = pnl_mat_create(1, path->n);
    int N = (path->n) - 1;
    for (int j = 0; j < path->n; j++) {
        poids[0, j] = 1 / (path->n);
    }

    for (int j = 0; j < path->n; j++) {

        for (int i = 0; i < path->n; i++) {
            sommeInterne = path[i, j] - strike;
        }

        sommeGenerale = sommeGenerale + poids[j]*(1 / (path->n)) * sommeInterne;
        sommeInterne = 0;
    }

    if (sommeGenerale > 0) {
        return sommeGenerale;
    } else {
        return 0.0;
    }
    
}