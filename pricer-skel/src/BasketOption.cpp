/* 
 * File:   BasketOption.cpp
 * Author: samotd
 * 
 * Created on September 16, 2016, 9:44 AM
 */

#include "BasketOption.h"

BasketOption::BasketOption() {
}

BasketOption::BasketOption(const BasketOption& orig) {
}
 double T_; /// maturité !
    int nbTimeSteps_; /// nombre de pas de temps de discrétisation
    int size_; 
BasketOption::BasketOption(double monT, int monNbTimeStep, int maSize, double monStrike){
   T_ = monT;
   nbTimeSteps_ = monNbTimeStep;
   size_ = maSize;
   strike = monStrike;
}

BasketOption::~BasketOption() {
}

virtual double BasketOption::payoff(const PnlMat* path) {
    double somme = 0;
    PnlMat* poids = pnl_mat_create(1, path->n);
    int N = (path->n) - 1;
    for (int j = 0; j < path->n; j++) {
        poids[0, j] = 1 / (path->n);
    }


    for (int i = 0; i < path->n; i++) {
        somme = somme + poids[i] * path[N, i] - strike;
    }

    if (somme > 0) {
        return somme;
    } else {
        return 0.0;
    }
} 

