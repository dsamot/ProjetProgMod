/* 
 * File:   BasketOption.cpp
 * Author: samotd
 * 
 * Created on September 16, 2016, 9:44 AM
 */

#include "BasketOption.h"
#include <iostream>

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

double BasketOption::payoff(const PnlMat* path) {
    double somme = 0;
    PnlMat* poids = pnl_mat_create(1, path->n);
    std::cout << "taille: " << path->n << std::endl;
    int N = (path->m) - 1;
    double payoffcoeff = 1 / (double)(path->n);
    for (int j = 0; j < path->n; j++) {
        pnl_mat_set(poids,0,j,payoffcoeff);
    }



    for (int i = 0; i < path->n; i++) {
        somme = somme + pnl_mat_get(poids,0,i) * pnl_mat_get(path,N,i);
    }
    somme -= strike;
    if (somme > 0) {
        return somme;
    } else {
        return 0.0;
    }
} 

