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
 

BasketOption::BasketOption(double monT, int monNbTimeStep, int maSize, double monStrike, PnlVect *payoffcoeff){
   T_ = monT;
   nbTimeSteps_ = monNbTimeStep;
   size_ = maSize;
   strike = monStrike;
   payoffCoeff_ = payoffcoeff;
}

BasketOption::~BasketOption() {
}

double BasketOption::payoff(const PnlMat* path) {
    double somme = 0;
    //PnlMat* poids = pnl_mat_create(1, path->n);
    //int N = (path->m) - 1;
    //std::cout << "nbTimesSteps : " << nbTimeSteps_  << "  N : " << N << std::endl;
    //double payoffcoeff = 1 / (double)(path->n);
    /*for (int j = 0; j < path->n; j++) {
        pnl_mat_set(poids,0,j,payoffcoeff);
    }*/

    for (int i = 0; i < path->n; i++) {
        somme = somme + pnl_vect_get(payoffCoeff_,i) * pnl_mat_get(path,nbTimeSteps_,i);
    }
    //pnl_mat_free(&poids);
    somme -= strike;
    if (somme > 0) {
        return somme;
    } else {
        return 0.0;
    }

} 

