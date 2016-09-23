#include "Market.hpp"
#include <iostream>



Market::Market() {
    this->sigma_ = pnl_vect_create(0);
    this->spot_ = pnl_vect_create(0);
    this->mu_ = pnl_vect_create(0);
    this->r_ = 0;
    maturity_ = 0;
    nbTimeSteps_ = 0;
    size_ = 0;
    CorrelationMat = pnl_mat_create(0,0);
    cholesky = 0;
}


Market::Market(PnlVect *sigma, PnlVect *spot, PnlVect *mu,  double rho, double maturity, int nbTimeSteps, int size, double r) {
    this->sigma_ = sigma;
    this->spot_ = spot;
    this->mu_ = mu;
    this->maturity_ = maturity;
    this->nbTimeSteps_ = nbTimeSteps;
    this->size_ = size;
    this->r_ = r;
    CorrelationMat = pnl_mat_create(size_,size_);
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            if (i == j) {
                pnl_mat_set(CorrelationMat,i,j,1);
            } else {
                pnl_mat_set(CorrelationMat,i,j,rho);

            }
        }
    }
    // Factorisation de cholesky de la matrice de correlation
    cholesky = pnl_mat_chol(CorrelationMat);

}
