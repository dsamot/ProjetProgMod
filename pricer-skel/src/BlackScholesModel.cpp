#include "BlackScholesModel.hpp"


BlackScholesModel::BlackScholesModel() {
    this->sigma_ = new PnlVect();
    this->spot_ = new PnlVect();
    size_ = 0;
    r_ = 0;
    rho_ = 0;
}

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot) {
    this->sigma_ = &sigma;
    this->spot_ = &spot;
    size_ = size;
    r_ = r;
    rho_ = rho;
}

void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng) {
    double simulBrownien;
    double pasTemps = T/(double) nbTimeSteps;
    double t;
    
    //Creation de la matrice de correlation
    PnlMat CorrelationMat = new PnlMat(size_);
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            if (i == j) {
                CorrelationMat[i,j] = sigma_[i];
            } else {
                CorrelationMat[i,j] = rho_;
            }
        }
    }
    
    // Factorisation de cholesky de la matrice de correlation
    PnlMat CholeskyMat = new PnlMat();
    CholeskyMat = pnl_mat_chol(CorrelationMat);
    
    // Remplissage de la premiere ligne de la matrice des chemins
    // avec les prix spot
    for (int d = 0; d < size_; d++) {
        path[0,d] = spot_[d];
    }
    
    PnlMat G = new PnlMat(size_,nbTimeSteps);
    for (int d = 0; d < size_; d++) {      
        for (int n = 1; n < nbTimeSteps + 1; n++) {
            G[d,n] = pnl_rng_normal(rng);
        }
    }
    
    double LG;
    PnlVect Ld;
    PnlVect Gn;
    
    for (int d = 0; d < size_; d++) {
        pnl_mat_get_row(Ld,CholeskyMat,d);       
        for (int n = 0; n < nbTimeSteps + 1; n++) {
            pnl_mat_get_col(Gn,G,n);
            LG = pnl_vect_scalar_prod(Ld,Gn);
            path[n+1,d] = path[n,d] * exp((r_ - (pow(sigma_[d],2)/2)) * pasTemps + sigma_[d] * sqrt(pasTemps) * LG);
        }
    }
    
}
