#include "BlackScholesModel.hpp"
#include <iostream>


BlackScholesModel::BlackScholesModel() {
    this->sigma_ = new PnlVect();
    this->spot_ = new PnlVect();
    size_ = 0;
    r_ = 0;
    rho_ = 0;
}

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot) {
    this->sigma_ = sigma;
    this->spot_ = spot;
    size_ = size;
    r_ = r;
    rho_ = rho;
}

void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng) {
    double simulBrownien;
    double pasTemps = T/(double) nbTimeSteps;
    double t;
    
    //Creation de la matrice de correlation
    PnlMat *CorrelationMat = pnl_mat_create(size_,size_);
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            if (i == j) {
                pnl_mat_set(CorrelationMat,i,j,1);
                //std::cout << " valeur : " << pnl_mat_get(CorrelationMat,i,j) << std::endl;
                //CorrelationMat[i,j] = sigma_[i];
            } else {
                pnl_mat_set(CorrelationMat,i,j,rho_);
                //std::cout << " valeur : " << pnl_mat_get(CorrelationMat,i,j) << std::endl;
                //CorrelationMat[i,j] = rho_;
            }
        }
    }

    // Factorisation de cholesky de la matrice de correlation
    //PnlMat *CholeskyMat = pnl_mat_create(size_,size_);
    int cholesky;
    cholesky = pnl_mat_chol(CorrelationMat);
    
    
    // Remplissage de la premiere ligne de la matrice des chemins
    // avec les prix spot
    for (int d = 0; d < size_; d++) {
        pnl_mat_set(path,0,d,pnl_vect_get(spot_,d));
        //path[0,d] = spot_[d];
    }

    //PnlMat G = new PnlMat(size_,nbTimeSteps);
    PnlMat *G = pnl_mat_create(size_,nbTimeSteps);
    for (int d = 0; d < size_; d++) {      
        for (int n = 0; n < nbTimeSteps; n++) {
            pnl_mat_set(G,d,n,pnl_rng_normal(rng));
            //G[d,n] = pnl_rng_normal(rng);
        }
    }

    //pnl_mat_print(G);
    double LG;
    // A Initialiser
    PnlVect *Ld = pnl_vect_create(size_);
    PnlVect *Gn = pnl_vect_create(size_);
    double exprExp;
    
    for (int d = 0; d < size_; d++) {
        pnl_mat_get_row(Ld,CorrelationMat,d);       
        for (int n = 1; n < nbTimeSteps; n++) {
            pnl_mat_get_col(Gn,G,n);
            LG = pnl_vect_scalar_prod(Ld,Gn);
            exprExp = (r_ - (pow(pnl_vect_get(sigma_,d),2)/2)) * pasTemps + pnl_vect_get(sigma_,d) * sqrt(pasTemps) * LG;
            pnl_mat_set(path,n,d,(pnl_mat_get(path,n-1,d) * exp(exprExp)));
        }
    }

    /*for (int n = 0; n < nbTimeSteps; n++) {
        //std::cout << "sous jacent d=0, n : " << n << " valeur : " << pnl_mat_get(path,n,0) << std::endl;
    }*/
    
}
