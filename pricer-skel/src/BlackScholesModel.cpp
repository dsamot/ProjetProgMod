#include "BlackScholesModel.hpp"
#include <iostream>


BlackScholesModel::BlackScholesModel() {
    this->sigma_ = pnl_vect_create(0);
    this->spot_ = pnl_vect_create(0);
    size_ = 0;
    r_ = 0;
    rho_ = 0;
    CorrelationMat = pnl_mat_create(0,0);
    cholesky = 0;
}

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot) {
    this->sigma_ = sigma;
    this->spot_ = spot;
    size_ = size;
    r_ = r;
    rho_ = rho;
    CorrelationMat = pnl_mat_create(size_,size_);
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            if (i == j) {
                pnl_mat_set(CorrelationMat,i,j,1);
            } else {
                pnl_mat_set(CorrelationMat,i,j,rho_);

            }
        }
    }
    // Factorisation de cholesky de la matrice de correlation
    cholesky = pnl_mat_chol(CorrelationMat);
}


void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng) {
    double pasTemps = T/(double) nbTimeSteps;    
    // Remplissage de la premiere ligne de la matrice des chemins
    // avec les prix spot
    for (int d = 0; d < size_; d++) {
        pnl_mat_set(path,0,d,pnl_vect_get(spot_,d));
    }

    // Creation d'une matrice D x (N+1) qui représente la suite de vecteurs gaussiens
    PnlMat *G = pnl_mat_create(size_,nbTimeSteps);
    for (int d = 0; d < size_; d++) {      
        for (int n = 0; n < nbTimeSteps; n++) {
            pnl_mat_set(G,d,n,pnl_rng_normal(rng));
        }
    }

    double LG;
    PnlVect *Ld = pnl_vect_create(size_);
    PnlVect *Gn = pnl_vect_create(size_);
    double exprExp;
    
    // Simulation du prix du sous-jacent à partir de t=0 
    for (int d = 0; d < size_; d++) {
        pnl_mat_get_row(Ld,CorrelationMat,d);       
        for (int n = 0; n < nbTimeSteps; n++) {
            pnl_mat_get_col(Gn,G,n);
            LG = pnl_vect_scalar_prod(Ld,Gn);
            exprExp = (r_ - (pow(pnl_vect_get(sigma_,d),2)/2.0)) * pasTemps + pnl_vect_get(sigma_,d) * sqrt(pasTemps) * LG;
            pnl_mat_set(path,n+1,d,(pnl_mat_get(path,n,d) * exp(exprExp)));
        }
    }
}



void BlackScholesModel::asset(PnlMat *path, double t, double T, int nbTimeSteps,
               PnlRng *rng, const PnlMat *past) {  
    // Remplissage d'une partie de la matrice des chemins
    // avec la matrice past contenant des donnees historiques
    pnl_mat_set_subblock(path,past,0,0);
    double u;
    double sTilde;
    double sSimul;
    double st;

    // Creation d'une matrice D x (N+1) qui représente la suite de vecteurs gaussiens
    PnlMat *G = pnl_mat_create(size_,nbTimeSteps);
    for (int d = 0; d < size_; d++) {      
        for (int n = 0; n < nbTimeSteps; n++) {
            pnl_mat_set(G,d,n,pnl_rng_normal(rng));
        }
    }

    double LG;
    PnlVect *Ld = pnl_vect_create(size_);
    PnlVect *Gn = pnl_vect_create(size_);
    double exprExp;

    bool prem = true;

    // Simulation du prix du sous-jacent à partir de la date t
    for (int d = 0; d < size_; d++) {
        st = pnl_mat_get(past,(past->m - 1),d);
        pnl_mat_get_row(Ld,CorrelationMat,d);
        for (int n = (past->m - 1); n < nbTimeSteps; n++) {
            pnl_mat_get_col(Gn,G,n);
            LG = pnl_vect_scalar_prod(Ld,Gn);

            if (prem) {
                u = (T/(double) nbTimeSteps) * (n + 1) - t;
            } else {
                u = (T/(double) nbTimeSteps);
            }
            sTilde = exp((r_ - (pow(pnl_vect_get(sigma_,d),2)/2.0)) * u + pnl_vect_get(sigma_,d) * sqrt(u) * LG);
            if (prem) {
                sSimul = st * sTilde;
            } else {
                sSimul = pnl_mat_get(path,n,d) * sTilde;
            }
            pnl_mat_set(path,n+1,d,sSimul);
            prem = false;
            //std::cout << "n+1: " << n+1 << std::endl;
       }
       prem = true;
    }

    pnl_mat_free(&G);
    pnl_vect_free(&Ld);
    pnl_vect_free(&Gn);
}

void BlackScholesModel::shiftAsset(PnlMat *shift_path, const PnlMat *path,
                    int d, double h, double t, double timestep) {
    double compteur = 0;
    double pas = 0;

    //Initialisation de la matrice shift_path avec les valeurs de path
    pnl_mat_set_subblock(shift_path,path,0,0);

    //Recuperation de la ligne de la matrice path qui correspond a la date t
    while (pas <= t) {
        compteur++;
        pas += timestep;
    }

    double valeurShift = 0;
    
    //On shift sur le sous-jacent d
    for (int i = compteur; i < path->m ; i++) {
        valeurShift = pnl_mat_get(path,compteur,d) * (1 + h);
        pnl_mat_set(shift_path,i,d,valeurShift);
    }
}