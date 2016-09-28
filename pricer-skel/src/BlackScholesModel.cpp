#include "BlackScholesModel.hpp"
#include "Market.hpp"
#include <iostream>

BlackScholesModel::~BlackScholesModel() {
}


BlackScholesModel::BlackScholesModel() {
    this->sigma_ = pnl_vect_create(0);
    this->spot_ = pnl_vect_create(0);
    size_ = 0;
    r_ = 0;
    rho_ = 0;
    CorrelationMat = pnl_mat_create(0,0);
    cholesky = 0;
}

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot, PnlVect *mu, int hedgingDateNb) {
    this->sigma_ = sigma;
    this->spot_ = spot;
    size_ = size;
    r_ = r;
    rho_ = rho;
    this->mu_ = mu;
    hedgingDateNb_ = hedgingDateNb;
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
    pnl_mat_free(&G);
    pnl_vect_free(&Ld);
    pnl_vect_free(&Gn);
}



void BlackScholesModel::asset(PnlMat *path, double t, double T, int nbTimeSteps,
               PnlRng *rng, const PnlMat *past) {  
    double pasTemps = T/(double) nbTimeSteps;
    int debutSimulation;
    double diff = t;
    // On regarde si t est une date de rebalancement

    while (diff > 0 ) {
        diff -= pasTemps;
    }
    if(diff == 0) {
        debutSimulation =past->m - 1;
    } else {
        debutSimulation = past->m - 2;
    }

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
    double volatilite;

    // Simulation du prix du sous-jacent à partir de la date t
    for (int d = 0; d < size_; d++) {
        st = pnl_mat_get(past,(past->m - 1),d);
        volatilite = pnl_vect_get(sigma_,d);
        pnl_mat_get_row(Ld,CorrelationMat,d);
        for (int n = debutSimulation; n < nbTimeSteps; n++) {
            pnl_mat_get_col(Gn,G,n);
            LG = pnl_vect_scalar_prod(Ld,Gn);

            if (prem) {
                u = (T/(double) nbTimeSteps) * (n + 1) - t;
            } else {
                u = (T/(double) nbTimeSteps);
            }
            sTilde = exp((r_ - (pow(volatilite,2)/2.0)) * u + volatilite * sqrt(u) * LG);

            if (prem) {
                sSimul = st * sTilde;
            } else {
                sSimul = pnl_mat_get(path,n,d) * sTilde;
            }
            
            pnl_mat_set(path,n+1,d,sSimul);
            prem = false;
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
        valeurShift = pnl_mat_get(path,i,d) * (1 + h);
        pnl_mat_set(shift_path,i,d,valeurShift);
    }
}

/*
void simul_market(int nbAssets, PnlMat *market)
{
    char *marketFile = "market.dat";
    PnlMat *market_from_file = pnl_mat_create_from_file(marketFile);
    pnl_mat_extract_subblock(market, market_from_file, 0, market_from_file->m, 0, nbAssets);
    pnl_mat_free(&market_from_file);
}
*/

void BlackScholesModel::simul_market(PnlMat *path, double T, PnlRng *rng) {
    double pasTemps = T/(double) hedgingDateNb_;    
    
    // Remplissage de la premiere ligne de la matrice des chemins
    // avec les prix spot
    for (int d = 0; d < size_; d++) {
        pnl_mat_set(path,0,d,pnl_vect_get(spot_,d));
    }

    // Creation d'une matrice D x (N+1) qui représente la suite de vecteurs gaussiens
    PnlMat *G = pnl_mat_create(size_,hedgingDateNb_);
    for (int d = 0; d < size_; d++) {      
        for (int n = 0; n < hedgingDateNb_; n++) {
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
        for (int n = 0; n < hedgingDateNb_; n++) {
            pnl_mat_get_col(Gn,G,n);
            LG = pnl_vect_scalar_prod(Ld,Gn);
            exprExp = (pnl_vect_get(mu_,d) - (pow(pnl_vect_get(sigma_,d),2)/2.0)) * pasTemps + pnl_vect_get(sigma_,d) * sqrt(pasTemps) * LG;
            pnl_mat_set(path,n+1,d,(pnl_mat_get(path,n,d) * exp(exprExp)));
        }
    }
    pnl_mat_free(&G);
    pnl_vect_free(&Ld);
    pnl_vect_free(&Gn);
}

