#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"


class Market
{
public:

    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales du sous-jacent
     PnlVect *mu_; 
    PnlMat *CorrelationMat;
    int cholesky;
    double maturity_;
    int nbTimeSteps_;
    int size_;
    double r_;

    /*!
     * \brief Constructeur par défaut
     */
    Market();
    
    /*!
     * \brief Constructeur 
     */
    Market(PnlVect *sigma, PnlVect *spot, PnlVect *mu,  double rho, double maturity, int nbTimeSteps, int size, double r_);
    


};


