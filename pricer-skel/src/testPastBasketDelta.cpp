#include <iostream>
#include <string>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"

#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"
#include "BasketOption.h"
#include "AsianOption.hpp"
#include "PerformanceOption.hpp"
using namespace std;

int main()
{
    
    double maturity = 2;
    int timeStepsNb = maturity*12;
    double interest = 0.04879;
    double corr = -1;
    int size = 2;
    int sample = 50000;
    double strike = 100.0;
    double steps = 0.1;
    double spotprice = 100;
    double vol = 0.2;

    PnlMat *past = pnl_mat_create(12 +1,size);
    PnlVect *sigma = pnl_vect_create(size);
    PnlVect *spot = pnl_vect_create(size);
    for (int i=0; i<size; i++) {
        pnl_vect_set(spot,i,spotprice);
        pnl_vect_set(sigma,i,vol);
    }
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot);
    
    // Creation des donnees historiques
    model->asset(past, maturity, 12, rng);
    BasketOption *basket = new BasketOption(maturity,timeStepsNb,size,strike);
    MonteCarlo *montecarlo = new MonteCarlo(model,basket,rng,steps,sample);

    PnlVect *delta = pnl_vect_create(size); 
    montecarlo->delta(past,maturity*0.5,delta);
    
    pnl_vect_free(&spot);
    pnl_vect_free(&sigma);
    pnl_mat_free(&past);
    pnl_rng_free(&rng);

    delete model;
    delete basket;
    delete montecarlo;
    return 0;
}
