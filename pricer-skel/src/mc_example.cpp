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
    
    // Data pour performance option
    double maturity = 2;
    int timeStepsNb = maturity*12;
    double interest = 0.04879;
    double corr = 0;
    int size = 2;
    int sample = 50000;
    double strike = 100.0;
    double steps = 0.1;
    double spotprice = 100;
    double vol = 0.2;

    PnlMat *past = pnl_mat_create(12 +1,size);
    PnlVect *sigma = pnl_vect_create(size);
    PnlVect *spot = pnl_vect_create(size);
    PnlVect *payoffCoeff = pnl_vect_create(size);
    for (int i=0; i<size; i++) {
        pnl_vect_set(spot,i,spotprice);
        pnl_vect_set(sigma,i,vol);
        pnl_vect_set(payoffCoeff,i,(1/(double)size));
    }
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot);
    
    // Creation des donnees historiques
    model->asset(past, maturity, 12, rng);
    BasketOption *basket = new BasketOption(maturity,timeStepsNb,size,strike,payoffCoeff);
    MonteCarlo *montecarlo = new MonteCarlo(model,basket,rng,steps,sample);
    /*double prix;
    double ic;
    // T1 - t
    montecarlo->price(past,365/(maturity*365) * maturity ,prix,ic);
    std::cout << "prix1 : " << prix << std::endl;
    


    maturity = 1;
    timeStepsNb = 365*maturity;
    for (int i=0; i<size; i++) {
        pnl_vect_set(spot,i,spotprice);
        pnl_vect_set(sigma,i,vol);
    }
    
    // Creation des donnees historiques
    basket = new BasketOption(maturity,timeStepsNb,size,strike);
    montecarlo = new MonteCarlo(model,basket,rng,steps,sample);

    montecarlo->price(prix,ic);
    std::cout << "prix2 : " << prix << std::endl;
    */
    PnlVect *delta = pnl_vect_create(size); 
    montecarlo->delta(past,maturity*0.5,delta);
    //montecarlo->delta(past,7,delta);
    
    pnl_vect_free(&spot);
    pnl_vect_free(&sigma);
    pnl_mat_free(&past);
    pnl_rng_free(&rng);

    delete model;
    delete basket;
    delete montecarlo;
    return 0;
}
