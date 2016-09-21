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
    /*
    // Data pour basket option
    double maturity = 3;
    int timeStepsNb = 1;
    double interest = 0.04879;
    double corr = 0.0;
    int size = 40;*/
    
    // Data pour performance option
    double maturity = 365;
    int timeStepsNb = 365/2;
    double interest = 0.04879;
    double corr = 0;
    int size = 2;
    int sample = 50000;
    double strike = 20.0;
    double steps = 0;
    double spotprice = 20;
    double vol = 0.2;
    PnlMat *path = pnl_mat_create(365/2 + 1,size);
    PnlMat *past = pnl_mat_create(365/2+1,size);
    PnlVect *sigma = pnl_vect_create(size);
    PnlVect *spot = pnl_vect_create(size);
    for (int i=0; i<size; i++) {
        pnl_vect_set(spot,i,spotprice);
        pnl_vect_set(sigma,i,vol);
    }
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    /*int M = 1E5;
    int dim = 2;
    pnl_rng_sseed(rng, time(NULL));

    double acc = 0., var = 0;

    for (int i = 0; i < M; i++)
    {
        pnl_vect_rng_normal(G, dim, rng);
        double tmp = pnl_vect_norm_two(G);
        acc += tmp;
        var += tmp * tmp;
    }

    acc /= M;
    var = var / M - acc * acc;

    cout << "E[||G||_2] = " << acc << endl;
    cout << "Var(||G||_2) = " << var << endl;

    pnl_vect_free(&G);
    pnl_rng_free(&rng);*/

    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot);
    
    // Creation des donnees historiques
    model->asset(past, maturity/2, maturity/2, rng);
    BasketOption *basket = new BasketOption(maturity,timeStepsNb,size,strike);
    MonteCarlo *montecarlo = new MonteCarlo(model,basket,rng,steps,sample);

    //BasketOption *basket = new BasketOption(maturity,timeStepsNb,size,strike);
    //AsianOption *asian = new AsianOption(maturity, timeStepsNb, size, strike);
    //PerformanceOption *performance = new PerformanceOption(maturity, timeStepsNb, size);
    //MonteCarlo *montecarlo = new MonteCarlo(model,performance,rng,steps,sample);
    double prix;
    double ic;
    //montecarlo->price(past,1,prix,ic);
    //pnl_mat_print(past);
    montecarlo->price(past,365/2,prix,ic);

    std::cout << "prix : " << prix << std::endl;
    std::cout << "ic : " << ic << std::endl;
    cout << "Prix S(Maturité)" << pnl_mat_get(past, 0, 365) << endl;
    return 0;
}
