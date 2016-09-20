#include <iostream>
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
    int timeStepsNb = 12;
    double interest = 0.03;
    double corr = 0.5;
    int size = 5;
    int sample = 50000;
    double strike = 100.0;
    double steps = 0;
    double spotprice = 100;
    double vol = 0.2;
    PnlMat *path = pnl_mat_create(timeStepsNb + 1,size);
    PnlVect *sigma = pnl_vect_create(size);
    PnlVect *spot = pnl_vect_create(size);
    for(int i = 0; i < size; i++) {
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
  //  BasketOption *basket = new BasketOption(maturity,timeStepsNb,size,strike);
    //AsianOption *asian = new AsianOption(maturity, timeStepsNb, size, strike);
   PerformanceOption *performance = new PerformanceOption(maturity, timeStepsNb, size);
    MonteCarlo *montecarlo = new MonteCarlo(model,performance,rng,steps,sample);
    double prix;
    double ic;
    montecarlo->price(prix,ic);
    std::cout << "prix : " << prix << std::endl;
    std::cout << "ic : " << ic << std::endl;
    //model->asset(path, T, nbTimeSteps, rng);
    return 0;
}
