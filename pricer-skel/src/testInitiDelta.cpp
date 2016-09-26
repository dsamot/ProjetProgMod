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
#include "parser.hpp"
#include "Market.hpp"
using namespace std;

int main(int argc, char **argv)
{
    double steps = 0.3;
    double maturity, interest, strike, corr;
    PnlVect *spot, *mu, *sigma, *divid, *payoffCoeff;
    string type;
    int size, timeStepsNb;
    size_t sample;

    char *infile = argv[1];
    Param *P = new Parser(infile);

    P->extract("option type", type);
    P->extract("maturity", maturity);
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("mu", mu, size);
    P->extract("volatility", sigma, size);
    P->extract("correlation", corr);
    P->extract("payoff coefficients", payoffCoeff, size);

    P->extract("interest rate", interest);
    if (P->extract("dividend rate", divid, size) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }
    P->extract("strike", strike);
    P->extract("sample number", sample);
    P->extract("timestep number", timeStepsNb);

    //Cas ou on connait deja des donnees historiques
    PnlMat *past = pnl_mat_create(1,size);
    for(int i = 0; i < size; i++) {
        pnl_mat_set(past,0,i,100);
    }

    PnlMat *path = pnl_mat_create(timeStepsNb + 1,size);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    Option *option;

    if(type == "basket") {
        option = new BasketOption(maturity, timeStepsNb, size, strike, payoffCoeff);
    } else if(type == "asian") {
        option = new AsianOption(maturity, timeStepsNb, size, strike, payoffCoeff);
    } else if(type == "performance") {
        option = new PerformanceOption(maturity, timeStepsNb, size, payoffCoeff);
    } else {
        exit(0);
    }


    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot);
    MonteCarlo *montecarlo = new MonteCarlo(model,option,rng,steps,sample);

    double prix;
    double ic;
    PnlVect *delta = pnl_vect_create(size);
    
    montecarlo->delta(past,0,delta);

    
    pnl_vect_free(&spot);
    pnl_vect_free(&sigma);
    pnl_vect_free(&divid);
    pnl_mat_free(&past);
    pnl_rng_free(&rng);
    delete P;
    delete model;
    delete option;
    delete montecarlo;

    exit(0);
}