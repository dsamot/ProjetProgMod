#include <iostream>
#include <string>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"


#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"
#include "BasketOption.hpp"
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
    int size, timeStepsNb, hedgingDateNumber;
    size_t sample;

    char *infile = argv[1];
    Param *P = new Parser(infile);

    // Extraction des parametres avec le parser
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("maturity", maturity);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", interest);
    P->extract("correlation", corr);
    P->extract("trend", mu, size);
    P->extract("strike", strike);
    P->extract("option type", type);
    P->extract("payoff coefficients", payoffCoeff, size);
    P->extract("timestep number", timeStepsNb);
    P->extract("hedging dates number", hedgingDateNumber);
    P->extract("fd step", steps);
    P->extract("sample number", sample);
    if (P->extract("dividend rate", divid, size) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }

    //Cas ou on connait deja des donnees historiques
    PnlMat *past = pnl_mat_create(1,size);
    for(int i = 0; i < size; i++) {
        pnl_mat_set(past,0,i,100);
    }

   
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

    hedgingDateNumber = 5;
    mu = pnl_vect_create_from_scalar(size,0.2);
    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot,mu,hedgingDateNumber);
    MonteCarlo *montecarlo = new MonteCarlo(model,option,rng,steps,sample);
    
    double prix;
    double ic;
    montecarlo->price(prix,ic);

    std::cout << "prix : " << prix << std::endl;
    std::cout << "ic : " << ic << std::endl;

    
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