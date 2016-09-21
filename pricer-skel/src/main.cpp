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

using namespace std;

int main(int argc, char **argv)
{
    double steps = 0;
    double maturity, interest, strike, corr;
    PnlVect *spot, *sigma, *divid;
    string type;
    int size, timeStepsNb;
    size_t sample;

    char *infile = argv[1];
    Param *P = new Parser(infile);

    P->extract("option type", type);
    P->extract("maturity", maturity);
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("correlation", corr);

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
        option = new BasketOption(maturity, timeStepsNb, size, strike);
    } else if(type == "asian") {
        option = new AsianOption(maturity, timeStepsNb, size, strike);
    } else if(type == "performance") {
        option = new PerformanceOption(maturity, timeStepsNb, size);
    } else {
        exit(0);
    }


    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot);
    MonteCarlo *montecarlo = new MonteCarlo(model,option,rng,steps,sample);
    //Cas data historiques
    //std::cout << "nbTimeSteps : " << (timeStepsNb *(maturity-maturity) / maturity) << std::endl;
    //model->asset(past, 1, 0, rng);

    double prix;
    double ic;
    montecarlo->price(prix,ic);

    std::cout << "prix : " << prix << std::endl;
    std::cout << "ic : " << ic << std::endl;
    
    //Cas data historiques
    montecarlo->price(past,0,prix,ic);

    std::cout << "prix histo : " << prix << std::endl;
    std::cout << "ic histo : " << ic << std::endl;

    pnl_vect_free(&spot);
    pnl_vect_free(&sigma);
    pnl_vect_free(&divid);
    delete P;

    exit(0);
}