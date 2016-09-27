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
    double steps = 0.1;
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

    // Gestion des cas interdits
    if (maturity < 0 ) {
        std::cerr << "La date de maturité est négative" << std::endl;
        return 1;
    } else if (size <= 0) {
        std::cerr << "La taille de l'option est négative ou nulle" << std::endl;
        return 1;
    } else if (strike < 0) {
        std::cerr << "Le strike de l'option a une valeur négative" << std::endl;
        return 1;
    } else if (sample <= 0) {
        std::cerr << "Le nombre d'échantillon pour la méthode de montecarlo est négatif ou nul" << std::endl;
        return 1;
    } else if (timeStepsNb <= 0) {
        std::cerr << "Le nombre de dates à laquelle on évalue le prix du sous-jacent est négatif ou nul" << std::endl;
        return 1;
    } else if (hedgingDateNumber <= 0) {
        std::cerr << "Le nombre de dates à laquelle on évalue le prix du portefeuille de couverture est négatif ou nul" << std::endl;
        return 1;
    }
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


    hedgingDateNumber = 7;
    mu = pnl_vect_create_from_scalar(size,0.2);
    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot,mu,hedgingDateNumber);
    MonteCarlo *montecarlo = new MonteCarlo(model,option,rng,steps,sample);
    //Cas data historiques
    //std::cout << "nbTimeSteps : " << (timeStepsNb *(maturity-maturity) / maturity) << std::endl;
    //model->asset(past, 1, 0, rng);


/*
    PnlMat* result;
    std::cout << "Debut" << std::endl;
    Market *myMarket = new Market(sigma, spot, mu,  corr, maturity, timeStepsNb, size, interest);
    std::cout << "Milieu" << std::endl;
    result = model->simul_market(*myMarket, rng);
    std::cout << "Fin" << std::endl;
    pnl_mat_print(result);*/


    double prix;
    double ic;
    double p0 ;
    double PnL;
    //montecarlo->price(p0, ic);
    PnlVect *V = pnl_vect_create(hedgingDateNumber + 1);
    montecarlo->profitAndLoss(V,PnL);
    //pnl_vect_print(V);
    std::cout << "PnL : " << PnL << std::endl;

    //model->profitLoss(*myMarket, *result, p0, montecarlo);


    //PnlVect *delta = pnl_vect_create(size);
    /*montecarlo->price(prix,ic);

    std::cout << "prix : " << prix << std::endl;
    std::cout << "ic : " << ic << std::endl;*/
    
    //Cas data historiques
    /*montecarlo->price(past,2*maturity/timeStepsNb,prix,ic);
    
    std::cout << "prix histo : " << prix << std::endl;
    std::cout << "ic histo : " << ic << std::endl;*/

   // montecarlo->delta(past,0,delta);

    
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