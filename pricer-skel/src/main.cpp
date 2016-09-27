#include <iostream>
#include <string>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include <time.h>
#include <string.h>


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
    clock_t tStart = clock();
    

    double steps = 0.1;
    double maturity, interest, strike, corr;
    PnlVect *spot, *mu, *sigma, *divid, *payoffCoeff;
    string type;
    int size, timeStepsNb, hedgingDateNumber;
    size_t sample;


    int argFileNb;

    if(argc == 2) {
        argFileNb = 1;
    } else if(argc == 3) {
        // if -c // 
        char *parametre = argv[1];
        if(strcmp(parametre,"-c") != 0){
            std::cerr << "Le paramètre passé n'est pas le bon"<< parametre << std::endl;
            return 1;
        }
        argFileNb = 2;

    } else {
        std::cerr << "Le nombre de paramètres passés au programme n'est pas le bon" << std::endl;
        return 1;
    }

    char *infile = argv[argFileNb];


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


    //hedgingDateNumber = 5;
    mu = pnl_vect_create_from_scalar(size,0.2);
    BlackScholesModel *model = new BlackScholesModel(size,interest,corr,sigma,spot,mu,hedgingDateNumber);
    MonteCarlo *montecarlo = new MonteCarlo(model,option,rng,steps,sample);

    double prix;
    double ic;
    double p0 ;
    double PnL;



    if(argc == 2) {
        //    std::cout << << std::endl;
        montecarlo->price(prix,ic);
        PnlVect *delta = pnl_vect_create(size);
        montecarlo->delta(past,0,delta);
        std::cout << "---------------- Prix et delta à t=0 -----------------" << std::endl;
        std::cout << "Prix à t=0: " << prix << " - Intervalle de confiance : " << ic << std::endl;
        std::cout << "delta à t=0:"<< std::endl;
        pnl_vect_print(delta);
        printf("\nTemps d'execution: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        return 1;
    }


    PnlVect *V = pnl_vect_create(hedgingDateNumber + 1);
    montecarlo->profitAndLoss(V,PnL);
    pnl_vect_print(V);
    std::cout << " le P&L du portefeuille de couverture est : " << PnL << std::endl;
    
    printf("\nTemps d'execution: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    pnl_vect_free(&V);
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