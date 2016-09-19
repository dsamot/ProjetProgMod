#include <stdio.h>
#include <math.h>
#include "MonteCarlo.hpp"
#include "parser.hpp"
#include "BasketOption.h"



    MonteCarlo::MonteCarlo() {
        this->mod_= new BlackScholesModel(); 
        //this->opt_ = new BasketOption(); 
        //this->rng_ = new PnlRng();
        this->fdStep_ =0.0; 
        this->nbSamples_ = 0.0;
    }

    MonteCarlo::MonteCarlo(BlackScholesModel *bsModel, Option* option,PnlRng* pnlRng,double steps,int nbSamples) {
        this->mod_= bsModel; 
        this->opt_ = option; 
        this->rng_ = pnlRng;
        this->fdStep_ = steps; 
        this->nbSamples_ = nbSamples;
    }


void MonteCarlo::price(double &prix, double &ic){
    double interet = mod_->r_ ;
    double maturite = opt_->T_;
    
    double termeExp = exp(-interet*maturite);
    double sommePayOff =0;
    double sommePayOffCarre = 0;
    int M = nbSamples_;
    for (int i=0; i< M; i++){
        PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_,mod_->size_);
        mod_->asset(path,maturite,opt_->nbTimeSteps_,rng_);
        double res =  opt_->payoff(path);
        std::cout << "res" << res << std::endl;
        sommePayOff += res;
        sommePayOffCarre += res*res;
    }
    double moyenne = sommePayOff/M;
    double variance = exp(-2*interet*maturite)*(sommePayOffCarre/M - pow(moyenne,2.0));
    ic = sqrt(variance/M)*2.0*1.96;
    prix = termeExp*moyenne;
}



void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& ic) {
    double interet = mod_->r_ ;
    double maturite = opt_->T_;
    
    double termeExp = exp(-interet*(maturite-t));
    double sommePayOff =0;
    double sommePayOffCarre = 0;
    int M = nbSamples_;
    for (int i=0; i< M; i++){
        PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_,mod_->size_);
        pnl_mat_set_subblock(path,past,(opt_->nbTimeSteps_ - past->m),(mod_->size_ - past->n));
        mod_->asset(path,maturite,opt_->nbTimeSteps_,rng_);
        double res =  opt_->payoff(path);
        sommePayOff += res;
        sommePayOffCarre += res*res;
    }
    double moyenne = (double)sommePayOff/(double)M;
    double variance = exp(-2*interet*(maturite-t))*(sommePayOffCarre/M - pow(moyenne,2.0));
    ic = sqrt(variance/M)*2.0*1.96;
    prix = termeExp*moyenne;
}

