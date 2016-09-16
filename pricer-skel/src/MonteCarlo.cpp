#include <stdio.h>
#include <math.h>
#include "MonteCarlo.hpp"
#include "parser.hpp"




    MonteCarlo::MonteCarlo() {
        this->mod_= new BlackScholesModel(); 
        this->opt_ = new Option(); 
        this->rng_ = new PnlRng();
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
    double interet = this.mod_->r_ ;
    double maturite = this.opt_->T_;
    
    double termeExp = exp(-interet*maturite);
    double sommePayOff =0;
    double sommePayOffCarre = 0;
    int M = this.nbSamples_;
    for (int i=0; i< M; i++){
        PnlMat *path = new PnlMat();
        this.mod_->asset(path,maturite,this.fdStep_,this.rng_);
        double res =  this.opt_->payoff(path);
        sommePayOff += res;
        sommePayOffCarre += res*res;
    }
    double moyenne = sommePayOff/M;
    double variance = exp(-2*interet*maturite)*(sommePayOffCarre/M - pow(moyenne,2.0));
    *ic = sqrt(variance/M)*2.0*1.96/;
    *prix = termeExp*moyenne;
}



void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& ic) {
    double interet = this.mod_->r_ ;
    double maturite = this.opt_->T_;
    
    double termeExp = exp(-interet*(maturite-t));
    double sommePayOff =0;
    double sommePayOffCarre = 0;
    int M = this.nbSamples_;
    for (int i=0; i< M; i++){
        PnlMat *path = new PnlMat(past);
        this.mod_->asset(path,maturite,this.fdStep_,this.rng_);
        double res =  this.opt_->payoff(path);
        sommePayOff += res;
        sommePayOffCarre += res*res;
    }
    double moyenne = sommePayOff/M;
    double variance = exp(-2*interet*(maturite-t))*(sommePayOffCarre/M - pow(moyenne,2.0));
    *ic = sqrt(variance/M)*2.0*1.96/;
    *prix = termeExp*moyenne;
}

