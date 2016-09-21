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
    double M = (double)nbSamples_;
    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1,mod_->size_);
    for (int i=0; i< M; i++){
        mod_->asset(path,maturite,opt_->nbTimeSteps_,rng_);
        double res =  opt_->payoff(path);
        //std::cout << "res" << res << std::endl;
        sommePayOff += res;
        sommePayOffCarre += res*res;
    }
    double moyenne = sommePayOff/M;
    double variance = (sommePayOffCarre/M - pow(moyenne,2.0));
    //double variance = exp(-interet*maturite)*(sommePayOffCarre/M - pow(moyenne,2.0));
    ic = sqrt(variance/M)*2.0*1.96*termeExp;
    std::cout << "Variance " <<  variance*exp(-2*interet*maturite) << std::endl;
    std::cout << "Standard Variation " <<  sqrt(variance*exp(-2*interet*maturite)/M) << std::endl;
    prix = termeExp*moyenne;
}



void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& ic) {
    double interet = mod_->r_ ;
    double maturite = opt_->T_;
    
    double termeExp = exp(-interet*(maturite-t));
    double sommePayOff =0;
    double sommePayOffCarre = 0;
    double M = (double)nbSamples_;
    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1,mod_->size_);
    for (int i=0; i< M; i++){
        mod_->asset(path,t,maturite,opt_->nbTimeSteps_,rng_, past);
        double res =  opt_->payoff(path);
        sommePayOff += res;
        sommePayOffCarre += res*res;
    }
    double moyenne = (double)sommePayOff/(double)M;
    double variance = exp(-2*interet*(maturite-t))*(sommePayOffCarre/M - pow(moyenne,2.0));
    ic = sqrt(variance/M)*2.0*1.96;
    prix = termeExp*moyenne;
}


void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta) {
    double M = (double)nbSamples_;
    double interet = mod_->r_ ;
    double maturite = opt_->T_;


    PnlMat *shift_path_up = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 
    PnlMat *shift_path_down = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 

    double sommeDiffPayOff = 0;

    for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
        for( int i=1; i<M; i++) {
            mod_->shiftAsset(shift_path_up, past, idAsset, fdStep_, t , opt_->nbTimeSteps_);
            mod_->shiftAsset(shift_path_down, past, idAsset, -fdStep_, t , opt_->nbTimeSteps_);

            sommeDiffPayOff += opt_->payoff(shift_path_up) - opt_->payoff(shift_path_down);
        }
        double facteurExp = exp(-interet*(maturite - t))/(M*2*pnl_mat_get(past, past->m, idAsset)*fdStep_);
        pnl_vect_set(delta,idAsset, facteurExp*sommeDiffPayOff);
        sommeDiffPayOff = 0;
    }
    pnl_vect_print(delta);
}


