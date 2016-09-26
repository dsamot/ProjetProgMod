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


MonteCarlo::~MonteCarlo() {
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
    pnl_mat_free(&path);
}



void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& ic) {

    if (t < 0) {
        std::cerr << "La date de pricing est négative!" << std::endl;
    } else if (t > opt_->T_) {
        std::cerr << "La date de pricing est supérieure à la maturité de l'option!" << std::endl;
    } else {
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
        pnl_mat_free(&path);
    }
}


void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta) {
    if (t < 0) {
        std::cerr << "La date de pricing est négative!" << std::endl;
    } else if (t > opt_->T_) {
        std::cerr << "La date de pricing est supérieure à la maturité de l'option!" << std::endl;
    } else {
        double M = (double)nbSamples_;
        double interet = mod_->r_ ;
        double maturite = opt_->T_;
        double tempdelta = 0;
        double facteurExp;
        double expo = exp(-interet*(maturite - t));

        PnlMat *shift_path_up = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 
        PnlMat *shift_path_down = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 

        PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_);
        //mod_->asset(path, t, maturite, opt_->nbTimeSteps_, rng_, past);
        double sommeDiffPayOff = 0;
        //pnl_mat_print(path);

        /*for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
            
            for(int i=0; i < M; i++) {
                //std::cout << "i: " << i << std::endl;
                mod_->asset(path, t, maturite, opt_->nbTimeSteps_, rng_, past);
                mod_->shiftAsset(shift_path_up, path, idAsset, fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                mod_->shiftAsset(shift_path_down, path, idAsset, -fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));


                sommeDiffPayOff += opt_->payoff(shift_path_up) - opt_->payoff(shift_path_down);
            }
            //std::cout << "sommeDiffPayOff: " << sommeDiffPayOff << std::endl;
            facteurExp = expo/(M*2*pnl_mat_get(past, (past->m - 1), idAsset)*fdStep_);
            //facteurExp = exp(-interet*(maturite - t))/(M*2*pnl_mat_get(past, (past->m - 1), idAsset)*fdStep_);
            tempdelta = facteurExp*sommeDiffPayOff;
            std::cout << "d: " << tempdelta << " indice: " << idAsset << std::endl;
            pnl_vect_set(delta,idAsset, tempdelta);
            sommeDiffPayOff = 0;
        }*/

        /*for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
            facteurExp = expo/(M*2*pnl_mat_get(past, (past->m - 1), idAsset)*fdStep_);
        }*/

        for(int i=0; i < M; i++) {
            mod_->asset(path, t, maturite, opt_->nbTimeSteps_, rng_, past);
            for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
                mod_->shiftAsset(shift_path_up, path, idAsset, fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                mod_->shiftAsset(shift_path_down, path, idAsset, -fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                sommeDiffPayOff = opt_->payoff(shift_path_up) - opt_->payoff(shift_path_down);
                pnl_vect_set(delta,idAsset,(pnl_vect_get(delta,idAsset) + sommeDiffPayOff));
            }
        }

        for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
            facteurExp = expo/(M*2*pnl_mat_get(past, (past->m - 1), idAsset)*fdStep_);
            pnl_vect_set(delta,idAsset,(pnl_vect_get(delta,idAsset)*facteurExp));
        }

        pnl_vect_print(delta);
        pnl_mat_free(&shift_path_up);
        pnl_mat_free(&shift_path_down);
        pnl_mat_free(&path);
    }
}


