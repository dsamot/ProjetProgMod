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
        double sommeDiffPayOff = 0;

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




void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *ic) {
    /*if (t < 0) {
        std::cerr << "La date de pricing est négative!" << std::endl;
    } else if (t > opt_->T_) {
        std::cerr << "La date de pricing est supérieure à la maturité de l'option!" << std::endl;
    } else {
        double M = (double)nbSamples_;
        double interet = mod_->r_ ;
        double maturite = opt_->T_;
        double tempdelta = 0;
        double facteur;

        PnlMat *shift_path_up = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 
        PnlMat *shift_path_down = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 

        PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_);
        double sommeDiffPayOff = 0;
        double sommeDiffPayOffSquared = 0;
        double diff = 0;
        for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
            for(int i=0; i < M; i++) {
                mod_->asset(path, t, maturite, opt_->nbTimeSteps_, rng_, past);
                mod_->shiftAsset(shift_path_up, path, idAsset, fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                mod_->shiftAsset(shift_path_down, path, idAsset, -fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                diff = opt_->payoff(shift_path_up) - opt_->payoff(shift_path_down);
                sommeDiffPayOff += diff;
                sommeDiffPayOffSquared += diff*diff;
            }
            //std::cout << "fdStep: " << fdStep_<< std::endl; 
            facteur = 1/(2*pnl_mat_get(past, (past->m - 1), idAsset)*fdStep_);
            tempdelta = exp(-interet*(maturite - t))/M*facteur*sommeDiffPayOff;
            pnl_vect_set(delta,idAsset, tempdelta);

            double variance = exp(-2*interet*(maturite-t))*facteur*facteur*(sommeDiffPayOffSquared/M - pow(sommeDiffPayOff/M,2));
            pnl_vect_set(ic,idAsset,sqrt(variance)*1.96/sqrt(M));
            sommeDiffPayOff = 0;
            sommeDiffPayOffSquared = 0;
        }
        pnl_mat_free(&shift_path_up);
        pnl_mat_free(&shift_path_down);
        pnl_mat_free(&path);
    }*/
    if (t < 0) {
        std::cerr << "La date de pricing est négative!" << std::endl;
    } else if (t > opt_->T_) {
        std::cerr << "La date de pricing est supérieure à la maturité de l'option!" << std::endl;
    } else {
        double M = (double)nbSamples_;
        double interet = mod_->r_ ;
        double maturite = opt_->T_;
        double tempdelta = 0;
        double factor;
        double expo = exp(-interet*(maturite - t));

        PnlMat *shift_path_up = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 
        PnlMat *shift_path_down = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_); 

        PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ +1, mod_->size_);
        double sommeDiffPayOff = 0;
        PnlVect sommeDiffPayOffSquared = pnl_vect_create(mod_->size_);

        for(int i=0; i < M; i++) {
            mod_->asset(path, t, maturite, opt_->nbTimeSteps_, rng_, past);
            for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
                mod_->shiftAsset(shift_path_up, path, idAsset, fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                mod_->shiftAsset(shift_path_down, path, idAsset, -fdStep_, t , (maturite / (double)opt_->nbTimeSteps_));
                sommeDiffPayOff = opt_->payoff(shift_path_up) - opt_->payoff(shift_path_down);

                pnl_vect_set(delta,idAsset,(pnl_vect_get(delta,idAsset) + sommeDiffPayOff));
                pnl_vect_set(sommeDiffPayOffSquared,idAsset,(pnl_vect_get(sommeDiffPayOffSquared,idAsset) + sommeDiffPayOff*sommeDiffPayOff));
            }
        }

        for (int idAsset=0; idAsset < mod_->size_; idAsset++) {
            factor = 1/(2*pnl_mat_get(past, (past->m - 1), idAsset)*fdStep_);
            double carreUn = pnl_vect_get(sommeDiffPayOffSquared,idAsset);
            double carreDeux = pnl_vect_get(delta,idAsset);
            pnl_vect_set(ic,idAsset,expo*factor*sqrt((carreUn/M - pow(carreDeux/M,2))/M)*1.96);
            pnl_vect_set(delta,idAsset,(pnl_vect_get(delta,idAsset)*expo*factor/M));
        }

        pnl_vect_print(delta);
        pnl_vect_print(ic);
        pnl_mat_free(&shift_path_up);
        pnl_mat_free(&shift_path_down);
        pnl_mat_free(&path);
    }
}
