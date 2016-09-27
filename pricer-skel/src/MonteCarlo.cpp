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

        //pnl_vect_print(delta);
        pnl_mat_free(&shift_path_up);
        pnl_mat_free(&shift_path_down);
        pnl_mat_free(&path);
    }
}




void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *ic) {
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
        PnlVect *sommeDiffPayOffSquared = pnl_vect_create(mod_->size_);

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

        //pnl_vect_print(delta);
        //pnl_vect_print(ic);
        pnl_mat_free(&shift_path_up);
        pnl_mat_free(&shift_path_down);
        pnl_mat_free(&path);
    }
}

void MonteCarlo::profitAndLoss(PnlVect *V, double &PnL) {
    PnlMat *marketPath = pnl_mat_create(mod_->hedgingDateNb_ + 1, opt_->size_);
    mod_->simul_market(marketPath, opt_->T_, rng_);
    //pnl_mat_print(marketPath);
    // Initialisation du portefeuille pour t=0
    double valeur;
    double prix;
    double ic;
    PnlVect *delta1 = pnl_vect_create(opt_->size_);
    PnlVect *delta2 = pnl_vect_create(opt_->size_);
    PnlVect *Stauxi = pnl_vect_create(opt_->size_);
    //PnlMat *past = pnl_mat_create(1,opt_->size_);
    PnlMat *past = pnl_mat_new();
    PnlMat *pastcopy = pnl_mat_new();
    pnl_mat_extract_subblock(past, marketPath, 0, 1, 0, marketPath->n);
    pnl_mat_get_row(Stauxi,marketPath,0);
    price(prix,ic);
    delta(past,0,delta1);
    valeur = prix - pnl_vect_scalar_prod(delta1,Stauxi);
    pnl_vect_set(V,0,valeur);
    double pas = opt_->T_ / (double) mod_->hedgingDateNb_;
    double pasConstatation = opt_->T_ / (double) opt_->nbTimeSteps_;
    PnlVect *sub = pnl_vect_create(opt_->size_);
    int cpt = 0; 

    // Calcul de la valeur du portefeuille à chaque instant
    for (int j = 1; j < mod_->hedgingDateNb_ + 1; j++) {
        if ((pas*j > pasConstatation*cpt) == 1) {
            pastcopy = pnl_mat_copy(past);
            pnl_mat_resize(past, (past->m + 1), opt_->size_);
            pnl_mat_set_subblock(past,pastcopy,0,0);
            cpt ++;
        }
        pnl_mat_get_row(sub,marketPath,j);
        pnl_mat_set_row(past,sub,(past->m - 1));
        std::cout << "pas*j------------  " << pas*j << std::endl;
        price(past,pas*j,prix,ic);
        std::cout << "prix------------  " << prix << std::endl;
        delta(past,pas*j,delta2);
        pnl_vect_minus_vect(delta2, delta1);
        pnl_mat_get_row(Stauxi,past,(past->m - 1));
        std::cout << "feeeeeeeeeeeeeeee   " << pnl_vect_scalar_prod(delta2,Stauxi) << std::endl;
        valeur = pnl_vect_get(V,j-1) * exp((mod_->r_ * pas)) - pnl_vect_scalar_prod(delta2,Stauxi);
        pnl_vect_set(V,j,valeur);
        //pnl_vect_print(V);
        pnl_vect_plus_vect(delta2, delta1);    
        delta1 = pnl_vect_copy(delta2);
    }

    std::cout << "eeeeeeeeeeeeeeee   " << pnl_vect_scalar_prod(delta2,Stauxi) << std::endl;
    PnL = pnl_vect_get(V,mod_->hedgingDateNb_) + pnl_vect_scalar_prod(delta2,Stauxi) - prix;
}
