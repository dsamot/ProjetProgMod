/* 
 * File:   AsianOption.cpp
 * Author: samotd
 * 
 * Created on September 16, 2016, 10:40 AM
 */

#include "AsianOption.hpp"
#include <iostream>

AsianOption::AsianOption() {
}

AsianOption::AsianOption(const AsianOption& orig) {
}

AsianOption::AsianOption(double monT, int monNbTimeStep, int maSize, double monStrike, PnlVect *payoffcoeff){
   T_ = monT;
   nbTimeSteps_ = monNbTimeStep;
   size_ = maSize;
   strike = monStrike;
   payoffCoeff_ = payoffcoeff;
}

AsianOption::~AsianOption() {
}

double AsianOption::payoff(const PnlMat* path) {
  double sommeGenerale = 0;
  double sommeInterne = 0;
  //PnlMat* poids = pnl_mat_create(1, path->n);
  int N = (path->m) - 1;
  //double payoffcoeff = 1 / (double)(path->n);
  /*for (int j = 0; j < path->n; j++) {
    pnl_mat_set(poids,0,j,payoffcoeff);
  }*/

  for (int j = 0; j < path->n; j++) {
    for (int i = 0; i < path->m; i++) {
        sommeInterne += pnl_mat_get(path,i,j);
    }
    sommeGenerale = sommeGenerale + pnl_vect_get(payoffCoeff_,j)*(1.0 / (double)(nbTimeSteps_ + 1)) * sommeInterne;
    sommeInterne = 0;
  }
  sommeGenerale = sommeGenerale - strike;
  //pnl_mat_free(&poids);
  if (sommeGenerale > 0) {
    return sommeGenerale;
  } else {
    return 0.0;
  }

}   