/* 
 * File:   AsianOption.hpp
 * Author: samotd
 *
 * Created on September 16, 2016, 10:40 AM
 */

#ifndef ASIANOPTION_HPP
#define	ASIANOPTION_HPP
#include "Option.hpp"

class AsianOption  : public Option{
public:
    AsianOption();
    AsianOption(double monT, int monNbTimeStep, int maSize, double monStrike);
    AsianOption(const AsianOption& orig);
    virtual ~AsianOption();
            /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @return phi(trajectoire)
     */
   double payoff(const PnlMat *path);
  
private:
    double strike;
};

#endif	/* ASIANOPTION_HPP */

