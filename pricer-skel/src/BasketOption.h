/* 
 * File:   BasketOption.h
 * Author: samotd
 *
 * Created on September 16, 2016, 9:44 AM
 */

#ifndef BASKETOPTION_H
#define	BASKETOPTION_H
#include "Option.hpp"

class BasketOption : public Option {
public:
    BasketOption();
    BasketOption(const BasketOption& orig);
    BasketOption(double monT, int monNbTimeStep, int maSize, double monStrike);
    virtual ~BasketOption();
        /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset.
     * @return phi(trajectoire)
     */
    virtual double payoff(const PnlMat *path) = 0;
private:
    double strike;
};

#endif	/* BASKETOPTION_H */

