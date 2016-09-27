/* 
 * File:   BasketOption.h
 * Author: samotd
 *
 * Created on September 16, 2016, 9:44 AM
 */

#ifndef BASKETOPTION_HPP
#define	BASKETOPTION_HPP
#include "Option.hpp"

/// \brief Basket Option: Extends Option
class BasketOption : public Option {
public:
    BasketOption();
    BasketOption(const BasketOption& orig);
    BasketOption(double monT, int monNbTimeStep, int maSize, double monStrike, PnlVect *payoffcoeff);
    virtual ~BasketOption();
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

#endif	/* BASKETOPTION_H */

