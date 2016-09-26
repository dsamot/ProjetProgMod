/* 
 * File:   PerformanceOption.hpp
 * Author: samotd
 *
 * Created on September 16, 2016, 11:06 AM
 */

#ifndef PERFORMANCEOPTION_HPP
#define	PERFORMANCEOPTION_HPP
 #include "Option.hpp"

class PerformanceOption : public Option{
public:
    PerformanceOption();
    PerformanceOption(const PerformanceOption& orig);
    virtual ~PerformanceOption();
    PerformanceOption(double monT, int monNbTimeStep, int maSize, PnlVect *payoffcoeff);
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

};

#endif	/* PERFORMANCEOPTION_HPP */

