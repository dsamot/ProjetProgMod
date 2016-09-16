/* 
 * File:   PerformanceOption.hpp
 * Author: samotd
 *
 * Created on September 16, 2016, 11:06 AM
 */

#ifndef PERFORMANCEOPTION_HPP
#define	PERFORMANCEOPTION_HPP

class PerformanceOption {
public:
    PerformanceOption();
    PerformanceOption(const PerformanceOption& orig);
    virtual ~PerformanceOption();
    PerformanceOption(double monT, int monNbTimeStep, int maSize);
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

};

#endif	/* PERFORMANCEOPTION_HPP */

