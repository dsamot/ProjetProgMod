#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"

#include "BlackScholesModel.hpp"

using namespace std;

int main()
{
    double T = 365;
    int nbTimeSteps = 365;
    PnlMat *path = pnl_mat_create(366,2);
    PnlVect *sigma = pnl_vect_create(2);
    PnlVect *spot = pnl_vect_create(2);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    /*int M = 1E5;
    int dim = 2;
    pnl_rng_sseed(rng, time(NULL));

    double acc = 0., var = 0;

    for (int i = 0; i < M; i++)
    {
        pnl_vect_rng_normal(G, dim, rng);
        double tmp = pnl_vect_norm_two(G);
        acc += tmp;
        var += tmp * tmp;
    }

    acc /= M;
    var = var / M - acc * acc;

    cout << "E[||G||_2] = " << acc << endl;
    cout << "Var(||G||_2) = " << var << endl;

    pnl_vect_free(&G);
    pnl_rng_free(&rng);*/
    BlackScholesModel *model = new BlackScholesModel(2,1.5,0,sigma,spot);
    model->asset(path, T, nbTimeSteps, rng);
    return 0;
}
