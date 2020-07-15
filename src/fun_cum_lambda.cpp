#include <stan/math.hpp>
#include "fun_cum_lambda.h"

void fun_cum_lambda(const Eigen::MatrixXd &lambda, Eigen::MatrixXd &cum_lambda,
                    const int &I, const int &N, const Eigen::MatrixXi &NA,
                    const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                    const Eigen::MatrixXd &H)
{
    cum_lambda.setZero();
   
    for (int c = 0; c < 2; c++)
    {
        for(int i=0; i<I; i++)
        {
            for(int k=0; k<N; k++)
            {
                if (NA(i,k) == 1)
                {
                    for (int g = 0; g < seg(i,k); g++) {
                        cum_lambda(c*I + i,k) += len(g) * lambda(c*I + i,g);
                    }
                    cum_lambda(c*I + i,k) += H(i,k) * lambda(c*I + i,seg(i,k));
                }
            }
        }
    }
}
