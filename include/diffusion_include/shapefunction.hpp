#ifndef _SHAPE_FUNCTION_H_
#define _SHAPE_FUNCTION_H_

#include<vector>

class ShapeFunction2D{
    public:
        static void C2D4_N(std::vector<double> &N,const double &g1,const double &g2)
        {
            N[0] = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
            N[1] = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
            N[2] = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
            N[3] = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
        }
        static void C2D4_dNdr(std::vector<std::vector<double>> &dNdr,const double &g1,const double &g2)
        {
          dNdr[0][0] = -2.5e-1 * (1e+0-g2);
          dNdr[0][1] = -2.5e-1 * (1e+0-g1);
          dNdr[1][0] =  2.5e-1 * (1e+0-g2);
          dNdr[1][1] = -2.5e-1 * (1e+0+g1);
          dNdr[2][0] =  2.5e-1 * (1e+0+g2);
          dNdr[2][1] =  2.5e-1 * (1e+0+g1);
          dNdr[3][0] = -2.5e-1 * (1e+0+g2);
          dNdr[3][1] =  2.5e-1 * (1e+0-g1);
        }
};

#endif