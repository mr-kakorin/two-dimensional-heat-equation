//
// Created by artoria on 5/13/20.
//

#include "CommonTools.h"

void write_out( double *T, int phi_size, int r_size, double dphi, double dr, const double R[2] )
{
    std::ofstream r_out;
    std::ofstream phi_out;
    std::ofstream t_out;
    r_out.open( "R.txt" );
    phi_out.open( "Phi.txt" );
    t_out.open( "T.txt" );
    for (int i = 0; i < phi_size; ++i)
    {
        double phi = i * dphi;
        for (int j = 0; j < r_size; ++j)
        {
            double r = j * dr + R[0];
            r_out << r << std::endl;
            phi_out << phi << std::endl;
            t_out << T[i * phi_size + j] << std::endl;
        }
    }
    r_out.close();
    phi_out.close();
    t_out.close();
}