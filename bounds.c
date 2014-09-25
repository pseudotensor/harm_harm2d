/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/


#include "decs.h"

/* bound array containing entire set of primitive variables */

void bound_prim( double prim[][N2+4][NPR] )
{
        int i,j,k ;
	void inflow_check(double *pr, int ii, int jj, int type );
        struct of_geom geom ;

        /* inner r boundary condition: u, gdet extrapolation */
        for(j=0;j<N2;j++) {
          PLOOP prim[-1][j][k] = prim[0][j][k] ;
          PLOOP prim[-2][j][k] = prim[0][j][k] ;
          pflag[-1][j] = pflag[0][j] ;
          pflag[-2][j] = pflag[0][j] ;

          PLOOP prim[N1  ][j][k] = prim[N1-1][j][k] ;
          PLOOP prim[N1+1][j][k] = prim[N1-1][j][k] ;
          pflag[N1  ][j] = pflag[N1-1][j] ;
          pflag[N1+1][j] = pflag[N1-1][j] ;

        }

        
        /* polar BCs */
        for(i=-2;i<=N1+1;i++) { 
          PLOOP prim[i][-1  ][k] = prim[i][0   ][k] ;
          PLOOP prim[i][-2  ][k] = prim[i][0   ][k] ;
          pflag[i][-1  ] = pflag[i][0   ] ;
          pflag[i][-2  ] = pflag[i][0   ] ;

          PLOOP prim[i][N2  ][k] = prim[i][N2-1][k] ;
          PLOOP prim[i][N2+1][k] = prim[i][N2-1][k] ;
          pflag[i][N2  ] = pflag[i][N2-1] ;
          pflag[i][N2+1] = pflag[i][N2-1] ;
        }

}

void inflow_check(double *pr, int ii, int jj, int type )
{
	    return ;

	}

