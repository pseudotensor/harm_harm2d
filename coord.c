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

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

/* should return boyer-lindquist coordinte of point */
void bl_coord(double *X, double *r, double *th)
{

  *r = X[1];
  *th = X[2];

	
	return ;
}

/* insert metric here */
void gcov_func(double *X, double gcov[][NDIM])
{
        int j,k ;
        double sth,cth,s2,rho2 ;
	double r,th ;
	double tfac,rfac,hfac,pfac ;

        DLOOP gcov[j][k] = 0. ;


    gcov[TT][TT] = -1.0;
    gcov[1][1] = 1.0;
    gcov[2][2] = 1.0;
    gcov[3][3] = 1.0;
}

/* some grid location, dxs */
void set_points()
{
        startx[1] = -1.0 ;
        startx[2] = 0. ;
        dx[1] = 2.0/N1 ;
        dx[2] = 1./N2 ;
}

void fix_flux(double F1[][N2+4][NPR], double F2[][N2+4][NPR])
{
	int i,j,k ;


	return;
}


void rescale(double *pr, int which, int dir, int ii, int jj, int face, struct of_geom *geom)
{
	double scale[NPR], r, th, X[NDIM];
	int k ;

    // do nothing to pr
}
