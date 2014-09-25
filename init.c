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

/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


void coord_transform(double *pr,int i, int j) ;

void init()
{
	int i,j ;
	double r,th,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax, lfish_calc(double rmax) ;

	/* some physics parameters */
	gam = 4./3. ;

	/* disk parameters (use fishbone.m to select new solutions) */
    a = 0.0 ;

    /* some numerical parameters */
    lim = DONOR ;
    failed = 0 ;	/* start slow */
    cour = 0.1 ;
    cour=0.3;
    //    dt = 1E-6 ;
    dt = 5E-5;
	R0 = 0.0 ;
    Rin = 1E10;
    Rout = 1E10+2.0;

        t = 0. ;
        hslope = 0.0 ;

        set_arrays() ;
        set_grid() ;


        /* output choices */
        //tf = 0.0255 ;
        tf=2.5;
	DTd = 2.5/100.0;	/* dumping frequency, in units of M */
	DTl = DTd ;	/* logfile frequency, in units of M */
	DTi = DTd ; 	/* image file frequ., in units of M */
	DTr = 100 ; 	/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
	defcon = 1. ;

#define FTYPE double


    FTYPE pleft[NPR], pright[NPR], P;


  
    //zero out initial conditions
    int pl;
    for(pl=RHO;pl<=B3;pl++) pleft[pl] = 0.;
    for(pl=RHO;pl<=B3;pl++) pright[pl] = 0.;
    
    int WHICHKOMI=1;
    
    //fast shock
    if(WHICHKOMI==1){
      //left state
      pleft[U1] = 25.0;
      pleft[U2] =  0.0;
      pleft[U3] =  0.0;
      pleft[B1] = 20.0;
      pleft[B2] = 25.02;
      pleft[B3] =  0.0;
      P = 1.0;
      pleft[RHO] = 1.;
      pleft[UU] = P/(gam-1);

      //right state
      pright[U1] = 1.091;
      pright[U2] =  0.3923;
      pright[U3] =  0.0;
      pright[B1] = 20.0;
      pright[B2] = 49.0;
      pright[B3] =  0.0;
      P = 367.5;
      pright[RHO] = 25.48;
      pright[UU] = P/(gam-1);
    }


	ZSLOOP(0,N1-1,0,N2-1) {
      coord(i,j,CENT,X) ;
      bl_coord(X,&r,&th) ;
      FTYPE x = r;

      if(x<=0){
        p[i][j][RHO] = pleft[RHO] ;
        p[i][j][UU] = pleft[UU] ;
        p[i][j][U1] = pleft[U1] ;
        p[i][j][U2] = pleft[U2] ;
        p[i][j][U3] = pleft[U3] ;
        p[i][j][B1] = pleft[B1] ;
        p[i][j][B2] = pleft[B2] ;
        p[i][j][B3] = pleft[B3] ;
      }      
      else{
        p[i][j][RHO] = pright[RHO] ;
        p[i][j][UU] = pright[UU] ;
        p[i][j][U1] = pright[U1] ;
        p[i][j][U2] = pright[U2] ;
        p[i][j][U3] = pright[U3] ;
        p[i][j][B1] = pright[B1] ;
        p[i][j][B2] = pright[B2] ;
        p[i][j][B3] = pright[B3] ;
      }      
	}

	bound_prim(p) ;


}

/* this version starts w/ BL 4-velocity and
 * converts to 3-velocities in modified
 * Kerr-Schild coordinates */

void coord_transform(double *pr,int ii, int jj)
{
        double X[NDIM],r,th,ucon[NDIM],trans[NDIM][NDIM],tmp[NDIM] ;
        double AA,BB,CC,discr ;
	double alpha,gamma,beta[NDIM] ;
	struct of_geom geom ;
	struct of_state q ;
	int i,j,k ;
}

double lfish_calc(double r)
{
  return(0.0);

}


