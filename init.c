/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Bell instability setup.

  Initialize fluid variables using Eqns. of [MVBM18], sect. 4.4. using a
  1D configuration.
  The configuration is then rotated around the y- and z-axis in a way similar
  to Mignone, Tzeferacos \& Bodo, JCP (2010).
  The initial phase is computed as \f$ \phi = k_0x_0 = \vec{k}\cdot\vec{x} \f$
  (which is invariant under rotations),  where
  \f$ \vec{k} \f$ and \f$ \vec{x} \f$ are the wavevector and coordinate
  vector in the rotated system, while \f$ k_0 \f$ and \f$ x_0 \f$ are
  the original (unrotated) 1D vectors.

  The amount of rotation depends on the dimensionality of the problem
  and is uniquely specified by the domain size in the three directions
  Lx, Ly, Lz:
  \f[
    \tan(\alpha) = \frac{L_x}{L_y} = \frac{k_y}{k_x} \,,\quad
    \tan(\beta)  = \frac{L_x}{L_z} = \frac{k_z}{k_x}
  \f]
  where \f$ k_x = 2\pi/L_x,\, k_y = 2\pi/L_y\f$ and \f$ k_z = 2\pi/L_z \f$.
  In such a way the modulus of the wavevector is

  \f[
  |\vec{k}| = \sqrt{k_x^2 + k_y^2 + k_z^2}
            = 2\pi\sqrt{\frac{1}{L_x^2} + \frac{1}{L_y^2} + \frac{1}{L_z^2}}
  \f]

  In order to have one wavelength, Ly and Lz must be chosen so that
  \f$ |\vec{k}| = 2\pi \f$.
  Thus in 1D, 2D and 3D we choose:

  <CENTER>
  Dim |        Lx         |         Ly         |      Lz
  ----|-------------------|--------------------|------------
   1  |           1       |          -         |      -
   2  | \f$ \sqrt{5} \f$  | \f$ \sqrt{5}/2 \f$ |      -
   3  |           3       |        3/2         |     3/2
  </CENTER>

  The 6 configurations come in pairs (serial/parallel) and test the
  instability for eps = 0.1 in 1D (Conf. #1/#2), in 2D (Conf. #3/#4)
  and 3D (Conf. #5/#6).
  
  
  \author A. Mignone (andrea.mignone@unito.it)

  \date   June 10, 2019
  
  \b References: \n
   - [MVBM18]"A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF
               THE MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.4 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef LINEAR_SETUP
  #define LINEAR_SETUP TRUE
#endif

/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*! 
 *
 *********************************************************************** */
{
  static int first_call = 1;
   
  double thetatorad = CONST_PI/180.;
  double B0      = 1.;// g_inputParam[USER_BFIELD];
  /*
  double theta   = thetatorad * g_inputParam[USER_BTHETA];
  double phi     = thetatorad * g_inputParam[USER_BPHI];
  */
  double vamp    = g_inputParam[USER_VNOISE];
 

  if (first_call){ 
    RandomSeed (prank,0);
    first_call = 0;
  }          

  v[RHO] = 1.0;
  v[PRS] = 1.0;

  v[VX1] = vamp*RandomNumber(-1,1);
  v[VX2] = vamp*RandomNumber(-1,1);
  v[VX3] = vamp*RandomNumber(-1,1);
 
  v[BX1] = B0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;


  first_call = 0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

return; 
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    double Bm2, vA, vAmax = 100.0;
    static Data_Arr Uc;
    RBox   domBox;

    DOM_LOOP(k,j,i){
      Bm2 =   d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i]
            + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i]
            + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i];
      vA = sqrt(Bm2/d->Vc[RHO][k][j][i]);
      if (vA > vAmax) {
        d->Vc[RHO][k][j][i] = sqrt(Bm2)/vAmax;
      }
    }
    #if HAVE_ENERGY
    g_smallPressure = 1.e-3;
    DOM_LOOP(k,j,i){
      if (d->Vc[PRS][k][j][i] < g_smallPressure) {
        d->Vc[PRS][k][j][i] = g_smallPressure;
      }
    }
    #endif
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

