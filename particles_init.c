/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Initialize particles for the Bell instability test-
 
 \authors A. Mignone (andrea.mignone@unito.it)\n
 
 \date    Aug 17, 2020
 
 \b References: \n
    - [MVBM18] "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF
                THE MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.4 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "rotate.h"

#ifndef LINEAR_SETUP
 #define LINEAR_SETUP TRUE
#endif

double uni_ran();

void MonoenergeticDistribution(double v[3], double vbst, double vcr);

/* ********************************************************************* */
void Particles_Init(Data *d, Grid *grid)
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{
  int i,j,k, np, dir;
  int np_glob = RuntimeGet()->Nparticles_glob;
  int np_cell = RuntimeGet()->Nparticles_cell;
  double xbeg[3], xend[3], vbeg[3], vend[3];
  double gamma, c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  Particle p;

/* --------------------------------------------------------------
    Global initialization
   -------------------------------------------------------------- */

  if (np_glob > 0){
    print ("! Use cell-by-cell initialization\n");
    QUIT_PLUTO(1);
  }

/* --------------------------------------------------------------
    Cell by cell initialization
   -------------------------------------------------------------- */

  if (np_cell > 0){

    static int first_call = 1;
    double rhoCR  = g_inputParam[USER_NCR_by_NGAS]/np_cell;
    double uCR    = g_inputParam[USER_UCR]; 
    double vCR;

    double gamma;
    gamma = sqrt( 1. + pow(uCR/PARTICLES_CR_C,2.) );
    vCR = sqrt(1. - 1./(gamma*gamma)) * PARTICLES_CR_C; // vCR in unit of vA;

    double pcr[3], vboost;

    vboost = g_inputParam[USER_VBOOST];
    double Lambda = 2.*CONST_PI/(0.5 * g_inputParam[USER_NCR_by_NGAS]  * vboost); // wavelength;

    DOM_LOOP(k,j,i){

      xbeg[IDIR] = grid->xl[IDIR][i]; xend[IDIR] = grid->xr[IDIR][i];
      xbeg[JDIR] = grid->xl[JDIR][j]; xend[JDIR] = grid->xr[JDIR][j];
      xbeg[KDIR] = grid->xl[KDIR][k]; xend[KDIR] = grid->xr[KDIR][k];

      /* -- Loop on particles -- */
  
      for (np = 0; np < np_cell; np++){
	MonoenergeticDistribution(pcr, vboost, vCR);
        Particles_LoadUniform(np, np_cell, xbeg, xend, p.coord);
        p.speed[IDIR] = pcr[0];
        p.speed[JDIR] = pcr[1];
        p.speed[KDIR] = pcr[2];

      /* -- Assign four-velocity -- */

	
        p.mass        = rhoCR*grid->dV[k][j][i];
        p.color       = 1;

        Particles_Insert (&p, d, PARTICLES_CREATE, grid);
      }
    }

    if (first_call){

      print ("> Particles_Init():\n");
      print ("  ------------------------------------ \n");
      print ("  normal               = [%f, %f, %f]\n", p.speed[IDIR]/uCR,
                                                        p.speed[JDIR]/uCR,
                                                        p.speed[KDIR]/uCR);
      print ("  gyration radius      = %8.3e\n", uCR/PARTICLES_CR_E_MC);
      print ("  most unstable lambda = %f\n", Lambda);
      print ("  ------------------------------------ \n");
    }
  }

  Particles_SetID(d->PHead);
}
/* ********************************************************************* */
void Particles_Inject(Data *data, Grid *grid)
/*!
 *  Sets user-defined boundary conditions on particles.
 *
 *  \param [in]  data    Pointer to the PLUTO data structure.
 *  \param [in]  grid    Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Particles_UserDefBoundary(Data *d, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{

  int    dir;
  double xbeg[3], xend[3];
  particleNode *curr = d->PHead, *next;
  Particle *p;
  
  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  if (side == X1_BEG){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X1_END){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_BEG){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_END){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_BEG){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_END){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

}


double uni_ran()
{
   double x;

   x = (1.0 * ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. ));

   return(x);
}

void MonoenergeticDistribution(double v[3], double vbst, double vcr)
{

  double gamcr, pcr;
  double c2 = PARTICLES_CR_C*PARTICLES_CR_C;

  gamcr = 1./sqrt(1. - vcr*vcr/c2);
  pcr = gamcr*vcr;

  double px_prime, py_prime, pz_prime;
  double mu, phi;

  double ran;
  ran = rand();

  mu  = 2.*uni_ran()-1.;
  phi = 2.*CONST_PI*uni_ran();

  px_prime = pcr * mu;
  py_prime = pcr * sqrt(1.-(mu*mu))*cos(phi);
  pz_prime = pcr * sqrt(1.-(mu*mu))*sin(phi);

  double px, py, pz;

  px = (px_prime + gamcr*vbst)/sqrt(1.- (vbst*vbst)/c2);
  py = py_prime;
  pz = pz_prime;

  v[0] = px;
  v[1] = py;
  v[2] = pz;

  return;
}

