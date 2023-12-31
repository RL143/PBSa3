#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

void set_parameters(struct Parameters *p_parameters)
/* Set the parameters of this simulation */
{
// The parameters first 5 parameters are only used for demonstration puprposes
  p_parameters->kT = 1.0;                                   //thermal energy
  p_parameters->mass = 1.0;                                 //mass of a particle
  p_parameters->epsilon = 1.0;                              //LJ interaction strength

// The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_dt_steps = 200;//10000;                        //number of time steps
  p_parameters->exclude_12_nb = 0;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 0;                          // 1-3 connected atoms exluded from non-bonded interactions    
  p_parameters->dt = 0.04;                                  //integration time step
  p_parameters->L = (struct Vec3D){10,10,10};//{8, 8, 20}; //box size
    p_parameters->r_cut = 1;                              //cut-off distance used for neigbor list
  p_parameters->r_shell = 0.4;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 5;                           //number of time steps in between pdb outputs

  p_parameters->gamma = 4.5;
  p_parameters->sigma = 3.0;
  p_parameters->rho = 3.0;                      // density of the system
  p_parameters->a_AA = 25;
  p_parameters->a_BB = 25;                      // maximum repulsion
  p_parameters->a_AB = 25;//37;
  p_parameters->c = 2;
  p_parameters->num_partA = 2000;//1920;                                             // number of A particles
  p_parameters->num_partB = 2000;//1920;                                             // number of B particles
  p_parameters->num_part = p_parameters->num_partA + p_parameters->num_partB; // total number of particles
  p_parameters->N_A = 1;                                                      // Number of particles A per chain
  p_parameters->N_B = 1;                                                      // Number of particles B per chain
  p_parameters->num_chains = p_parameters->num_part / p_parameters->N_A;      // number of chains
  p_parameters->grcount = 0.0;
  p_parameters->counter = 0.0;
  strcpy(p_parameters->rad_filename, "radial_distribution.txt");

  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file


  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
