/******************************************************************************/
/*                                                    			              */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 4.0, 18/9/2018    			              */
/******************************************************************************/

/*
 * For the 2023 PBS assigment the code needs to be extended
 *
 * -Implement a Berendsen thermostat in dynamics.c
 * -Implement bonds in initialise_bonds in file initialise.c
 * -Initialize vectors.type such that particles get the proper type 
 * -Implement the bonded and non-bonded force in force.c. (Make the forces type dependent)
 * -Change the particle position initialisation such that it takes into account bond lengths and angles
 * -Implement the needed on-the-fly data analysis
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"
// #include "density.h"

/**
 * @brief main The main of the MD code. After initialization,
 * a velocity-Verlet scheme is executed for a specified number of time steps.
 *
 * @return int 0 if successful
 */
int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    size_t step;
    double Ekin, Epot, time;
    FILE *Energy = NULL;
    FILE *histogram = NULL;
    FILE *RDF = NULL;
    int num_bin = 100;
    double bin[num_bin]; // Number of bins for the simulation
    for (int i = 0; i < num_bin; ++i)
        bin[i] = 0.0;  

    set_parameters(&parameters);
    alloc_memory(&parameters, &vectors, &nbrlist);
    if (parameters.load_restart == 1)
    {
        load_restart(&parameters, &vectors);
        initialise_structure(&parameters, &vectors, &nbrlist);
        step = 0;
        time = 0.0;
    }
    else
    initialise(&parameters, &vectors, &nbrlist, &step, &time);
    boundary_conditions(&parameters, &vectors);
    build_nbrlist(&parameters, &vectors, &nbrlist);
    
    Epot = calculate_forces(&parameters, &nbrlist, &vectors);
    record_trajectories_pdb(1, &parameters, &vectors, time);

    double vmax = sqrt(2 * (10 * parameters.kT) / parameters.mass);
    double binsize = vmax / num_bin;
    struct Vec3D L = parameters.L;
    double dbin = 0.5 * L.x / (double)Nbins_radial;                

    // Initialize
    for (size_t i = 0; i < (double)Nbins_radial; i++) // Loop over bins
    {
        vectors.grbin[i] = 0.0; // Initialize radial distribution to 0
    }

    while (step < parameters.num_dt_steps) // start of the velocity-Verlet loop
    {
        step++;
        time += parameters.dt;
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);
        update_positions(&parameters, &nbrlist, &vectors);
        boundary_conditions(&parameters, &vectors);
        update_nbrlist(&parameters, &vectors, &nbrlist);
        Epot = calculate_forces(&parameters, &nbrlist, &vectors);
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);
        
        if (step > 3 * parameters.num_dt_steps / 4){
          //  Radial_distribution_function(&parameters, &vectors, dbin);
            density_function(&parameters, &vectors, step);
        }
        if (step % parameters.num_dt_pdb == 0)
            record_trajectories_pdb(0, &parameters, &vectors, time);
        if (step % parameters.num_dt_restart == 0)
            save_restart(&parameters, &vectors);
        if (step > 3 * parameters.num_dt_steps / 4)
            histogram_generation(&parameters, &vectors, bin, binsize, num_bin);
            printf("Step %zu, Time %f, Epot %f, Ekin %f, Etot %f\n", step, time, Epot, Ekin, Epot + Ekin);
        if (step % parameters.num_dt_pdb == 0)
            record_trajectories_pdb(0, &parameters, &vectors, time);
        if (step % parameters.num_dt_restart == 0)
            save_restart(&parameters, &vectors);
        Energy = fopen("Energy.csv", "a+");
        if (Energy == NULL){
            printf("Filecouldnotopencorrectly\n");
            exit(1);
        }
        fprintf(Energy, "%E, %E, %E, %E\n", time, Epot, Ekin, Epot + Ekin);
        fclose(Energy);
    }

    histogram = fopen("Histogram.csv", "a+");
    if (histogram == NULL)
    {
        printf("Filecouldnotopencorrectly\n");
        exit(1);
    }
    for (int i = 0; i < num_bin; i++)
        fprintf(histogram, "%d, %E\n", i, bin[i]);
    fclose(histogram);


    double volume;
    double rho_rdf = (parameters.num_part - 1) / (L.x * L.y * L.z); // Calculate average particle density (other particles are N-1)
    double *grbin = vectors.grbin;
  
    RDF = fopen("RDF.csv", "a+");
    if (RDF == NULL)
    {
        printf("Filecouldnotopencorrectly\n");
        exit(1);
    }
    fprintf(RDF, "%f %f\n", 0.0, 0.0);
    for (size_t ibin = 0; ibin < Nbins_radial - 1 ; ibin++)
    {
        volume = (4.0 / 3.0) * PI * (((ibin + 1) * (ibin + 1) * (ibin + 1)) - (ibin * ibin * ibin)) * (dbin * dbin * dbin); 
        grbin[ibin] = grbin[ibin] / (rho_rdf * volume * parameters.grcount * parameters.num_part);  // Normalization over mass, steps, and particles
        fprintf(RDF, "%f %f\n", ((double)ibin + 0.5) * dbin, grbin[ibin]);
    }
    fclose(RDF);

    save_restart(&parameters, &vectors);
    free_memory(&vectors, &nbrlist);

    return 0;
}
