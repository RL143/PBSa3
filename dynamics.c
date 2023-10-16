#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr;
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        dr.x = p_vectors->v[i].x * p_parameters->dt;
        dr.y = p_vectors->v[i].y * p_parameters->dt;
        dr.z = p_vectors->v[i].z * p_parameters->dt;
        p_vectors->dr[i] = dr;
        p_vectors->r[i].x += dr.x; //updating positions
        p_vectors->r[i].y += dr.y;
        p_vectors->r[i].z += dr.z;
        p_nbrlist->dr[i].x += dr.x; // update displacements with respect to neighbor list creation are updated
        p_nbrlist->dr[i].y += dr.y;
        p_nbrlist->dr[i].z += dr.z;
        p_nbrlist->dr[i].sq = (p_nbrlist->dr[i].x) * (p_nbrlist->dr[i].x) + (p_nbrlist->dr[i].y) * (p_nbrlist->dr[i].y) + (p_nbrlist->dr[i].z) * (p_nbrlist->dr[i].z);
    }
}

double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
/* Update velocities for half a time step using forces. This function returns the kinetic energy.
   v(t+dt/2) = v(t) + 0.5*F(t)/m*dt, or
   v(t+dt) = v(t+dt/2)+0.5*F(t+dt)/m*dt
*/
{
    double Ekin = 0.0;
    const double factor = 0.5 * p_parameters->mass * p_parameters->dt;
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x += factor * p_vectors->f[i].x;
        p_vectors->v[i].y += factor * p_vectors->f[i].y;
        p_vectors->v[i].z += factor * p_vectors->f[i].z;
        Ekin += p_vectors->v[i].x * p_vectors->v[i].x + p_vectors->v[i].y * p_vectors->v[i].y + p_vectors->v[i].z * p_vectors->v[i].z;
    }
    Ekin = 0.5 * Ekin * p_parameters->mass;
    return Ekin;
}

void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* Apply boundary conditions. In this of periodic BCs case particles are put back in the box */
{
    struct Vec3D invL;

    invL.x = 1.0 / p_parameters->L.x;
    invL.y = 1.0 / p_parameters->L.y;
    invL.z = 1.0 / p_parameters->L.z;

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->r[i].x -= p_parameters->L.x * floor(p_vectors->r[i].x * invL.x);
        p_vectors->r[i].y -= p_parameters->L.y * floor(p_vectors->r[i].y * invL.y);
        p_vectors->r[i].z -= p_parameters->L.z * floor(p_vectors->r[i].z * invL.z);
    }
}


/*Additional code for verifications*/
void histogram_generation(struct Parameters *p_parameters, struct Vectors *p_vectors, double *bin, double binsize, int num_bins)
/* generate histogram of velocity profile*/
{
    int bin_number;
    for (int m = 0; m < p_parameters->num_part; m++)
    {
        double vx = p_vectors->v[m].x;
        double vy = p_vectors->v[m].y;
        double vz = p_vectors->v[m].z;

        double magnitude = sqrt(vx * vx + vy * vy + vz * vz);
        bin_number = floor(magnitude / binsize);

        if (bin_number < num_bins)
            if (bin_number > 0)
                bin[bin_number] += 1;
    }
}
/*
void Radial_distribution_function(struct Parameters *p_parameters, struct Vectors *p_vectors, double dbin)
{
    int i = 0, j = 0;
    double ibin;
    struct Vec3D *r = p_vectors->r; // position
    struct Vec3D rij;
    struct Vec3D L = p_parameters->L;
    double rijabs_sq; 
    double dbin_sq = dbin * dbin; 

    // Updating
    for (i = 0; i < (p_parameters->num_part - 1); i++)
    {
        for (j = i + 1; j < p_parameters->num_part; j++)    
        {
            // Distance between the particles for each coordinate
            rij.x = r[i].x - r[j].x;
            rij.y = r[i].y - r[j].y;
            rij.z = r[i].z - r[j].z;

            // Apply periodic boundary condition for each coordinate
            rij.x = rij.x - L.x * floor((rij.x / L.x) + 0.5);
            rij.y = rij.y - L.y * floor((rij.y / L.y) + 0.5);
            rij.z = rij.z - L.z * floor((rij.z / L.z) + 0.5);

            rijabs_sq = (rij.x * rij.x) + (rij.y * rij.y) + (rij.z * rij.z);

            if (rijabs_sq < dbin_sq)
            {
                ibin = floor(sqrt(rijabs_sq) / dbin);
                p_vectors->grbin[(int)ibin] += 2.0;
            }
        }
    }
    p_parameters->grcount += 1.0;
}*/


void Radial_distribution_function(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, double *rdf, double dr)
{
    size_t num_part = p_parameters->num_part;
    double volume = p_parameters->L.x * p_parameters->L.y * p_parameters->L.z;

    for (size_t i = 0; i < num_part; ++i)
    {
        for (size_t j = 0; j < p_nbrlist->num_nbrs; ++j)
        {
            double r = sqrt(p_nbrlist->nbr[j].rij.sq);

            if (r < p_parameters->r_cut)
            {
                size_t bin = (size_t)(r / dr);
                rdf[bin] += 2.0; // Factor of 2 to account for pairs (i, j) and (j, i)
            }
        }
    }
/*
    // Normalize RDF
    for (size_t i = 0; i < num_part; ++i)
    {
        double r = (i + 0.5) * dr;
        double shell_volume = 4.0 / 3.0 * PI * (pow(r + dr, 3) - pow(r, 3));
        rdf[i] /= (shell_volume * num_part * volume);
    }*/
}


/*void density_function(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step)
{
    struct Vec3D L = p_parameters->L;
    double Maxbin = L.x; // measure positions between 0 and L.x
    int id;
    double binwidth = Maxbin / ((double)Nbins_density);
    double chi;
    struct Vec3D *r = p_vectors->r; 

    for (int i = 0; i < p_parameters->num_part; i++)
    {
        id = floor(r[i].x / binwidth);

        if (0 <= id && id < Nbins_density) // Assigning the particles A and B to different density
        {
            if (p_vectors->type[i] == 0)
                p_vectors->DAbin[id] += 1.0;
            else if (p_vectors->type[i] == 1)
                p_vectors->DBbin[id] += 1.0;
            else
                printf("type is not defined.\n");
        }
    }
    p_parameters->counter += 1.0;

    if (step == p_parameters->num_dt_steps) // Call once at end of simulation
    {
        FILE *density = NULL;
        density = fopen("density.csv", "a+");
        if (density == NULL)
        {
            printf("Filecouldnotopencorrectly\n");
            exit(1);
        }

        for (int i = 0; i < Nbins_density; i++) // Density and chi calculation
        {
            p_vectors->DAbin[i] = p_vectors->DAbin[i] / (p_parameters->counter * binwidth * L.y * L.z * (double)p_parameters->N_A);
            p_vectors->DBbin[i] = p_vectors->DBbin[i] / (p_parameters->counter * binwidth * L.y * L.z * (double)p_parameters->N_B);
            double phi_A = p_vectors->DAbin[i] / (p_vectors->DAbin[i] + p_vectors->DBbin[i]);
            chi = log((1.0 - phi_A) / phi_A) / (p_parameters->N_A * (1.0 - 2.0 * phi_A));
            fprintf(density, "%lf %lf %lf %lf %lf \n", ((double)i + 0.5) * binwidth, p_vectors->DAbin[i], p_vectors->DBbin[i], p_vectors->DAbin[i] + p_vectors->DBbin[i], chi);
        }
        fclose(density);
    }
}*/