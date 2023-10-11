#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "random.h"

double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;
    for (size_t i = 0; i < num_part; i++)
        // initialize the forces to zero
        f[i] = (struct Vec3D){0.0, 0.0, 0.0}; /*initialize forces to zero*/

    double Epot = calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_angle(p_parameters, p_vectors);
    Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);
    return Epot;
}

double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0;
    struct Bond * bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi={0};
    for(size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x*floor(rij.x/L.x+0.5); //apply minimum image convenction for bonded particles
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y*floor(rij.y/L.y+0.5); 
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z*floor(rij.z/L.z+0.5);

        /*
            Here provide the force calculation
        */

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;
    }
    return Epot;
}

double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0;
    struct Angle * angles = p_vectors->angles;
    size_t num_angles = p_vectors->num_angles;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj;
    struct Vec3D fi={0}, fk={0};
    for(size_t q = 0; q < num_angles; ++q)
    {
        size_t i = angles[q].i;
        size_t j = angles[q].j;
        size_t k = angles[q].k;

        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x*floor(rij.x/L.x+0.5); //apply minimum image convenction for bonded particles
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y*floor(rij.y/L.y+0.5); 
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z*floor(rij.z/L.z+0.5);

        rkj.x = r[k].x - r[j].x;
        rkj.x = rkj.x - L.x*floor(rkj.x/L.x+0.5); //apply minimum image convenction for bonded particles
        rkj.y = r[k].y - r[j].y;
        rkj.y = rkj.y - L.y*floor(rkj.y/L.y+0.5); 
        rkj.z = r[k].z - r[j].z;
        rkj.z = rkj.z - L.z*floor(rkj.z/L.z+0.5);

        /*
            Here provide the force calculation
        */

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= (fi.x + fk.x);
        f[j].y -= (fi.y + fk.y);
        f[j].z -= (fi.z + fk.z);
        f[k].x += fk.x;
        f[k].y += fk.y;
        f[k].z += fk.z;
    }
    return Epot;
}

double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0;
    struct Dihedral * dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    for(size_t q = 0; q < num_dihedrals; ++q)
    {
        size_t i = dihedrals[q].i;
        size_t j = dihedrals[q].j;
        size_t k = dihedrals[q].k;
        size_t l = dihedrals[q].l;

    }
    return Epot;
}


double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
/* Compute non-bonded forces on particles using the pairs in a neighbor list.
This function returns the total potential energy of the system. */
{
    struct Vec3D fC, fD, fR;
    double r_cutsq, sigmasq, sr2, sr6, sr12, fr, prefctr, aij, Fr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;
    sigmasq = p_parameters->sigma * p_parameters->sigma;
    double epsilon = p_parameters->epsilon;

    double Epot = 0.0;//, Epot_cutoff;
    //sr2 = sigmasq / r_cutsq;
    //sr6 = sr2 * sr2 * sr2;
    //sr12 = sr6 * sr6;
    //Epot_cutoff = sr12 - sr6;

    for (size_t k = 0; k < num_nbrs; k++)
    {
        // for each pair in the neighbor list compute the pair forces
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        if (rij.sq < r_cutsq)
        // Compute forces if the distance is smaller than the cutoff distance
        {
            //Load maximum repulsion parameter
            aij = p_parameters->aij;
            
            //Storing reoccuring calculations in double
            Fr= aij*(1-sqrt(rij.sq))/sqrt(rij.sq);

            //Calculate conservative force in the x,y and z directions
            fC.x = Fr* rij.x;
            fC.y = Fr* rij.y;
            fC.z = Fr* rij.z;
            
            //Calculate the dissipative force in the x,y and z directions
            fD.x = -p_parameters->gamma*pow(1-sqrt(rij.sq),2)*(rij.x/sqrt(rij.sq)*(p_vectors->v[i].x-p_vectors->v[j].x))*(rij.x/sqrt(rij.sq));
            fD.y = -p_parameters->gamma*pow(1-sqrt(rij.sq),2)*(rij.y/sqrt(rij.sq)*(p_vectors->v[i].y-p_vectors->v[j].y))*(rij.y/sqrt(rij.sq));
            fD.z = -p_parameters->gamma*pow(1-sqrt(rij.sq),2)*(rij.z/sqrt(rij.sq)*(p_vectors->v[i].z-p_vectors->v[j].z))*(rij.z/sqrt(rij.sq));

            //Calculate the random force in the x,y and z directions
            double factor = p_parameters->sigma*sqrt(1-sqrt(rij.sq))*pow(p_parameters->dt,-0.5)*generate_uniform_random();
            fR.x = factor*(rij.x/sqrt(rij.sq));
            fR.y = factor*(rij.y/sqrt(rij.sq));
            fR.z = factor*(rij.z/sqrt(rij.sq));

            Epot += -aij*sqrt(rij.sq) + 0.5* aij* rij.sq - aij/2;
        }

       else
       {
            fC.x = 0;
            fC.y = 0;
            fC.z = 0;
            fD.x = 0;
            fD.y = 0;
            fD.z = 0;
            fR.x = 0;
            fR.y = 0;
            fR.z = 0;
        }
            
        //Add to overall forces
        f[i].x += fC.x;// + fD.x; //+ fR.x;
        f[i].y += fC.y;// + fD.y;  //+ fR.y;
        f[i].z += fC.z;// + fD.z;  //+ fR.z;
        f[j].x -= fC.x;// + fD.x;   //+ fR.x;
        f[j].y -= fC.y;// + fD.y;   //+ fR.y;
        f[j].z -= fC.z;// + fD.z;   //+ fR.z;

        //Calculate contribution to potential energy
        //Epot += -aij*sqrt(rij.sq) + 0.5* aij* rij.sq - aij/2;

    }
    return Epot;
}
