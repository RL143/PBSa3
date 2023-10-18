/* Overall changes for all forces are implemented here*/
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
        // Initialize the forces to zero
        f[i] = (struct Vec3D){0.0, 0.0, 0.0}; /*initialize forces to zero*/
    //double Epot = 0;    // initialize for only non-concervative forces
    double Epot = calculate_conservative_force(p_parameters, p_nbrlist, p_vectors);
    calculate_dissipative_force(p_parameters, p_nbrlist, p_vectors);
    calculate_random_force(p_parameters, p_nbrlist, p_vectors);
    Epot += calculate_spring_force(p_parameters, p_nbrlist, p_vectors);
    return Epot;
}

double calculate_conservative_force(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D fC;
    double r_cutsq,  a, Fc;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;
    double Epot = 0.0;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;

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
            if (p_vectors->type[i] == 0 && p_vectors->type[j] == 0)
            {
                a = p_parameters->a_AA;
            }
            else if (p_vectors->type[i] == 1 && p_vectors->type[j] == 1)
                a = p_parameters->a_BB;
            else
                a = p_parameters->a_AB;
            
            //Storing reoccuring calculations as factors in double
            Fc= a*(1-sqrt(rij.sq))/sqrt(rij.sq);

            //Calculate conservative force in the x,y and z directions
            fC.x = Fc* rij.x;
            fC.y = Fc* rij.y;
            fC.z = Fc* rij.z;

            f[i].x += fC.x;
            f[i].y += fC.y;
            f[i].z += fC.z;
            f[j].x -= fC.x;
            f[j].y -= fC.y;
            f[j].z -= fC.z;

            Epot += -a*sqrt(rij.sq) + 0.5* a* rij.sq + a/2;
        }
    }
    return Epot;
}

double calculate_dissipative_force(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D fD;
    double r_cutsq, Fd;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;
    double Epot = 0.0;
    double gamma = p_parameters->gamma;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;

    for (size_t k = 0; k < num_nbrs; k++)
    {
        // for each pair in the neighbor list compute the pair forces
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        if (rij.sq < r_cutsq)
        // Compute forces if the distance is smaller than the cutoff distance
        {
            //Storing reoccuring calculations as factors in double
            double inp = ((p_vectors->v[i].x-p_vectors->v[j].x)*rij.x+(p_vectors->v[i].y-p_vectors->v[j].y)*rij.y
                            +(p_vectors->v[i].z-p_vectors->v[j].z)*rij.z)/sqrt(rij.sq);

            Fd = -gamma*pow(1-sqrt(rij.sq),2);

            //Calculate the dissipative force in the x,y and z directions
            fD.x = Fd *inp*(rij.x/sqrt(rij.sq));
            fD.y = Fd *inp*(rij.y/sqrt(rij.sq));
            fD.z = Fd *inp*(rij.z/sqrt(rij.sq));

            f[i].x += fD.x;
            f[i].y += fD.y;
            f[i].z += fD.z;
            f[j].x -= fD.x;
            f[j].y -= fD.y;
            f[j].z -= fD.z;
        }
    }
}

double calculate_random_force(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D fR;
    double r_cutsq, Fr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;
    double Epot = 0.0;
    double sigma = p_parameters->sigma;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;

    for (size_t k = 0; k < num_nbrs; k++)
    {
        // for each pair in the neighbor list compute the pair forces
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        if (rij.sq < r_cutsq)
        // Compute forces if the distance is smaller than the cutoff distance
        {
            //Storing reoccuring calculations as factors in double
            Fr = sigma*(1-sqrt(rij.sq))*(1/sqrt(p_parameters->dt))*generate_uniform_random();

            //Calculate the random force in the x,y and z directions
            fR.x = Fr*(rij.x/sqrt(rij.sq));
            fR.y = Fr*(rij.y/sqrt(rij.sq));
            fR.z = Fr*(rij.z/sqrt(rij.sq));

            f[i].x += fR.x;
            f[i].y += fR.y;
            f[i].z += fR.z;
            f[j].x -= fR.x;
            f[j].y -= fR.y;
            f[j].z -= fR.z;
        }
    }
}

double calculate_spring_force(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{ 
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    int C = p_parameters->c;
    struct Vec3D rij;
    struct Vec3D fS = {0};
    double Epot = 0.0;
    for (size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        fS.x = -C * rij.x;
        fS.y = -C * rij.y;
        fS.z = -C * rij.z;

        Epot -= 0.5 * fS.x * rij.x;
        Epot -= 0.5 * fS.y * rij.y;
        Epot -= 0.5 * fS.z * rij.z;

        f[i].x += fS.x;
        f[i].y += fS.y;
        f[i].z += fS.z;
        f[j].x -= fS.x;
        f[j].y -= fS.y;
        f[j].z -= fS.z;
    }
    return Epot;
}