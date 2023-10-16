#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/**
 * @brief Update particle positions by using velocities.
 * @param[in] p_parameters member: dt
 * @param[out] p_nbrlist used members: dr
 * @param[in,out] p_vectors members r, dr, v
 */
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Update velocities for half a time step using forces.
 * @param[in] p_parameters used members: mass, dt
 * @param[in] p_nbrlist
 * @param[in, out] p_vectors used members: v, f
 * @return double kinetic energy
 */
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Apply boundary conditions: particles folded back in periodic box.
 * @param[in] p_parameters used members: L
 * @param[in, out] p_vectors used members: r
 * @todo include neighbor list and get box size from there
 */
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Apply thermostat by manipolating particle velocities
 * @param[in] p_parameters parameters for the thermostat
 * @param[in, out] p_vectors used members: v
 * @param[in] Ekin current kinetic energy
 */

void histogram_generation(struct Parameters *p_parameters, struct Vectors *p_vectors, double *bin, double binsize, int num_bin);
void Radial_distribution_function(struct Parameters *p_parameters, struct Vectors *p_vectors, double dbin);
void density_function(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step);

#endif /* DYNAMICS_H_ */