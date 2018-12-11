//  running.h
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#ifndef running_h
#define running_h

#include <stdio.h>

void run_simulation(void);

void calculate_external_forces_on_particles(void);
void calculate_pairwise_forces(void);

void move_particles(void);
void fold_particle_back_PBC(int i);

void distance_squared_folded_PBC(double x0,double y0,double x1,double y1,
        double *r2_return, double *dx_return,double *dy_return);

void write_cmovie_frame(void);

void rebuild_Verlet_list(void);
void check_Verlet_rebuild_condition_and_set_flag(void);
void calculate_pairwise_forces_verlet(void);

void run_simple_time(void);
void run_verlet_time(void); 
void run_verlet_simulation(void);
void run_verlet_cell_time(void);
void run_verlet_cell(void);

void update_verlet_list_from_cell(int one, int two, int same);
void rebuild_Verlet_list_with_cell(void);
void rebuild_Verlet_cell(void);

//PINNING
void calculate_pinning_forces(void);
void rebuild_pinning_grid(void);
void calculate_neighbour_pinning_sites(int i);
void calculate_pinning_forces_slow(void);
void run_verlet_simulation_slow(void);

void write_statistics();

#endif /* running_h */
