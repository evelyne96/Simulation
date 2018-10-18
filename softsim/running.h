//
//  running.h
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//

#ifndef running_h
#define running_h

#include <stdio.h>

void run_simulation(void);

void calculate_external_forces_on_pinningsites(void);
void move_pinningsites(void);
void fold_pinningsite_back_PBC(int i);
void rebuild_pinning_grid(void);

void calculate_external_forces_on_particles(void);
void calculate_pairwise_forces(void);
void calculate_pairwise_forces_simple(void);
void calculate_forces_between(int i, int j);
void calculate_pairwise_forces_Verlet_lookup_cell();
void check_Verlet_rebuild_condition_and_set_flag(void);
void rebuild_Verlet_list(void);
void move_particles(void);
void fold_particle_back_PBC(int i);

void distance_squared_folded_PBC(double x0,double y0,double x1,double y1,
        double *r2_return, double *dx_return,double *dy_return);

void write_cmovie_frame(void);

void adjust_pinningsite_directions(void);
void rotate_pinningsite_directions(void);

void go_through_linked_Verlet_list(int one, int two, int same);

//functions for testing the program
void test_program_by_coloring(void);

//test how much Verlet runs
void run_test_verlet_vs_time(void);

//verlet cell
void rebuild_Verlet_cell(void);
void rebuild_Verlet_list_with_cell(void);

#endif /* running_h */
