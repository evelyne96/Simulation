
//  globaldata.h
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#ifndef globaldata_h
#define globaldata_h

#include <stdio.h>
#include "list.h"

struct global_struct
    {
    //PARTICLE
    double SX, SY;
    double halfSX,halfSY;
    
    int     N_particles;
    
    double  *particle_x;                        //particle x position
    double  *particle_y;                        //particle y position
    double  *particle_fx;                       //particle total force x direction
    double  *particle_fy;                       //particle total force y direction
    int     *particle_color;                    //particle color
    double  *particle_direction;                //direction of external force

    double  particle_driving_force;             //external driving force
    double  particle_screening_length;          //inter-partcile force screening length
    double particle_screening_wavevector;
    
    double dt;          //simulation step length
    
    int total_time;     //total running time
    int echo_time;      //echo to screen
    int movie_time;     //write to file
    int time;           //current time step


    //Verlet
    double Verlet_cutoff_distance;
    double Verlet_cutoff_distance_squared;
    double Verlet_intershell_squared;

    int N_Verlet;
    int N_Verlet_max; //initial allocation + later, longest allocation
    int *Verletlisti;
    int *Verletlistj;
    int flag_to_rebuild_Verlet;

    double *particle_dx_so_far;
    double *particle_dy_so_far;

    List **Verlet_cell_list;
    int Nx_Verlet_cell_list;
    int Ny_Verlet_cell_list;
    int Verlet_cell_list_size;
    int Verlet_cell_grid_capacity;
    
    FILE *moviefile;

    //PINNINGSITES

    int N_pinningsites;
    double *pinningsite_x;
    double *pinningsite_y;
    double *pinningsite_fx;
    double *pinningsite_fy;
    int *pinningsite_color;
    double *pinningsite_direction_x;
    double *pinningsite_direction_y;
    double *pinningsite_dx_so_far;
    double *pinningsite_dy_so_far;
    double *pinningsite_R;
    double *pinningsite_fmax;
    
    double pinning_lattice_constant;
    double pinning_driving_force;
    int pinning_direction_change;

    double pinningsite_setradius;
    double pinningsite_grid_dx;
    double pinningsite_grid_dy;
    int Nx_pinningsite_grid;
    int Ny_pinningsite_grid;
    int **pinningsite_grid;

    double *particle_dist_so_far;
    };

extern struct global_struct global;

#endif /* globaldata_h */
