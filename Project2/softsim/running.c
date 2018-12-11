//  running.h
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#include "running.h"
#include "globaldata.h"
#include "timing.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "list.h"
#include "initializer.h"

#define PI 3.1415926535

//runs the simulation for the required steps
void run_simulation()
{

global.total_time = 100000;
// global.total_time = 40000;
// global.echo_time = 5000;
global.echo_time = 1000;
global.movie_time = 100;
global.N_particles = 400;
global.SX = 60;
global.SY = 60;
global.halfSX = 30;
global.halfSY = 30;
// get_starting_time();

init_particles();
init_pinningsites();
for(global.time=0;global.time<global.total_time;global.time++)
    {
    calculate_external_forces_on_particles();
    calculate_pairwise_forces();
    calculate_pinning_forces();
   
    move_particles();
        
    //echo time
    if (global.time % global.echo_time == 0)
        {
        printf("Timestep: %d / %d\n",global.time,global.total_time);
        fflush(stdout);
        }
    
    // movie write time
    if (global.time % global.movie_time == 0)
        write_cmovie_frame();
        
    }

    global.particle_driving_force += 0.01;
    
// get_finishing_time();
// echo_running_time();
}


void run_verlet_simulation_slow()
{
global.total_time = 100000;
global.echo_time = 1000;
global.movie_time = 100;
global.N_particles = 500;
global.SX = 60;
global.SY = 60;
global.halfSX = 30;
global.halfSY = 30;

init_particles();
init_pinningsites();

rebuild_pinning_grid();
rebuild_Verlet_list();//build for the first time

for(global.time=0;global.time<global.total_time;global.time++)
    {
    calculate_external_forces_on_particles();
    calculate_pairwise_forces_verlet();
    calculate_pinning_forces_slow();
    move_particles();
    
    check_Verlet_rebuild_condition_and_set_flag();
    
    //one particle moved enough to rebuild the Verlet list
    if (global.flag_to_rebuild_Verlet==1)
        rebuild_Verlet_list();

    //echo time
    if (global.time % global.echo_time == 0)
        {
        printf("Timestep: %d / %d\n",global.time,global.total_time);
        fflush(stdout);
        }
    
    //movie write time
    if (global.time % global.movie_time == 0)
        write_cmovie_frame();

    }
}
void run_verlet_simulation()
{
global.total_time = 1000000;
global.echo_time = 1000;
global.movie_time = 1000;
global.N_particles = 500;
global.SX = 60;
global.SY = 60;
global.halfSX = 30;
global.halfSY = 30;

init_particles();
init_pinningsites();
rebuild_pinning_grid();
rebuild_Verlet_list();//build for the first time

for(global.time=0;global.time<global.total_time;global.time++)
    {
    calculate_external_forces_on_particles();
    calculate_pairwise_forces_verlet();
    calculate_pinning_forces();
    move_particles();
    
    check_Verlet_rebuild_condition_and_set_flag();
    
    //one particle moved enough to rebuild the Verlet list
    if (global.flag_to_rebuild_Verlet==1)
        rebuild_Verlet_list();

    //echo time
    if (global.time % global.echo_time == 0)
        {
        printf("Timestep: %d / %d\n",global.time,global.total_time);
        fflush(stdout);
        }
    
    //movie write time
    if (global.time % global.movie_time == 0) {
        write_cmovie_frame();
        write_statistics();
    }
    global.particle_driving_force += 0.005;
    }
}

void calculate_pinning_forces_slow() {
    int i, j, pi, pj;   // i = particle, j = pinning_site
    double r,r2,f;
    double dx,dy;

    for(i=0;i<global.N_particles;i++) {
        for(j=0;j<global.N_pinningsites;j++) {
            distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
                                        global.pinningsite_x[j],global.pinningsite_y[j],&r2,&dx,&dy);
        
            r = sqrt(r2);
                if(r < global.pinningsite_R[j]) {
                    f = global.pinningsite_fmax[j] * (r / global.pinningsite_R[j]);
                    global.particle_fx[i] = f * (dx / global.pinningsite_R[j]);
                    global.particle_fy[i] = f * (dy / global.pinningsite_R[j]);
                }
        }
    }  
}

void calculate_pinning_forces() {
    int i, j, pi, pj;   // i = particle, j = pinning_site
    double r,r2,f;
    double dx,dy;
    
    for(i=0;i<global.N_particles;i++) {
        calculate_neighbour_pinning_sites(i);
    }        
}

void update_particle_force_with_pinning(int i, int psi)
{
    double r,r2,f;
    double dx,dy;
    distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
                                        global.pinningsite_x[psi],global.pinningsite_y[psi],&r2,&dx,&dy);
        
    r = sqrt(r2);
    if(r < global.pinningsite_R[psi]) {
                f = global.pinningsite_fmax[psi] * (r / global.pinningsite_R[psi]);
                global.particle_fx[i] += f * (dx / global.pinningsite_R[psi]);
                global.particle_fy[i] += f * (dy / global.pinningsite_R[psi]);
    }
}

void calculate_neighbour_pinning_sites(int i)
{
    int pi, pj;


    pi = (int) global.particle_x[i]/global.pinningsite_grid_dx;
    pj = (int) global.particle_y[i]/global.pinningsite_grid_dy;
    
    int psi = global.pinningsite_grid[pi][pj];
    update_particle_force_with_pinning(i, psi);

    //above
    int pi_up = pi - 1;
    if (pi_up < 0 ) pi_up += global.pinningsite_grid_dx;
    psi = global.pinningsite_grid[pi_up][pj];
    update_particle_force_with_pinning(i, psi);

    //below
    int pi_down = pi + 1;
    if (pi_down >= global.pinningsite_grid_dx) pi_down -= global.pinningsite_grid_dx;
    psi = global.pinningsite_grid[pi_down][pj];
    update_particle_force_with_pinning(i, psi);

    //right
    int py_right = pj + 1;
    if (py_right >= global.pinningsite_grid_dy) py_right -= global.pinningsite_grid_dy;
    psi = global.pinningsite_grid[pi][py_right];
    update_particle_force_with_pinning(i, psi);

    //left
    int py_left = pj - 1;
    if (py_right < 0) py_left += global.pinningsite_grid_dy;
    psi = global.pinningsite_grid[pi][py_left];
    update_particle_force_with_pinning(i, psi);

    //right up
    psi = global.pinningsite_grid[pi_up][py_right];
    update_particle_force_with_pinning(i, psi);

    //left up
    psi = global.pinningsite_grid[pi_up][py_left];
    update_particle_force_with_pinning(i, psi);

    //right down
    psi = global.pinningsite_grid[pi_down][py_right];
    update_particle_force_with_pinning(i, psi);

    //left down
    psi = global.pinningsite_grid[pi_down][py_left];
    update_particle_force_with_pinning(i, psi);
}


void calculate_external_forces_on_particles()
{
int i;

for(i=0;i<global.N_particles;i++)
    {
        //increase particle driving force slowly with time
    global.particle_fx[i] += global.particle_direction[i] * global.particle_driving_force;
    }
}

void calculate_pairwise_forces()
{
int i,j;
double r,r2,f;
double dx,dy;


for(i=0;i<global.N_particles;i++)
    for(j=i+1;j<global.N_particles;j++)
    {
    //calculate the distance between the particles
    distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&r2,&dx,&dy);
    
    //the particles are close enough to interact
    if (r2<16.0)
        {

        r = sqrt(r2);
        if (r<0.2)
            {
            // printf("WARNING:PARTICLES TOO CLOSE. LOWER CUTOFF FORCE USED\n");
            f = 100.0;
            }
        else
            {
            //calculate the force
            f = 1/r2 * exp(-r / global.particle_screening_length);
            }
            
        //projection to the x,y axes
        f = f/r;
        
        global.particle_fx[i] -= f*dx;
        global.particle_fy[i] -= f*dy;
    
        global.particle_fx[j] += f*dx;
        global.particle_fy[j] += f*dy;
        }
   
    }
}

void calculate_pairwise_forces_verlet()
{
int ii,i,j;
double r,r2,f;
double dx,dy;

for(ii=0;ii<global.N_Verlet;ii++)
    {
    //obtain the i,j from the Verlet list
    i = global.Verletlisti[ii];
    j = global.Verletlistj[ii];
    
    //perform the pairwise force calculation
    distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&r2,&dx,&dy);
    
    //non-tabulated version
    //try to not divide just multiply division is costly
    r = sqrt(r2);
    if (r<0.2)
        {
        // printf("WARNING:PARTICLES TOO CLOSE. LOWER FORCE USED %d - %d \n", i, j);
        f = 100.0;
        }
    else
        {
        f = 1/r2 * exp(-r * global.particle_screening_wavevector);
        // if (r<1.0) printf("r=%lf f=%lf\n",r,f);
        }
        
    //division and multiplication for projection to the x,y axes
    
    f = f/r;
    //if (r<1.0) printf("%f/r=lf fx=%lf fy=%lf\n\n",f,f*dx,f*dy);
    
    global.particle_fx[i] -= f*dx;
    global.particle_fy[i] -= f*dy;
    
    global.particle_fx[j] += f*dx;
    global.particle_fy[j] += f*dy;
    }

}

//calculates the shortest distance squared between 2 points in a PBC configuration
//this is squared because I want to save on sqrt with the lookup table
//also used by the Verlet rebuild flag check where I check how much a particle moved
void distance_squared_folded_PBC(double x0,double y0,double x1,double y1, double *r2_return, double *dx_return,double *dy_return)
{
double dr2;
double dx,dy;

dx = x1 - x0;
dy = y1 - y0;

//PBC fold back
//if any distance is larger than half the box
//the copy in the neighboring box is closer
if (dx > global.halfSX) dx -= global.SX;
if (dx <= -global.halfSX) dx += global.SX;
if (dy > global.halfSY) dy -= global.SY;
if (dy <= -global.halfSY) dy += global.SY;

dr2 = dx*dx + dy*dy;

*r2_return = dr2;
*dx_return = dx;
*dy_return = dy;
}


//moves the particles one time step
void move_particles()
{
int i;
double dx,dy;

for(i=0;i<global.N_particles;i++)
    {
    dx = global.particle_fx[i] * global.dt;
    dy = global.particle_fy[i] * global.dt;
    
    global.particle_x[i] += dx;
    global.particle_y[i] += dy;

    global.particle_dx_so_far[i] += dx;
    global.particle_dy_so_far[i] += dy;
    
    //if (i==300) printf("%lf\n",global.particle_x[i]);
  
    fold_particle_back_PBC(i);
    
    global.particle_fx[i] = 0.0;
    global.particle_fy[i] = 0.0;
    }
}

void fold_particle_back_PBC(int i)
{

//fold back the particle into the pBC simulation box
//assumes it did not jump more thana  box length
//if it did the simulation is already broken anyhow

if (global.particle_x[i]<0) global.particle_x[i] += global.SX;
if (global.particle_y[i]<0) global.particle_y[i] += global.SY;
if (global.particle_x[i]>=global.SX) global.particle_x[i] -= global.SX;
if (global.particle_y[i]>=global.SY) global.particle_y[i] -= global.SY;
}

void check_Verlet_rebuild_condition_and_set_flag()
{
int i;
double dr2;

global.flag_to_rebuild_Verlet = 0;

for(i=0;i<global.N_particles;i++)
        {
        dr2 = global.particle_dx_so_far[i] * global.particle_dx_so_far[i] +
                global.particle_dy_so_far[i] * global.particle_dy_so_far[i];

        if (dr2>=global.Verlet_intershell_squared)
            {
            global.flag_to_rebuild_Verlet = 1;
            break; //exit the for cycle
            }
        }

}

//rebuilds the Verlet list
void rebuild_Verlet_list()
{
int i,j;
double dr2,dx,dy;
double estimation;

//initialize the Verlet list for the first time
if (global.N_Verlet_max==0)
    {
    //we are building the Verlet list for the first time in the simulation
    estimation = global.N_particles/(double)global.SX/(double)global.SY;
    estimation *= PI * global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
    global.N_Verlet_max = (int)estimation * global.N_particles / 2;
    global.Verletlisti = (int *) malloc(global.N_Verlet_max * sizeof(int));
    global.Verletlistj = (int *) malloc(global.N_Verlet_max * sizeof(int));
    }

//build the Verlet list
global.N_Verlet = 0;

for(i=0;i<global.N_particles;i++)
    {
    for(j=i+1;j<global.N_particles;j++)
        {
        distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&dr2,&dx,&dy);
            
        if (dr2<36.0)
            {
            global.Verletlisti[global.N_Verlet] = i;
            global.Verletlistj[global.N_Verlet] = j;
            
            global.N_Verlet++;
            if (global.N_Verlet>=global.N_Verlet_max)
                {
                global.N_Verlet_max = (int)(1.1*global.N_Verlet);
                global.Verletlisti = (int *) realloc(global.Verletlisti ,global.N_Verlet_max * sizeof(int));
                global.Verletlistj = (int *) realloc(global.Verletlistj ,global.N_Verlet_max * sizeof(int));
                }
            }
        }
    }

global.flag_to_rebuild_Verlet = 0;
for(i=0;i<global.N_particles;i++)
    {
    global.particle_dist_so_far[i] += global.particle_dx_so_far[i];
    global.particle_dx_so_far[i] = 0.0;
    global.particle_dy_so_far[i] = 0.0;
    }
}


void rebuild_pinning_grid()
{
int i,j;
int gi,gj;

if (global.pinningsite_grid==NULL)
    {
    //build the pinningsite grid for the first time;
    global.Nx_pinningsite_grid = (int) (global.SX/global.pinningsite_setradius) + 1;
    global.Ny_pinningsite_grid = (int) (global.SY/global.pinningsite_setradius) + 1;
    printf("Pinning sites grid is %d x %d\n",global.Nx_pinningsite_grid,global.Ny_pinningsite_grid);
    global.pinningsite_grid_dx = global.SX/global.Nx_pinningsite_grid;
    global.pinningsite_grid_dy = global.SX/global.Ny_pinningsite_grid;
    printf("Pinning sites cell is %.2lf x %.2lf\n", global.pinningsite_grid_dx, global.pinningsite_grid_dy);
    printf("Pinning sites radius is = %.2lf\n",global.pinningsite_setradius);
    
    //initialize the grid
    global.pinningsite_grid = (int **) malloc(global.Nx_pinningsite_grid * sizeof(int *));
    for(i=0;i<global.Nx_pinningsite_grid;i++)
        global.pinningsite_grid[i] = (int *) malloc(global.Ny_pinningsite_grid * sizeof(int));
    }
    
//always do this - zero the values
for(i=0;i<global.Nx_pinningsite_grid;i++)
    for(j=0;j<global.Ny_pinningsite_grid;j++)
        global.pinningsite_grid[i][j] = -1;
    
//always do this - fill up the values
for(i=0;i<global.N_pinningsites;i++)
    {
    gi = (int) global.pinningsite_x[i]/global.pinningsite_grid_dx;
    gj = (int) global.pinningsite_y[i]/global.pinningsite_grid_dy;
    //this cannot fail (hopefully)
    global.pinningsite_grid[gi][gj] = i;
    }

}


void write_cmovie_frame()
{
int i;
float floatholder;
int intholder;

intholder = global.N_particles + global.N_pinningsites;
fwrite(&intholder,sizeof(int),1,global.moviefile);

intholder = global.time;
fwrite(&intholder,sizeof(int),1,global.moviefile);

for (i=0;i<global.N_pinningsites;i++)
    {
    intholder = global.pinningsite_color[i];
    fwrite(&intholder,sizeof(int),1,global.moviefile);
    intholder = i;//ID
    fwrite(&intholder,sizeof(int),1,global.moviefile);
    floatholder = (float)global.pinningsite_x[i];
    fwrite(&floatholder,sizeof(float),1, global.moviefile);
    floatholder = (float)global.pinningsite_y[i];
    fwrite(&floatholder,sizeof(float),1, global.moviefile);
    floatholder = global.pinningsite_R[i];//cum_disp, cmovie format
    fwrite(&floatholder,sizeof(float),1,global.moviefile);
    }

for (i=0;i<global.N_particles;i++)
    {
    intholder = global.particle_color[i];
    fwrite(&intholder,sizeof(int),1,global.moviefile);
    intholder = i;//ID
    fwrite(&intholder,sizeof(int),1,global.moviefile);
    floatholder = (float)global.particle_x[i];
    fwrite(&floatholder,sizeof(float),1, global.moviefile);
    floatholder = (float)global.particle_y[i];
    fwrite(&floatholder,sizeof(float),1, global.moviefile);
    floatholder = 1.0;//cum_disp, cmovie format
    fwrite(&floatholder,sizeof(float),1,global.moviefile);
    }
}

void write_statistics() {
    int i;
    double distance = 0.0;

    for(i=0;i<global.N_particles;i++) {
        distance += ((-1) * global.particle_dist_so_far[i]);
    }

    distance /= global.N_particles;
    distance /= global.time;

    FILE *testf;
    testf = fopen("test/avg_velocity.csv","a");
    fprintf(testf,"%f,%d\n",distance, global.time);
    fclose(testf);
}
