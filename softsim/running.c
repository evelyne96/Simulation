//
//  running.c
//  softsim
//
//  Created by András Libál on 7/19/18.
//  Copyright © 2018 András Libál. All rights reserved.
//
#include "running.h"
#include "globaldata.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "list.h"

#define PI 3.1415926535

//runs the simulation for the required steps
void run_simulation()
{

global.total_time = 100000;
global.echo_time = 1000;
global.movie_time = 100;

rebuild_Verlet_list();//build for the first time

for(global.time=0;global.time<global.total_time;global.time++)
    {
    calculate_external_forces_on_particles();
    calculate_pairwise_forces();
    move_particles();
    
    check_Verlet_rebuild_condition_and_set_flag();
    
    //one particle moved enough to rebuild the Verlet list
    if (global.flag_to_rebuild_Verlet==1)
        rebuild_Verlet_list();
    // //movie write time
    // if (global.time % global.movie_time == 0) {
    //     test_program_by_coloring();
    //     write_cmovie_frame();
    // }
    }
}

void run_simulation_Verlet_lookup_cell()
{
global.total_time = 100000;
global.echo_time = 1000;
global.movie_time = 100;

rebuild_Verlet_cell();//build for the first time
// rebuild_Verlet_list_with_cell();
rebuild_Verlet_list();

for(global.time=0;global.time<global.total_time;global.time++)
    {
    calculate_external_forces_on_particles();
    calculate_pairwise_forces();
    move_particles();
    
    check_Verlet_rebuild_condition_and_set_flag();
    
    //one particle moved enough to rebuild the Verlet list
    if (global.flag_to_rebuild_Verlet==1)
    {
        // rebuild_Verlet_cell();
        // rebuild_Verlet_list_with_cell();
        rebuild_Verlet_list();
    }

    // //movie write time
    // if (global.time % global.movie_time == 0)
    //     write_cmovie_frame();
    // }
    }
}

void run_simulation_simple()
{

global.total_time = 100000;
global.echo_time = 1000;
global.movie_time = 100;

for(global.time=0;global.time<global.total_time;global.time++)
    {
    
    calculate_external_forces_on_particles();
    calculate_pairwise_forces_simple();
    move_particles();
    }
}

//change the direction of the pinning sites
//they will move on a circle
void rotate_pinningsite_directions()
{
double theta,sint,cost;
double newx,newy;
int i;

//theta  = 3.1415/180.0*0.005;
theta = global.pinning_driving_force * global.dt * 0.2857;
sint = sin(theta);
cost = cos(theta);

//printf("%lf %lf\n",sint,cost);

//rotation matrix
//cos(theta)*x -sin(theta)*y
//sin(theta)*x +cos(theta)*y

for(i=0;i<global.N_particles;i++)
    {
     newx = global.particle_direction_x[i] * cost - global.particle_direction_y[i] * sint;
     newy = global.particle_direction_x[i] * sint + global.particle_direction_y[i] * cost;
    
    global.particle_direction_x[i] = newx;
    global.particle_direction_y[i] = newy;
    }

//printf("%lf %lf\n",global.particle_direction_x[i],global.particle_direction_y[i]);

}


//change the direction of the pinning sites
//they will move on a triangle
void adjust_pinningsite_directions()
{
int i,k;

if ( (global.pinningsite_dx_so_far[0]*global.pinningsite_dx_so_far[0] +
      global.pinningsite_dy_so_far[0]*global.pinningsite_dy_so_far[0])
      >= global.pinning_lattice_constant * global.pinning_lattice_constant * 1.0)
      {
      //change directions now
      
      for(k=0;k<global.N_pinningsites;k++)
        switch ( (global.pinningsite_color[k]-2 + global.pinning_direction_change)%3 )
            {
            case 0: {
                    global.pinningsite_direction_x[k] = - 0.5;
                    global.pinningsite_direction_y[k] = - sqrt(3)/2.0;
                    break;
                    }
            case 1: {
                    global.pinningsite_direction_x[k] = 1.0;
                    global.pinningsite_direction_y[k] = 0.0;
                    break;
                    }
            case 2: {
                    global.pinningsite_direction_x[k] = - 0.5;
                    global.pinningsite_direction_y[k] = + sqrt(3)/2.0;
                    break;
                    }
            }
      global.pinning_direction_change++;
      for(i=0;i<global.N_pinningsites;i++)
        {
        global.pinningsite_dx_so_far[i] = 0.0;
        global.pinningsite_dy_so_far[i] = 0.0;
        }
      }
}

void calculate_external_forces_on_pinningsites()
{
int i;

for(i=0;i<global.N_pinningsites;i++)
    {
    global.pinningsite_fx[i] += global.pinningsite_direction_x[i] * global. pinning_driving_force;
    global.pinningsite_fy[i] += global.pinningsite_direction_y[i] * global. pinning_driving_force;
    }
}



//moves the pinning sites one time step
void move_pinningsites()
{

int i;
double dx,dy;

for(i=0;i<global.N_pinningsites;i++)
    {
    dx = global.pinningsite_fx[i] * global.dt;
    dy = global.pinningsite_fy[i] * global.dt;
    
    global.pinningsite_x[i] += dx;
    global.pinningsite_y[i] += dy;
    
    global.pinningsite_dx_so_far[i] += dx;
    global.pinningsite_dy_so_far[i] += dy;
    
    fold_pinningsite_back_PBC(i);
    
    global.pinningsite_fx[i] = 0.0;
    global.pinningsite_fy[i] = 0.0;
    }
}

void fold_pinningsite_back_PBC(int i)
{

//fold back the pinningsite into the pBC simulation box
//assumes it did not jump more thana  box length
//if it did the simulation is already broken anyhow

if (global.pinningsite_x[i]<0) global.pinningsite_x[i] += global.SX;
if (global.pinningsite_y[i]<0) global.pinningsite_y[i] += global.SY;
if (global.pinningsite_x[i]>=global.SX) global.pinningsite_x[i] -= global.SX;
if (global.pinningsite_y[i]>=global.SY) global.pinningsite_y[i] -= global.SY;
}

void calculate_external_forces_on_particles()
{
int i;

for(i=0;i<global.N_particles;i++)
    {
    global.particle_fx[i] += global.particle_direction_x[i] * global. particle_driving_force;
    //global.particle_fy[i] += global.particle_direction_y[i] * global. particle_driving_force;
    //printf("ext_fx= %lf\n",global.particle_direction_x[i] * global. particle_driving_force);
    }
}

void calculate_pairwise_forces()
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
        f = 1/r2 * exp(-r * global.partile_particle_screening_wavevector);
        //if (r<1.0) printf("r=%lf f=%lf\n",r,f);
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

void calculate_pairwise_forces_simple()
{
int i,j;
double r,r2,f;
double dx,dy;

for(i=0;i<global.N_particles;i++)
    for(j=i+1;j<global.N_particles;j++)
    {
    
    //perform the pairwise force calculation
    distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&r2,&dx,&dy);
    
    //non-tabulated version
    //try to not divide just multiply division is costly
    r = sqrt(r2);
    if (r<0.2)
        {
        printf("WARNING:PARTICLES TOO CLOSE. LOWER FORCE USED\n");
        f = 100.0;
        }
    else
        {
        f = 1/r2 * exp(-r * global.partile_particle_screening_wavevector);
        //if (r<1.0) printf("r=%lf f=%lf\n",r,f);
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

void calculate_pairwise_forces_Verlet_lookup_cell()
{
int i,j;

for(i=0;i<global.Verlet_cell_list_size;i++)
{
    int up = i - global.Ny_Verlet_cell_list;
    int right = i + 1;
    if (up < 0) up += global.Verlet_cell_list_size;
    if (right == global.Verlet_cell_list_size) right = 0;
    int up_right = up + 1;
    int down_right = 0;
    if ((i+1) % global.Ny_Verlet_cell_list == 0) {
        up_right -= global.Ny_Verlet_cell_list;
        down_right = i + 1;
    }

    if (i >= (global.Verlet_cell_list_size - global.Ny_Verlet_cell_list)) {
        down_right = i - (global.Nx_Verlet_cell_list - 1) * global.Ny_Verlet_cell_list + 1;
        if (i == (global.Verlet_cell_list_size - 1)) down_right = 0;
    }

    go_through_linked_Verlet_list(i, i, 1);
    go_through_linked_Verlet_list(i, up_right, 0);
    go_through_linked_Verlet_list(i, right, 0);
    go_through_linked_Verlet_list(i, up, 0);
    go_through_linked_Verlet_list(i, down_right, 0);
}
}

void go_through_linked_Verlet_list(int one, int two, int same)
{
  List *current = global.Verlet_cell_list[one];
  List *other = global.Verlet_cell_list[two];
  if(current->head == NULL || other->head == NULL)  return;

  Node *start_node = current->head;
  Node *start_other = other->head;

  if (same)
  {
    start_other = start_other->next;
  }

  Node *other_start_node = start_other;
  
  
  for(; start_node != NULL; start_node = start_node->next) 
  {
      for(; other_start_node != NULL; other_start_node = other_start_node->next)
      {
         int i = start_node->data;
         int j = other_start_node->data;
        //  printf("i:%d,j:%d\n",i,j);
         if (i == j) continue;
         calculate_forces_between(i, j);
      }
      if (same)
      {
           start_other = start_other->next;
      }
      other_start_node = start_other;
  }
}

void update_verlet_list_from_cell(int one, int two, int same)
{
  List *current = global.Verlet_cell_list[one];
  List *other = global.Verlet_cell_list[two];
  if(current->head == NULL || other->head == NULL)  return;

  Node *start_node = current->head;
  Node *start_other = other->head;

  if (same)
  {
    start_other = start_other->next;
  }

  Node *other_start_node = start_other;
  double dr2,dx,dy;

  for(; start_node != NULL; start_node = start_node->next) 
  {
      for(; other_start_node != NULL; other_start_node = other_start_node->next)
      {
         int i = start_node->data;
         int j = other_start_node->data;
        //  printf("i:%d,j:%d\n",i,j);
         if (i == j) continue;
          distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&dr2,&dx,&dy);
                if (dr2<36.0)
            {
            global.Verletlisti[global.N_Verlet] = i;
            global.Verletlistj[global.N_Verlet] = j;
            
            global.N_Verlet++;
            if (global.N_Verlet>=global.N_Verlet_max)
                {
                printf("Verlet list reallocated from %d\n",global.N_Verlet_max);
                global.N_Verlet_max = (int)(1.1*global.N_Verlet);
                global.Verletlisti = (int *) realloc(global.Verletlisti ,global.N_Verlet_max * sizeof(int));
                global.Verletlistj = (int *) realloc(global.Verletlistj ,global.N_Verlet_max * sizeof(int));
                printf("New Verlet list max size = %d\n",global.N_Verlet_max);
                }
            }
      }
      other_start_node = start_other;
   }
}

void calculate_forces_between(int i, int j)
{
    double r,r2,f;
    double dx,dy;
        //perform the pairwise force calculation
    distance_squared_folded_PBC(global.particle_x[i],global.particle_y[i],
            global.particle_x[j],global.particle_y[j],&r2,&dx,&dy);
    
    //non-tabulated version
    //try to not divide just multiply division is costly
    r = sqrt(r2);
    if (r<0.2)
        {
        printf("WARNING:PARTICLES TOO CLOSE. LOWER FORCE USED %.f %.f\n", r2, r);
        f = 100.0;
        }
    else
        {
        f = 1/r2 * exp(-r * global.partile_particle_screening_wavevector);
        //if (r<1.0) printf("r=%lf f=%lf\n",r,f);
        }
        
    //division and multiplication for projection to the x,y axes
    
    f = f/r;
    //if (r<1.0) printf("%f/r=lf fx=%lf fy=%lf\n\n",f,f*dx,f*dy);
    
    global.particle_fx[i] -= f*dx;
    global.particle_fy[i] -= f*dy;
    
    global.particle_fx[j] += f*dx;
    global.particle_fy[j] += f*dy;
}


void check_Verlet_rebuild_condition_and_set_flag()
{
int i;
double dr2;

//check if any particle moved (potentially) enough to enter the inner Verlet shell
//coming from the outside, if one did, rebuild the Verlet lists
global.flag_to_rebuild_Verlet = 0;

for(i=0;i<global.N_particles;i++)
        {
        dr2 = global.particle_dx_so_far[i] * global.particle_dx_so_far[i] +
                global.particle_dy_so_far[i] * global.particle_dy_so_far[i];
        
        //if (i==0) printf("%d %lf\n",global.time,dr2);fflush(stdout);
        
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
    global.particle_dx_so_far[i] = 0.0;
    global.particle_dy_so_far[i] = 0.0;
    }
}

void rebuild_Verlet_list_with_cell()
{
int i,j;
double dr2,dx,dy;
double estimation;

//initialize the Verlet list for the first time
if (global.N_Verlet_max==0)
    {
    estimation = global.N_particles/(double)global.SX/(double)global.SY;
    // printf("System density is %.3lf\n",estimation);
    
    estimation *= PI * global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
    // printf("Particles in a R = %.2lf shell = %lf\n", global.Verlet_cutoff_distance,estimation);
    
    global.N_Verlet_max = (int)estimation * global.N_particles / 2;
    // printf("Estimated N_Verlet_max = %d\n",global.N_Verlet_max);
    
    global.Verletlisti = (int *) malloc(global.N_Verlet_max * sizeof(int));
    global.Verletlistj = (int *) malloc(global.N_Verlet_max * sizeof(int));
    }

//build the Verlet list
global.N_Verlet = 0;
for(i=0;i<global.Verlet_cell_list_size;i++)
{
    int up = i - global.Ny_Verlet_cell_list;
    int right = i + 1;
    if (up < 0) up += global.Verlet_cell_list_size;
    if (right == global.Verlet_cell_list_size) right = 0;
    int up_right = up + 1;
    int down_right = 0;
    if ((i+1) % global.Ny_Verlet_cell_list == 0) {
        up_right -= global.Ny_Verlet_cell_list;
        down_right = i + 1;
    }

    if (i >= (global.Verlet_cell_list_size - global.Ny_Verlet_cell_list)) {
        down_right = i - (global.Nx_Verlet_cell_list - 1) * global.Ny_Verlet_cell_list + 1;
        if (i == (global.Verlet_cell_list_size - 1)) down_right = 0;
    }
    update_verlet_list_from_cell(i, i, 1);
    update_verlet_list_from_cell(i, up, 0);
    update_verlet_list_from_cell(i, right, 0);
    update_verlet_list_from_cell(i, up_right, 0);
    update_verlet_list_from_cell(i, down_right, 0);

}

global.flag_to_rebuild_Verlet = 0;
for(i=0;i<global.N_particles;i++)
    {
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
    
    /*if ((gi>=global.Nx_pinningsite_grid)||(gj>=global.Ny_pinningsite_grid))
        {
        printf("Out of bounds pinningsite\n");
        exit(1);
        }*/
    //this cannot fail (hopefully)
    global.pinningsite_grid[gi][gj] = i;
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
    
    fold_particle_back_PBC(i);

    //  int gridx = ((int)floor(global.particle_x[i]) / global.Nx_Verlet_cell_list);
    //  int gridy = ((int)floor(global.particle_y[i]) / global.Ny_Verlet_cell_list);
    //  int pos = (gridx * global.Ny_Verlet_cell_list) + gridy;
    //  global.particle_color[i] = (gridy+gridx) %4;
    
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

void write_cmovie_frame()
{
int i;
float floatholder;
int intholder;

//implement the color testing
test_program_by_coloring();

//legacy cmovie format for del-plot

intholder = global.N_pinningsites + global.N_particles;
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

//TESTS


void test_program_by_coloring()
{
int ii;
int i,j;

//testing the Verlet list by coloring one particle's neightbors one color
/*
for(i=0;i<global.N_particles;i++)
    global.particle_color[i] = 2;

for(ii=0;ii<global.N_Verlet;ii++)
    {
    i = global.Verletlisti[ii];
    j = global.Verletlistj[ii];
    if (i==150) global.particle_color[j] = 3;
    if (j==150) global.particle_color[i] = 3;
    }
//printf("\n");
*/

//testing the pinning grid by coloring the pins
for(i=0;i<global.Nx_pinningsite_grid;i++)
    for(j=0;j<global.Ny_pinningsite_grid;j++)
        {
        if (global.pinningsite_grid[i][j]!=-1)
            {
            global.pinningsite_color[global.pinningsite_grid[i][j]] = i%2 + 2;
            //printf("%d %d (%d)-> %d\n",i,j,global.pinningsite_grid[i][j], i%2 + 2);
            //printf("%d\n",global.pinningsite_color[global.pinningsite_grid[i][j]]);
            }
        }

}

void rebuild_Verlet_cell() 
{
    if (global.Verlet_cell_list==NULL)
    {
        //build the Verlet cell grid for the first time;
        global.Nx_Verlet_cell_list = (int) floor(global.SX/global.Verlet_cutoff_distance) ;
        global.Ny_Verlet_cell_list = (int) floor(global.SY/global.Verlet_cutoff_distance) ;
        printf("With floor Nx:%d,Ny:%d\n",global.Nx_Verlet_cell_list, global.Ny_Verlet_cell_list);
        global.Verlet_cell_list_size = global.Nx_Verlet_cell_list * global.Ny_Verlet_cell_list; 
        global.Verlet_cell_list = malloc(global.Verlet_cell_list_size * sizeof(List));
    }  

    for(int i=0;i<global.Verlet_cell_list_size;i++)
    {
        destroy(global.Verlet_cell_list[i]);
        global.Verlet_cell_list[i] = makelist();
    }

    for(int i=0;i<global.N_particles;i++)
    {
        int gridx = (int) floor(global.particle_x[i] / global.Nx_Verlet_cell_list);
        int gridy = (int) floor(global.particle_y[i] / global.Ny_Verlet_cell_list);
        int pos = (gridx * global.Ny_Verlet_cell_list) + gridy;
        add(i, global.Verlet_cell_list[pos]);
    }
}

//TESTS FOR COMPARISION
void run_test_verlet_vs_time()
{
    FILE *testf = fopen("test/verlet_vs_time.csv","wt");

    global.SX = 42.42;
    global.SY = global.SX;
    for(global.N_particles=200;global.N_particles<400;global.N_particles+=100)
    {
        global.halfSX = global.SX / 2.0;
        global.halfSY = global.SY / 2.0;
        printf("%d\n",global.N_particles);
        clock_t before = clock();
        // run_simulation_simple();
        run_simulation_Verlet_lookup_cell();
        // run_simulation();
        clock_t after = clock();
        double time_spent = (double)(after - before) / CLOCKS_PER_SEC;
        fprintf(testf,"%d,%f\n",global.N_particles,time_spent);
        printf("After simulation\n");
        global.SX = round( sqrt(((global.N_particles+100)*global.SX*global.SY)/global.N_particles));
        global.SY = global.SX;
    }
    fclose(testf);
}
