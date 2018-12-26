//  initializer.c
//  Simulation Methods course, 2018
//  First Assignment: Molecular Dynamics (Brownian Dynamics) Simulation

#include "initializer.h"
#include "globaldata.h"
#include <math.h>

//init the common variables in the simulation
void init_simulation()
{
//timestep length
global.dt = 0.001;
printf("Timestep (dt) = %lf\n",global.dt);
global.particle_screening_length = 4.0;
printf("Screening length = %lf\n",global.particle_screening_length);
global.particle_driving_force = 0.5;
printf("Driving force on particles = %lf\n",global.particle_driving_force);
global.particle_screening_wavevector =  1.0/global.particle_screening_length;
printf("Particle particle screening wavevector = %lf\n",global.particle_screening_wavevector);

//zero everything so rebuild Verlet can find this the first time
global.Verletlisti = NULL;
global.Verletlistj = NULL;
global.N_Verlet = 0;
global.N_Verlet_max = 0;

//zero verlet cell list
global.Verlet_cell_list = NULL;
global.Nx_Verlet_cell_list = 0;
global.Ny_Verlet_cell_list = 0;

global.Verlet_cutoff_distance = 1.5 * global.particle_screening_length;
global.Verlet_cutoff_distance_squared = global.Verlet_cutoff_distance * global.Verlet_cutoff_distance;
global.Verlet_intershell_squared = global.Verlet_cutoff_distance - global.particle_screening_length;
global.Verlet_intershell_squared = global.Verlet_intershell_squared / 2.0;
global.Verlet_intershell_squared *= global.Verlet_intershell_squared;

printf("Verlet cutoff distance = %.2lf\n",global.Verlet_cutoff_distance);
printf("Verlet cutoff distance squared = %.2lf\n",global.Verlet_cutoff_distance_squared);
printf("Half of Verlet intershell distance = %.2lf\n",sqrt(global.Verlet_intershell_squared));
printf("Half of Verlet intershell distance squared = %.2lf\n",global.Verlet_intershell_squared);


//zero everythiong so rebuild_pinning_grid can find this the first time
global.pinningsite_grid = NULL;
global.Nx_pinningsite_grid = 0;
global.Ny_pinningsite_grid = 0;


//the dirving force for the pinning sites needs to be set
global.pinning_driving_force = 0.5;
printf("Pinning site driving force = %lf\n",global.pinning_driving_force);

//radius of pinning sites
global.pinningsite_setradius = 0.5;
printf("Pinning site radius set to = %.2lf\n",global.pinningsite_setradius);

//the dirving force for the pinning sites needs to be set
global.particle_driving_force = 0.5;
printf("Particle driving force = %lf\n",global.pinning_driving_force);

}

//init the simulation box
void init_simulation_box()
{
printf("Initializing the simulation box\n");
// global.SX = 60.0;
// global.SY = 60.0;
global.SX = 41.0;
global.SY = 41.0;
global.halfSX = global.SX/2.0;
global.halfSY = global.SY/2.0;

printf("SX = %lf SY = %lf\n",global.SX,global.SY);
}


void init_particles()
{
// global.N_particles = 800;


global.particle_x = (double *)malloc(global.N_particles*sizeof(double));
global.particle_y = (double *)malloc(global.N_particles*sizeof(double));
global.particle_fx = (double *)malloc(global.N_particles*sizeof(double));
global.particle_fy = (double *)malloc(global.N_particles*sizeof(double));
global.particle_color = (int *)malloc(global.N_particles*sizeof(int));

global.particle_dx_so_far = (double *)malloc(global.N_particles*sizeof(double));
global.particle_dy_so_far = (double *)malloc(global.N_particles*sizeof(double));

global.particle_dist_so_far = (double *)malloc(global.N_particles*sizeof(double));

global.particle_direction = (double *)malloc(global.N_particles*sizeof(double));

init_particles_randomly();

printf("Particles initialized\n");
printf("N_particles = %d\n",global.N_particles);
}

//calculates the shortest distance between 2 points in a PBC configuration
//only used in the random deposition chekc which happens once at initialization
//so this is not so time crucial left the square root inside
double distance_folded_PBC(double x0,double y0,double x1,double y1)
{
double r;
double dx,dy;

dx = x1 - x0;
dy = y1 - y0;

//PBC fold back
//if any distance is lrger than half the box
//the copy in the neighboring box is closer
if (dx > global.halfSX) dx -=global.SX;
if (dx <= -global.halfSX) dx +=global.SX;
if (dy > global.halfSY) dx -=global.SY;
if (dy <= -global.halfSY) dx +=global.SY;

r = sqrt( dx*dx + dy*dy );

return r;
}


void init_particles_randomly()
{
int i,j;
double x_try,y_try;
double r_min;
double dr;
int overlap;
int N_trials;

r_min = 0.25;

for(i=0;i<global.N_particles;i++)
    {
    x_try = 0.0;
    y_try = 0.0;
    
    //check overlap with previous particles
    //assume there is overlap to get into the cycle
    overlap = 1;
    N_trials = 0;
    
    while ((overlap==1)&&(N_trials<global.N_particles))
        {
   
        //attempt to place the particle
        x_try = global.SX * rand()/(RAND_MAX+1.0);
        y_try = global.SY * rand()/(RAND_MAX+1.0);
 
        //assume this was good
        overlap = 0;
        
        for(j=0;j<i;j++)
            {
            //calculate distance
            dr = distance_folded_PBC(x_try, y_try, global.particle_x[j], global.particle_y[j]);
            if (dr < r_min)
                {
                overlap = 1; //found overlap
                N_trials++;
                break; //no need to check with other particles
                }
            }
        }
        
    if (N_trials==global.N_particles * global.N_particles)
        {
        printf("Can't place particles randomly, (system too dense) quitting\n");
        exit(1);
        }
      
    global.particle_x[i] = x_try;
    global.particle_y[i] = y_try;
    global.particle_fx[i] = 0.0;
    global.particle_fy[i] = 0.0;
    
    global.particle_color[i] = 3;
    global.particle_direction[i] = -1.0;
  
    // if (rand()/(RAND_MAX+1.0)<0.5)
    //     {
    //     global.particle_direction[i] = - 1.0;
    //     global.particle_color[i] = 2;
    //     }
    // else
    //     {
    //     global.particle_direction[i] = + 1.0;
    //     global.particle_color[i] = 3;
    //     }
    
    }

printf("Random arrangement of particles initialized\n");
printf("N_particles = %d placed\n",global.N_particles);
}

void init_files()
{
global.moviefile = fopen("particles.mvi","wb");
if (global.moviefile == NULL)
    {
    printf("Could not create/open movie file\n");
    exit(2);
    }
}


//PINNINGSITES
void init_pinningsites_square_lattice()
{
int i,j,k;
int N_rows,N_columns;
double lattice_const;

N_rows = (int)sqrt((double)global.N_pinningsites);
N_columns = N_rows;

lattice_const = global.SX/(double)N_rows;
global.pinning_lattice_constant = lattice_const;

k=0;

for(i=0;i<N_rows;i++)
    for(j=0;j<N_columns;j++)
        {
        global.pinningsite_x[k] = (i + 0.5) * lattice_const;
        global.pinningsite_y[k] = (j + 0.5) * lattice_const;
        global.pinningsite_dx_so_far[k] = 0.0;
        global.pinningsite_dy_so_far[k] = 0.0;
        global.pinningsite_fx[k] = 0.0;
        global.pinningsite_fy[k] = 0.0;
        global.pinningsite_R[k] = global.pinningsite_setradius;
        global.pinningsite_fmax[k] = 0.5;
        global.pinningsite_color[k] = 2.0;
        k++;
        }

printf("pinningsites initialized\n");
printf("N_pinningsites = %d\n",k);
printf("Square lattice\nRows = %d\nColumns = %d\nLattice constant = %.2lf\n",N_rows,N_columns, lattice_const);
}

void init_pinningsites()
{
// global.N_pinningsites = 576;
global.N_pinningsites = 200;
global.pinningsite_x = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_y = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_fx = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_fy = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_color = (int *)malloc(global.N_pinningsites*sizeof(int));
global.pinningsite_direction_x = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_direction_y = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_dx_so_far = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_dy_so_far = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_R = (double *)malloc(global.N_pinningsites*sizeof(double));
global.pinningsite_fmax = (double *)malloc(global.N_pinningsites*sizeof(double));

init_pinningsites_square_lattice();
global.pinning_direction_change = 1;
}



