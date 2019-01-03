//
//  main.c
//  dla_fractal


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_grid 1000
#define PI 3.141592653579892

int grid[N_grid][N_grid];

double R_max; // maximum radius of the crystal
double sticking_prob;

int randomwalker_x, randomwalker_y; // the random walker's position
double randomwalker_R;  //random walker's distance from the center

int t;

FILE *moviefile;

int N_particles;

void write_cmovie(void);

void init_grid()
{
    int i,j;
    
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            grid[i][j] = 0;
    
    grid[N_grid/2][N_grid/2] = 1;
    
    R_max = 5.0;
    N_particles = 1;
}


void init_randomwalker()
{
    double theta;
    
    if (( (double)N_grid/2 + 3.0*R_max ) > (double) N_grid )
        {
            printf("Cannot create any more random walkers");
            exit(1);
        }
    
    randomwalker_R = (R_max + 5.0);
    theta = 2*PI*rand()/(RAND_MAX+1.0);
    
    randomwalker_x = N_grid/2 + (int)floor(randomwalker_R * cos(theta));
    randomwalker_y = N_grid/2 + (int)floor(randomwalker_R * sin(theta));
    
}

void move_randomwalker()
{
    double r;
    int tempx, tempy;
    
    tempx = randomwalker_x;
    tempy = randomwalker_y;
    
    r = rand()/(RAND_MAX+1.0);
    
    if (r<0.25) tempx++;        //right
    else if (r<0.5) tempy++;    //up
    else if (r<0.75) tempx--;   //left
    else tempy--;               //down
    
    if (grid[tempx][tempy]==0)
        {
            randomwalker_x = tempx;
            randomwalker_y = tempy;
        }
}

void does_it_stick()
{
    double R_max_candidate2;
    double dx,dy;
    
    double r;
    
    r = rand()/(RAND_MAX+1.0);
  
    if (r<sticking_prob)
    //if it has a neighbor grid site that is occupied
    //nearest neighbor (right, left,up,down)
    // grid[randomwalker_x+1][randomwalker_y+1]
    //would be next nearest neighbor
    if ((grid[randomwalker_x+1][randomwalker_y]==1)||
        (grid[randomwalker_x-1][randomwalker_y]==1)||
        (grid[randomwalker_x][randomwalker_y+1]==1)||
        (grid[randomwalker_x][randomwalker_y-1]==1))
    {
        
        //printf("sticking at t=%d\n",t);
        
        //particle got stuck on the growing crystal
        grid[randomwalker_x][randomwalker_y] = 1;
        N_particles++;
        write_cmovie();
        
        dx = (randomwalker_x-(N_grid/2));
        dy = (randomwalker_y-(N_grid/2));
        
        R_max_candidate2 = dx*dx + dy*dy;
        
        if (R_max_candidate2 > R_max*R_max)
            R_max = sqrt(R_max_candidate2);
        
        init_randomwalker();
    }

}


int the_random_walker_is_lost()
{
    double dx,dy, dr2;
    
    dx  = randomwalker_x - (N_grid/2.0);
    dy  = randomwalker_y - (N_grid/2.0);

    dr2 = dx*dx + dy*dy;
    
    if (dr2>9.0*R_max*R_max) return 1;
    else return 0;
}

void write_cmovie()
{
    int i,j,k;
    float floatholder;
    int intholder;
    
    intholder = N_particles;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    intholder = t;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    k=0;
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            if (grid[i][j])
            {
            intholder = 3; //color of spin
            //else intholder = 3 + (klaszterszam[i][j]%10);
            fwrite(&intholder,sizeof(int),1,moviefile);
            intholder = k++;            //spin ID
            fwrite(&intholder,sizeof(int),1,moviefile);
            floatholder = (float)i;
            fwrite(&floatholder,sizeof(float),1, moviefile);
            floatholder = (float)j;
            fwrite(&floatholder,sizeof(float),1,moviefile);
            floatholder = 1.0;//cum_disp, cmovie format
            fwrite(&floatholder,sizeof(float),1,moviefile);
        }
}

double avg_distance() {
    double sum_dist = 0.0, N2 = N_grid/2;
    for(int i=0;i<N_grid;i++)
        for(int j=0;j<N_grid;j++) 
            if(grid[i][j]==1)
                sum_dist += sqrt( (N2-i)*(N2-i) + (N2-j)*(N2-j));

    return sum_dist/N_particles;
}

int main(int argc, const char * argv[]) {
  
    printf("Diffusion Limited Aggregation\n");

    moviefile = fopen("fraktal.mvi","wb");
    // sticking_prob = 0.5;
    sticking_prob = 0.1;

    FILE *testf;
    
    while(sticking_prob<=1.0) {
    sticking_prob += 0.1;    

        init_grid();
        init_randomwalker();
        
        char buffer[1024];
        snprintf(buffer, sizeof(buffer), "test/res%f.csv", sticking_prob);
        testf = fopen(buffer,"wt");
        
        for(t=0;t<10000;t++)
            {
                if (t%10==0) printf("t = %d\n",t);
                move_randomwalker();
                if (the_random_walker_is_lost()) init_randomwalker();
                does_it_stick();
                if (t%100==0) fprintf(testf,"%d,%f\n",t,avg_distance());
            }
        // write_cmovie();
    }
    fclose(moviefile);
    fclose(testf);
    return 0;
}
