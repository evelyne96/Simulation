//
//  main.c
//  softsim
//
//  Created by András Libál on 7/18/18.
//  Copyright © 2018 András Libál. All rights reserved.
//
// first problem: tim,e crystal simulation

#include <stdio.h>
#include "initializer.h"
#include "running.h"

int main(int argc, const char * argv[]) {
    
    printf("Soft Matter Simulator\n");
    

    // read_init_file("parameters/parameters.txt");
   
    init_simulation();
    init_simulation_box();
    // init_pinningsites();
    init_particles();
    init_files();
    
    // run_simulation();
    
    run_test_verlet_vs_time();

    return 0;
}
