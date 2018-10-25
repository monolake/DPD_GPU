#ifndef SIM_RUN
#define SIM_RUN

const int timestep=100;
const double dt=0.005;
const int force_a0=25;
const double force_gamma=4.5;
const double force_sigma=3.0;
const double rc=1.0;

#include <chrono>
#include <random>
#include "calc_neighbor_list.h"
#include <iostream>
using namespace std;

void clear_force(double** f, int N);
void compute_force(double** r, double** v, double** f, double rand_num, int N, OUTPUT_struct* cell_list);
#endif 
