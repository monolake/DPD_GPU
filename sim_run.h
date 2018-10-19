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


void run();
void compute_force(double** r, double** v, double** f, double rand_num, int N);
#endif 
