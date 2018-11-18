#ifndef SIM_RUN
#define SIM_RUN

const int timestep=100;
const double dt=0.005;
const int force_a0=18.75;
const double force_gamma=4.5;
const double force_sigma=3.0;
const double rc=1.0;
const int N = 4000;
#include <chrono>
#include <random>
#include "calc_neighbor_list.h"
#include <iostream>
#include <fstream>
#include "random_mars.h"
#include <map>
#include <vector>
using namespace std;

void clear_force(double** f, int N);
void cell_force(double** r, double** v, double** f, RanMars* random, int N, OUTPUT_struct* cell_list, int cell_id, int i);
void compute_force(double** r, double** v, double** f, RanMars * random, int N, OUTPUT_struct* cell_list, double len_cell,
double box_lenx, double box_leny, double box_lenz );
void compute_force_vector(double** r, double** v, double** f, RanMars * random, int N, map<int, vector<int>> cell_list, double len_cell,
double box_lenx, double box_leny, double box_lenz); 
void compute_force_std(double** r, double** v, double** f, RanMars * random, int N);
void writeDump(ofstream& outputfile, double** r, double** v, int ntimestep);
void pbc(double** r);
double computeKE(double** v);
#endif 
