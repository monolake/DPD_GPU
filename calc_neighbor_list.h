#ifndef CALC_NEIGHBOR_LIST
#define CALC_NEIGHBOR_LIST
#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <math.h>
#include    <iostream>
#include    <fstream>

const int       CONST_MAX_RAND = 1000;


typedef struct OUTPUT_struct{
    int             cell_num;
    int             col_num; // return max #partices in a cell
    int*            cell_partc_num_res; // return num of particles in each cell
    int**           cell_partc_list_res; // return cell list
    
    int**           cell_id2loc;
}OUTPUT_struct;

OUTPUT_struct* func_partc_incell_stat(double** partc_pos, int partc_num, double len_cell, double len_x, double len_y, double len_z);
void func_print_output_info(OUTPUT_struct*  output_info);
#endif
