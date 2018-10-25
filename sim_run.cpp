#include "sim_run.h"

void clear_force(double **f, int N) {
    for (int i = 0; i < N; i++) {
        f[i][0] = 0;
        f[i][1] = 0;
        f[i][2] = 0;    
    }
}

void compute_force(double** r, double** v, double** f, double rand_num, int N, OUTPUT_struct* cell_list) {

    double m = 1.0;
    int num_cx, num_cy, num_cz;
    double len_x = 2.0, len_y = 2.0, len_z = 5.0;
    double len_cell = 1.0;
    num_cx = (int) ceil(len_x/len_cell);
    num_cy = (int) ceil(len_y/len_cell);
    num_cz = (int) ceil(len_z/len_cell);
    for(int i = 0; i < N; i++) {
        double x = r[i][0];
        double y = r[i][1];
        double z = r[i][2];
        int idx = (int) floor(x / len_cell);
        int idy = (int) floor(y / len_cell);
        int idz = (int) floor(z / len_cell);
        int cell_id = idz + idy * num_cz + idx * num_cy * num_cz;
        // loop through 27 boxes
        for (int l = -1; l < 2; l++) {
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
                    int cell_id = (idz+n) + (idy+m) * num_cz + (idx+l) * num_cy * num_cz;
                    int* cell_partc_list_res = cell_list->cell_partc_list_res[cell_id];
                    int part_num = cell_list->cell_partc_num_res[cell_id];
                    for (int j = 0; j < part_num; j++) {
                        int part_id = cell_partc_list_res[j];
                        if (i != part_id) {
                            double delx, dely, delz;
                            delx = r[i][0] - r[part_id][0];
                            dely = r[i][1] - r[part_id][1];
                            delz = r[i][2] - r[part_id][2];
                            double r;
                            r = sqrt(delx * delx + dely * dely + delz * delz);
                            if(r < rc) {
                                double fpair;
                                double wr;
                                wr = 1 - r / rc;
                                fpair = force_a0 * wr;

                                double delvx, delvy, delvz;
                                delvx = v[i][0] - v[part_id][0];
                                delvy = v[i][1] - v[part_id][1];
                                delvz = v[i][2] - v[part_id][2];

                                double dot;
                                dot = (delx * delvx + dely * delvy + delz * delvz) / r;
                                fpair -= force_gamma * wr * wr * dot;

                                //randum - gaussian random number with zero mean and unit variance
                                fpair += force_sigma * wr * rand_num * 1 / sqrt(dt);

                                f[i][0] += delx / r * fpair;
                                f[i][1] += dely / r * fpair;
                                f[i][2] += delz / r * fpair;
                            }
                        } 
                    }
                }
            }
        }  
   }
}

