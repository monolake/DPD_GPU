#include "sim_run.h"
#include "random_mars.h"

void clear_force(double **f, int N) {
    for (int i = 0; i < N; i++) {
        f[i][0] = 0;
        f[i][1] = 0;
        f[i][2] = 0;    
    }
}

void cell_force(double** r, double** v, double** f, RanMars* random, int N, OUTPUT_struct* cell_list, 
int cell_id, int i) {

    int* cell_partc_list_res = cell_list->cell_partc_list_res[cell_id];
    int part_num = cell_list->cell_partc_num_res[cell_id];
    
    for (int j = 0; j < part_num; j++) {
        int part_id = cell_partc_list_res[j];
	if (part_id > i) {
	    double delx, dely, delz;
	    delx = r[i][0] - r[part_id][0];
	    dely = r[i][1] - r[part_id][1];
	    delz = r[i][2] - r[part_id][2];
	    if (delx < - 5)
		delx = delx + 10;
	    else if (delx > 5)
		delx = delx - 10;
	    if (dely < -5)
		dely = dely + 10;
	    else if (dely > 5)
		dely = dely - 10;
	    if (delz < -5)
		delz = delz + 10;
	    else if (delz > 5)
		delz = delz - 10;
	    double rr;
	    rr = sqrt(delx * delx + dely * dely + delz * delz);
	    if(rr < rc) {
		double fpair;
		double wr;
		wr = 1 - rr / rc;
		fpair = force_a0 * wr;

		double delvx, delvy, delvz;
		delvx = v[i][0] - v[part_id][0];
		delvy = v[i][1] - v[part_id][1];
		delvz = v[i][2] - v[part_id][2];

		double dot;
		dot = (delx * delvx + dely * delvy + delz * delvz) / rr;
		fpair -= force_gamma * wr * wr * dot;

		fpair += force_sigma * wr * random->gaussian() * 1 / sqrt(dt);

		f[i][0] += delx * fpair / rr;
		f[i][1] += dely * fpair / rr;
		f[i][2] += delz * fpair / rr;
		
		f[part_id][0] -= delx * fpair / rr;
		f[part_id][1] -= dely * fpair / rr;
		f[part_id][2] -= delz * fpair / rr;
	    }
	}
   }
}

void compute_force(double** r, double** v, double** f, RanMars * random, int N, OUTPUT_struct* cell_list , double len_cell,
double len_x, double len_y, double len_z) {

    double m = 1.0;
    int num_cx, num_cy, num_cz;
    num_cx = (int) floor(len_x/len_cell);// + 1;
    num_cy = (int) floor(len_y/len_cell);// + 1;
    num_cz = (int) floor(len_z/len_cell);// + 1;
    for(int i = 0; i < N - 1; i++) {
        double x = r[i][0];
        double y = r[i][1];
        double z = r[i][2];
        int idx = (int) floor(x / len_cell);
        int idy = (int) floor(y / len_cell);
        int idz = (int) floor(z / len_cell);
        //int cell_id = idz + idy * num_cz + idx * num_cy * num_cz;
        // loop through 27 boxes
	int l, m, n;
	for (int n = -1; n < 2; n++) { 
		for (int m = -1; m < 2; m++) {
			for (int l = -1; l < 2; l++) {
			    int newidz = (idz + n + num_cz) % num_cz;
			    int newidy = (idy + m + num_cy) % num_cy;
			    int newidx = (idx + l + num_cx) % num_cx;
			    int cell_id = newidz + newidy * num_cz + newidx * num_cy * num_cz;
   			    cell_force(r, v, f, random, N, cell_list, cell_id, i); 
			}
		}
	}
   }
}
void compute_force_vector(double** r, double** v, double** f, RanMars * random, int N, map<int,vector<int>> cell_list, double len_cell,
double len_x, double len_y, double len_z) {

    double m = 1.0;
    int num_cx, num_cy, num_cz;
    num_cx = (int) floor(len_x/len_cell) + 1;
    num_cy = (int) floor(len_y/len_cell) + 1;
    num_cz = (int) floor(len_z/len_cell) + 1;
    for(int i = 0; i < N - 1; i++) {
        double x = r[i][0];
        double y = r[i][1];
        double z = r[i][2];
        int idx = (int) floor(x / len_cell);
        int idy = (int) floor(y / len_cell);
        int idz = (int) floor(z / len_cell);
        // loop through 27 boxes
        for (int l = -1; l < 2; l++) {
            for (int m = -1; m < 2; m++) {
                for (int n = -1; n < 2; n++) {
                    
                    int newidz = (idz + n + num_cz) % num_cz;
                    int newidy = (idy + m + num_cy) % num_cy;
                    int newidx = (idx + l + num_cx) % num_cx;
                    
                    int cell_id = newidz + newidy * num_cz + newidx * num_cy * num_cz;
                    // cout << "particle i " << i << " " << cell_list[cell_id].size() << endl;
                    for (int j = 0; j < cell_list[cell_id].size(); j++) {
                        int part_id = cell_list[cell_id][j];
                        if (part_id > i) {
                            double delx, dely, delz;
                            delx = r[i][0] - r[part_id][0];
                            dely = r[i][1] - r[part_id][1];
                            delz = r[i][2] - r[part_id][2];
                            if (delx < - 5)
                                delx = delx + 10;
                            else if (delx > 5)
                                delx = delx - 10;
                            if (dely < -5)
                                dely = dely + 10;
                            else if (dely > 5)
                                dely = dely - 10;
                            if (delz < -5)
                                delz = delz + 10;
                            else if (delz > 5)
                                delz = delz - 10;
                            double rr;
                            rr = sqrt(delx * delx + dely * dely + delz * delz);
                            if(rr < rc) {
                                double fpair;
                                double wr;
                                wr = 1 - rr / rc;
                                fpair = force_a0 * wr;

                                double delvx, delvy, delvz;
                                delvx = v[i][0] - v[part_id][0];
                                delvy = v[i][1] - v[part_id][1];
                                delvz = v[i][2] - v[part_id][2];

                                double dot;
                                dot = (delx * delvx + dely * delvy + delz * delvz) / rr;
                                fpair -= force_gamma * wr * wr * dot;

                                fpair += force_sigma * wr * random->gaussian() * 1 / sqrt(dt);

                                f[i][0] += delx * fpair / rr;
                                f[i][1] += dely * fpair / rr;
                                f[i][2] += delz * fpair / rr;
                                
                                f[part_id][0] -= delx * fpair / rr;
                                f[part_id][1] -= dely * fpair / rr;
                                f[part_id][2] -= delz * fpair / rr;
       
                            }
                        } 
                    }
                }
            }
        }  
   }
}
void compute_force_std(double** r, double** v, double** f, RanMars * random, int N) {
    double m = 1.0;
    for(int i = 0; i < N - 1; i++) {
            for (int j = i+1; j < N; j++) {
                if (i != j) {
                    double delx, dely, delz;
                    delx = r[i][0] - r[j][0];
                    if (delx < - 5)
                        delx = delx + 10;
                    else if (delx > 5)
                        delx = delx - 10;
                    dely = r[i][1] - r[j][1];
                    if (dely < -5)
                        dely = dely + 10;
                    else if (dely > 5)
                        dely = dely - 10;
                    delz = r[i][2] - r[j][2];
                    if (delz < -5)
                        delz = delz + 10;
                    else if (delz > 5)
                        delz = delz - 10;
                    double rr;
                    rr = sqrt(delx * delx + dely * dely + delz * delz);
                    if(rr < rc) {
                        double fpair;
                        double wr;
                        wr = 1 - rr / rc;
                        fpair = force_a0 * wr;

                        double delvx, delvy, delvz;
                        delvx = v[i][0] - v[j][0];
                        delvy = v[i][1] - v[j][1];
                        delvz = v[i][2] - v[j][2];

                        double dot;
                        dot = (delx * delvx + dely * delvy + delz * delvz) / rr;
                        fpair -= force_gamma * wr * wr * dot;

                        fpair += force_sigma * wr * random->gaussian() * 1 / sqrt(dt);

                        f[i][0] += delx * fpair / rr;
                        f[i][1] += dely * fpair / rr;
                        f[i][2] += delz * fpair / rr;
                        
                        f[j][0] -= delx * fpair / rr;
                        f[j][1] -= dely * fpair / rr;
                        f[j][2] -= delz * fpair / rr;
                    }
                }
            }
   }
}
void writeDump(ofstream& outputfile, double** r, double** v, int ntimestep) {

    outputfile << "ITEM: TIMESTEP" << endl;
    outputfile << ntimestep << endl;
    outputfile << "ITEM: NUMBER OF ATOMS" << endl;
    outputfile << N << endl;
    outputfile << "ITEM: BOX BOUNDS pp pp pp" << endl;
    outputfile <<"-100 100" << endl;
    outputfile <<"-100 100" << endl;
    outputfile <<"-100 100" << endl;
    outputfile << "ITEM: ATOMS id type x y z vx vy vz" << endl;
    for(int i = 0; i < N; i++) {
        outputfile << i+1 << " " << 1 << " " << r[i][0] << " " << r[i][1] << " " << r[i][2] << " " << \
        v[i][0] << " " << v[i][1] << " " << v[i][2] << endl;
    }
}

void pbc(double** r) {
    for (int i = 0; i < N; i++) {
        if (r[i][0] < 0)
            r[i][0] += 10;
        if (r[i][0] > 10)
            r[i][0] -= 10;
        if (r[i][1] < 0)
            r[i][1] += 10;
        if (r[i][1] > 10)
            r[i][1] -= 10;
        if (r[i][2] < 0)
            r[i][2] += 10;
        if (r[i][2] > 10)
            r[i][2] -= 10;
    }
}

double computeKE(double** v) {
    double m = 1.0;
    double ke = 0.0;
    for (int i = 0; i < N; i++) {
        ke += 0.5 * m * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
    }
    return ke;
}


