#include	<stdio.h>
#include 	<stdlib.h>
#include 	<string.h>
#include	<math.h>
#include	<iostream>
#include	<fstream>
#include    "sim_run.h"
#include    "calc_neighbor_list.h"
#include    "random_mars.h"
#include    "cuda.h"
#include <map>
#include <vector>
using namespace std;

void cpu_dpd(double** r, double**f, double** v, const char* select, RanMars* random, OUTPUT_struct* output_info, double len_cell, double len_x, double len_y, double len_z) {
        ofstream outputfile;
        outputfile.open("dump.md", ios::out);

        //init force compute
        map<int, vector<int>> cell_list;
        if (!strcmp(select, "sijun"))
            compute_force(r, v, f, random, N, output_info, len_cell, len_x, len_y, len_z);
        else if (!strcmp(select, "vector")) {
            int cell_nx = 10 / len_cell + 1;
            int cell_ny = 10 / len_cell + 1;
            int cell_nz = 10 / len_cell + 1;
                    for (int i = 0; i < N; i++) {
                int cellx = (int) r[i][0] / len_cell;
                int celly = (int) r[i][1] / len_cell; 
                int cellz = (int) r[i][2] / len_cell;
                int cellid = cellz + celly * cell_nz + cellx * cell_nz * cell_ny;
                cell_list[cellid].push_back(i); 
            }
            compute_force_vector(r, v, f, random, N, cell_list, len_cell, len_x, len_y, len_z);
        }
        else if (!strcmp(select, "base")) 
            compute_force_std(r, v, f, random, N);

        //writeDump(outputfile, r, v, 0);

        int ntimestep = 5000;
        double m = 1.0;
        for (int i = 0; i <= ntimestep; i++) {

            //half integration
            for(int j = 0; j < N; j++) {
                v[j][0] += 0.5 * f[j][0] * dt;
                v[j][1] += 0.5 * f[j][1] * dt;
                v[j][2] += 0.5 * f[j][2] * dt;

                r[j][0] += v[j][0] * dt;
                r[j][1] += v[j][1] * dt;
                r[j][2] += v[j][2] * dt;
            }
            pbc(r); 
            //force computation
            clear_force(f, N);
            if (!strcmp(select, "sijun"))
                compute_force(r, v, f, random, N, output_info, len_cell, len_x, len_y, len_z);
            else if (!strcmp(select, "vector"))
                compute_force_vector(r, v, f, random, N, cell_list, len_cell, len_x, len_y, len_z);
            else if (!strcmp(select, "base"))
                compute_force_std(r, v, f, random, N);

            //full integration
            for(int j = 0; j < N; j++) {
                v[j][0] += 0.5 * f[j][0] * dt;
                v[j][1] += 0.5 * f[j][1] * dt;
                v[j][2] += 0.5 * f[j][2] * dt;
            }
             if(i % 1 == 0) {
                double ke = computeKE(v);
                cout << i << " temp is " << ke * 2 / (3 * 4000 * 1) << endl;
                //writeDump(outputfile, r, v, i);
            }
        }
}
__device__ void compute_force_gpu(int id, double* r, double* f, double* v, RanMars * random, int N, int* cell_list, int* cell_list_count,
double len_cell, double len_x, double len_y, double len_z, int avg_num_cell) {
   double m = 1.0;
    int num_cx, num_cy, num_cz;
    num_cx = (int) floor(len_x/len_cell) + 1;
    num_cy = (int) floor(len_y/len_cell) + 1;
    num_cz = (int) floor(len_z/len_cell) + 1;
    double x = r[id * 3];
    double y = r[id * 3 + 1];
    double z = r[id * 3 + 2];
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
                int num_neigh_particle = cell_list_count[cell_id];
                for (int j = 0; j < num_neigh_particle; j++) {
                    int part_id = cell_list[cell_id * avg_num_cell + j];
                    if (part_id > id) {
                        double delx, dely, delz;
                        delx = x - r[part_id * 3];
                        dely = y - r[part_id * 3 + 1];
                        delz = z - r[part_id * 3 + 2];
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
                            delvx = v[id*3] - v[part_id*3];
                            delvy = v[id*3+1] - v[part_id*3+1];
                            delvz = v[id*3+2] - v[part_id*3+2];

                            double dot;
                            dot = (delx * delvx + dely * delvy + delz * delvz) / rr;
                            fpair -= force_gamma * wr * wr * dot;

                            fpair += force_sigma * wr * random->gaussian() * 1 / sqrt(dt);

                            f[id*3] += delx * fpair / rr;
                            f[id*3+1] += dely * fpair / rr;
                            f[id*3+2] += delz * fpair / rr;

                            f[part_id*3] -= delx * fpair / rr;
                            f[part_id*3+1] -= dely * fpair / rr;
                            f[part_id*3+2] -= delz * fpair / rr;

                        }
                    }
                }
            }
        }
    }
}

__global__ void iteration(double* r, double* f, double* v, int* cell_list, int* cell_list_count, RanMars* random, int avg_num_cell) {
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    if (id < N) {
            int i;
            // init integration
            for(i = 0; i < 3; i++) {
                v[id * 3 + i] += 0.5 * f[id * 3 + i] * dt;
                r[id * 3 + i] += v[id * 3 + i] * dt;
            } 
            // do periodic boundary condition
            for (i = 0; i < 3; i++) {
                if (r[id * 3 + i] < 0)
                    r[id * 3 + i] += 10;
                else if (r[id * 3 + i] > 10)
                    r[id * 3 + i] -= 10;
            }
            //force computation
            for (i = 0; i < 3; i++)
                f[id * 3 + i] = 0;
            double len_cell = 2.0;
            double len_x = 10.0;
            double len_y = 10.0;
            double len_z = 10.0;
            compute_force_gpu(id, r, f, v, random, N, cell_list, cell_list_count, len_cell, len_x, len_y, len_z, avg_num_cell);

            //full integration
            for (i = 0; i < 3; i++) {
                v[id * 3 + i] += 0.5 * f[id * 3 + i] * dt;
            }
 
    }
}

void gpu_dpd(double** r, double**f, double** v, const char* select, RanMars* random, OUTPUT_struct* output_info, double len_cell, double len_x, double len_y, double len_z) {

        //allocate r,f,v on GPU
        int size = N * 3 * sizeof(double);
        double* cpu_r = (double *) malloc(size);
        double* cpu_f = (double *) malloc(size);
        double* cpu_v = (double *) malloc(size);
        int i = 0;
        int count = 0;
        for(; i < N; i++) {
            cpu_r[count] = r[i][0];
            cpu_f[count] = f[i][0];
            cpu_v[count] = v[i][0];
            cpu_r[count+1] = r[i][1];
            cpu_f[count+1] = f[i][1];
            cpu_v[count+1] = v[i][1];
            cpu_r[count+2] = r[i][2];
            cpu_f[count+2] = f[i][2];
            cpu_v[count+2] = v[i][2];
            count += 3; 
        }
        double* dev_r;
        double* dev_f;
        double* dev_v;
        size_t pitch;
        cudaMalloc((void **) dev_r, size);
        cudaMalloc((void **) dev_f, size);
        cudaMalloc((void **) dev_v, size);
        cudaMemcpy(dev_r, r, size, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_f, f, size, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_v, v, size, cudaMemcpyHostToDevice);

        // build cell list on cpu 
        int ncell_x = 10 / len_cell;
        int ncell_y = 10 / len_cell;
        int ncell_z = 10 / len_cell;
        int num_cell = ncell_x * ncell_y * ncell_z;
        int avg_num_particle = N / num_cell + 1;
        int* cell_list = (int *) malloc(sizeof(int) * num_cell * avg_num_particle);
        int* cell_list_count = (int *) malloc(sizeof(int) * num_cell);
        for (i = 0; i < num_cell * avg_num_particle; i++) {
            cell_list[i] = 0;
        }
        for (i = 0; i < num_cell; i++) {
            cell_list_count[i] = 0;
        }

        i = 0;
        for(; i < N; i++) {
            int cellx = (int) r[i][0] / len_cell;
            int celly = (int) r[i][1] / len_cell;
            int cellz = (int) r[i][2] / len_cell;
            int cellid = cellz + celly * ncell_z + cellx * ncell_z * ncell_y;
            int row = cellid;
            int col = cell_list_count[cellid];
            int index = row * avg_num_particle + col;
            cell_list[index] = i;
            cell_list_count[i]++; 
        }

        int* dev_cell_list;
        int* dev_cell_list_count;
        int size_cell_list = sizeof(int) * num_cell * avg_num_particle;
        cudaMalloc((void **) &dev_cell_list, size_cell_list);
        cudaMalloc((void **) &dev_cell_list_count, sizeof(int) * num_cell);
        cudaMemcpy(dev_cell_list, cell_list, size_cell_list, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_cell_list_count, cell_list_count, sizeof(int) * num_cell, cudaMemcpyHostToDevice);
 
        ofstream outputfile;
        outputfile.open("dump.md", ios::out);

        //init force compute
        //map<int, vector<int>> cell_list;
        //if (!strcmp(select, "sijun"))
        //    compute_force(r, v, f, random, N, output_info, len_cell, len_x, len_y, len_z);
        //else if (!strcmp(select, "vector")) {
        //    int cell_nx = 10 / len_cell + 1;
        //    int cell_ny = 10 / len_cell + 1;
        //    int cell_nz = 10 / len_cell + 1;
        //            for (int i = 0; i < N; i++) {
        //        int cellx = (int) r[i][0] / len_cell;
        //        int celly = (int) r[i][1] / len_cell; 
        //        int cellz = (int) r[i][2] / len_cell;
        //        int cellid = cellz + celly * cell_nz + cellx * cell_nz * cell_ny;
        //        cell_list[cellid].push_back(i); 
        //    }
        //    compute_force_vector(r, v, f, random, N, cell_list, len_cell, len_x, len_y, len_z);
        //}
        //else if (!strcmp(select, "base")) 
        //    compute_force_std(r, v, f, random, N);

        //writeDump(outputfile, r, v, 0);

        int ntimestep = 5000;
        double m = 1.0;
        for (int i = 0; i <= ntimestep; i++) {
            iteration<<<N/1024+1, N>>>(dev_r, dev_f, dev_v, dev_cell_list, dev_cell_list_count, random, avg_num_particle);
            cudaMemcpy(r, dev_r, size, cudaMemcpyDeviceToHost);
            //build new cell list
         int ncell_x = 10 / len_cell;
        int ncell_y = 10 / len_cell;
        int ncell_z = 10 / len_cell;
        int num_cell = ncell_x * ncell_y * ncell_z;
        int avg_num_particle = N / num_cell + 1;
        int* cell_list = (int *) malloc(sizeof(int) * num_cell * avg_num_particle);
        int* cell_list_count = (int *) malloc(sizeof(int) * num_cell);
        for (i = 0; i < num_cell * avg_num_particle; i++) {
            cell_list[i] = 0;
        }
        for (i = 0; i < num_cell; i++) {
            cell_list_count[i] = 0;
        }

        i = 0;
        for(; i < N; i++) {
            int cellx = (int) r[i][0] / len_cell;
            int celly = (int) r[i][1] / len_cell;
            int cellz = (int) r[i][2] / len_cell;
            int cellid = cellz + celly * ncell_z + cellx * ncell_z * ncell_y;
            int row = cellid;
            int col = cell_list_count[cellid];
            int index = row * avg_num_particle + col;
            cell_list[index] = i;
            cell_list_count[i]++; 
        }

           
            if(i % 1 == 0) {
                double ke = computeKE(v);
                cout << i << " temp is " << ke * 2 / (3 * 4000 * 1) << endl;
                //writeDump(outputfile, r, v, i);
            }
        }
        free(cpu_r);
        free(cpu_f);
        free(cpu_v);
        free(cell_list);
        free(cell_list_count);
        cudaFree(dev_r);
        cudaFree(dev_f);
        cudaFree(dev_v);
        cudaFree(dev_cell_list);
        cudaFree(dev_cell_list_count);
}

void next_func(FILE* fptr) {
    
	char		str_buff[256];
	
	//while (!feof(fptr)){
	//	fscanf(fptr, "%s", str_buff);
	//}
	
	for (int loop=0; loop<14; loop++){
		
		fgets(str_buff, 256, fptr);
		
		printf("%d %s", loop, str_buff);
		
	}
	
}


int load_func(FILE* fptr, double* outptr) {
    
	char		str_buff[256];
	double		b,c,d;
	int			a;
	
	fscanf(fptr, "%d", &a);
	fscanf(fptr, "%d", &a);
	
	if (!feof(fptr)){
		fscanf(fptr, "%lf", &b);
		fscanf(fptr, "%lf", &c);
		fscanf(fptr, "%lf", &d);
		outptr[0] = b;
		outptr[1] = c;
		outptr[2] = d;
		return 1;
	}
	else{
		return 0;
	}
	
}

int main(int argc, char* argv[])
{

    int type_of_device = 0; // 0 - CPU; 1 - GPU
    const char* select = argv[1];
	FILE*		file_ptr;
	
	char		str_input[5];
	double		result[3];
	int			count;
	int			end;
	int			flag;
    double** partc_pos_res;
    int pos_index;
    int partc_num_def = 4000;
    int i;
	
	// read input position file
	count = 0;
	if((file_ptr = fopen("4000_new.txt","r")) == NULL){
		printf("Cannt open the file!");
		exit(1);
	}
	
	next_func(file_ptr);
	

    // allocate memory for partc_pos_res
    partc_pos_res = (double**)malloc(partc_num_def*sizeof(double*));
    for (i=0; i<partc_num_def; i++){
        partc_pos_res[i] = (double*)malloc(3*sizeof(double));
        
    }

    // fill in the particle positions into partc_pos_res
	while (!feof(file_ptr)){
		// count = count+1;
		flag = load_func(file_ptr, result);
		if (flag){

			// printf("Line:%d, %.2f, %.2f, %.2f\n", count, result[0], result[1], result[2]);
            for (pos_index=0; pos_index<3; pos_index++){
                partc_pos_res[count][pos_index] = result[pos_index];
            }
			
		}
        count = count+1;
		
	}


	fclose(file_ptr);
	//scanf("%d",&end);

    // build cell list
    int     loop;
    unsigned int seed;
    
    //input_cube
    double  len_cell;
    double  len_x, len_y, len_z;
    
    //input_partc
    int     partc_num;
    double**    partc_pos;
    
    int     end_flag;
    
    //output
    OUTPUT_struct*  output_info;
    
    seed = 10;
    srand(seed);
    
    //INPUT information:
    len_cell = 2.0;
    len_x = 10.0;
    len_y = 10.0;
    len_z = 10.0;
    
    // partc_num = 10;
    // partc_pos = (double**)malloc(partc_num*sizeof(double*));
    // for(loop=0; loop<partc_num; loop++){
    //     partc_pos[loop] = (double*)malloc(3*sizeof(double));
    //     partc_pos[loop][0] = (double)(rand()%CONST_MAX_RAND)/(CONST_MAX_RAND-1)*len_x;
    //     partc_pos[loop][1] = (double)(rand()%CONST_MAX_RAND)/(CONST_MAX_RAND-1)*len_y;
    //     partc_pos[loop][2] = (double)(rand()%CONST_MAX_RAND)/(CONST_MAX_RAND-1)*len_z;
        
    // }
    
    
    
    output_info =  func_partc_incell_stat(partc_pos_res, partc_num_def, len_cell, len_x, len_y, len_z);
    
    
    //PRINT basic information
    //if (1){
        
        
    //}
    
    //PRINT output information
    //if (1){
    //    func_print_output_info(output_info);
    //}
    
    //scanf("%d",&end_flag);
    // finished building cell list

    RanMars * random = new RanMars(34387);
     
    int N = 4000;
    double** r = new double* [N];
    double** v = new double* [N];
    double** f = new double* [N];


    for(int i = 0; i < N; i++) {
        r[i] = new double[3];
        v[i] = new double[3];
        f[i] = new double[3];
    }

    for (int i = 0; i < N; i++) {
        r[i][0] = partc_pos_res[i][0]; 
        r[i][1] = partc_pos_res[i][1];
        r[i][2] = partc_pos_res[i][2];
    }

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < 3; j++) {
            v[i][j] = 0;
            f[i][j] = 0;
        }
    }

    if (!type_of_device) {
        cpu_dpd(r, f, v, select, random, output_info, len_cell, len_x, len_y, len_z);   
    } 
    else {
        gpu_dpd(r, f, v, select, random, output_info, len_cell, len_x, len_y, len_z);   
    }

    //memory release
    for(int i = 0; i < N; i++) {
        delete(r[i]);
        delete(v[i]);
        delete(f[i]);
    }

    delete(r);
    delete(v);
    delete(f);
    delete(random);
    cout << "position is read, f is computed" << endl;
 
}





