#include	<stdio.h>
#include 	<stdlib.h>
#include 	<string.h>
#include	<math.h>
#include	<iostream>
#include	<fstream>
#include    "sim_run.h"
#include    "calc_neighbor_list.h"
#include    "random_mars.h"
//#include    "cuda.h"
#include    "curand_kernel.h"
#include <map>
#include <vector>
using namespace std;

double verify_f(double** f) {
    double sum_f = 0;
    for (int j = 0; j < N; j++) {
        sum_f += f[j][0];
        sum_f += f[j][1];
        sum_f += f[j][2];
    }
    return sum_f;
}

#if __CUDA_ARCH__ < 600
__device__ double atomicAdd_30(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

void cpu_dpd(double** r, double**f, double** v, const char* select, RanMars* random, OUTPUT_struct* output_info, double len_cell, double len_x, double len_y, double len_z) {
        ofstream outputfile;
        outputfile.open("dump_cpu.md", ios::out);
        //init force compute
        if (!strcmp(select, "sijun"))
            compute_force(r, v, f, random, N, output_info, len_cell, len_x, len_y, len_z);
        else if (!strcmp(select, "base")) 
            compute_force_std(r, v, f, random, N);
        //writeDump(outputfile, r, v, 0);
        int ntimestep = 100;//5000;
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
            if (!strcmp(select, "sijun")) {
                compute_force(r, v, f, random, N, output_info, len_cell, len_x, len_y, len_z);
                int partc_num_def = 4000;
                output_info =  func_partc_incell_stat(r, partc_num_def, len_cell, len_x, len_y, len_z);
            }
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
                writeDump(outputfile, r, f, i);
                //cout << "verify force on cpu at step " << i << " is " << verify_f(f) << endl;
            }
        }
}

__global__ void print_kernel() {
    printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
}

__global__ void init(double *r, double *v, double *f) {
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    if (id < N) {
            int i;
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
    }
}

__global__ void f_clear(double *f) {
     int id = threadIdx.x + blockDim.x * blockIdx.x;
    if (id < N) {
           //force computation
            int i;
            for (i = 0; i < 3; i++)
                f[id * 3 + i] = 0;
 }  
}

__global__ void iteration(curandState *state, double* r, double* f, double* v, int* cell_list, int* cell_list_count, int avg_num_cell) {
    //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
    int id = threadIdx.x + blockDim.x * blockIdx.x;
        //printf("blockIdx.x %d threadIdx.x %d\n", blockIdx.x, threadIdx.x); 
    if (id < N) {
                        int i;
            if (id == 0) {
                double sum_f = 0;
                for (int j = 0; j < 3 * N; j++) {
                    sum_f += f[j];
                }
                //printf("sum f is %f \n", sum_f);
            }
            double len_cell = 2.0;
            double len_x = 10.0;
            double len_y = 10.0;
            double len_z = 10.0;
            double m = 1.0;
            int num_cx, num_cy, num_cz;
            num_cx = (int) floor(len_x/len_cell);
            num_cy = (int) floor(len_y/len_cell);
            num_cz = (int) floor(len_z/len_cell);
            //printf("number of cell is %d \n", num_cx * num_cy * num_cz);
            double x = r[id * 3];
            double y = r[id * 3 + 1];
            double z = r[id * 3 + 2];

            // which cell in 3d the particle is in
            int idx = (int) floor(x / len_cell);
            int idy = (int) floor(y / len_cell);
            int idz = (int) floor(z / len_cell);
            curand_init(34387, 0, 0, &state[id]);
            // loop through 27 boxes
            for (int n = -1; n < 2; n++) {
                for (int m = -1; m < 2; m++) {
                    for (int l = -1; l < 2; l++) {
                        int newidz = (idz + n + num_cz) % num_cz;
                        int newidy = (idy + m + num_cy) % num_cy;
                        int newidx = (idx + l + num_cx) % num_cx;
                        int cell_id = newidz + newidy * num_cz + newidx * num_cy * num_cz;
                        if (id == 3999) {
                            //printf("l %d m %d n %d cell_id %d \n", l, m, n, cell_id);
                        }
                        int num_neigh_particle = cell_list_count[cell_id];
                         if (id == 3999) {
                        //    printf("l %d m %d n %d num_neigh_particle %d \n", l, m, n, num_neigh_particle);
                        }
                        for (int j = 0; j < num_neigh_particle; j++) {
                            int part_id = cell_list[cell_id * avg_num_cell + j];
                            if (part_id > id) {
                                if (part_id == id) continue; 
                                double delx, dely, delz;
                                delx = x - r[part_id * 3];
                                dely = y - r[part_id * 3 + 1];
                                delz = z - r[part_id * 3 + 2];

                               if (delx < -5)
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
//printf("blockIdx %d threadIdx %d part_id %d %f %f %f %f %f %f %f\n", blockIdx.x, threadIdx.x, part_id, r[part_id * 3], r[part_id * 3 + 1], r[part_id *3 + 2], x, y, z);
                                double rr;
                                rr = sqrt(delx * delx + dely * dely + delz * delz);
                                                                if(rr < rc) {
                                    //printf("id %f part_id %f \n", id, part_id);
                                if (id == 3999) {
                                    //printf(" id %d part_id %d rr %f \n", id, part_id, rr);
                                }
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
                                    fpair += force_sigma * wr * curand_normal_double(&state[id]) * 1 / sqrt(dt);
                                    f[id*3] += delx * fpair / rr;
                                    f[id*3+1] += dely * fpair / rr;
                                    f[id*3+2] += delz * fpair / rr; 
                                    //f[part_id*3] -= delx * fpair / rr;
                                    //f[part_id*3+1] -= dely * fpair / rr;
                                    //f[part_id*3+2] -= delz * fpair / rr;
                                    atomicAdd_30(&(f[part_id*3]), -delx * fpair /rr);
                                    atomicAdd_30(&(f[part_id*3+1]), -dely * fpair / rr);
                                    atomicAdd_30(&(f[part_id*3+2]), -delz * fpair / rr);
                               }
                            }
                        }
                //                if (id == 3999) 
                //                    printf("force of each cell is %f %f %f\n", f[id*3], f[id*3+1], f[id*3+2]);

                   }
                }
            }
   }
}


__global__ void post_int(double *v, double *f) {
     int id = threadIdx.x + blockDim.x * blockIdx.x;
    //printf("blockIdx.x %d threadIdx.x %d\n", blockIdx.x, threadIdx.x); 
    if (id < N) {
            //full integration
            int i;
            for (i = 0; i < 3; i++) {
                v[id * 3 + i] += 0.5 * f[id * 3 + i] * dt;
            }
 
 }   
}

void gpu_dpd(double** r, double**f, double** v, const char* select, RanMars* random, OUTPUT_struct* output_info, double len_cell, double len_x, double len_y, double len_z) {

        compute_force(r, v, f, random, N, output_info, len_cell, len_x, len_y, len_z);
        //allocate r,f,v on GPU
        int size = N * 3 * sizeof(double);
        double* cpu_r = (double *) malloc(size);
        double* cpu_f = (double *) malloc(size);
        double* cpu_v = (double *) malloc(size);
        int i;
        int count = 0;
        for(i = 0; i < N; i++) {
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
        if (count == 3 * N) {
            printf("count right\n");
            //return;
            double residual_r = 0;
            double residual_f = 0;
            double residual_v = 0;
            int j;
            for (j = 0; j < 3 * N; j++) {
                residual_r += (cpu_r[j] - r[j/3][j%3]);
                residual_f += (cpu_f[j] - f[j/3][j%3]);
                residual_v += (cpu_v[j] - v[j/3][j%3]);
            }
            printf("residual %f %f %f\n", residual_r, residual_f, residual_v);
        }
        double* dev_r;
        double* dev_f;
        double* dev_v;
        cudaMalloc((void **) &dev_r, size);
        cudaMalloc((void **) &dev_f, size);
        cudaMalloc((void **) &dev_v, size);
        cudaMemcpy(dev_r, cpu_r, size, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_f, cpu_f, size, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_v, cpu_v, size, cudaMemcpyHostToDevice);
        int** cell_partc_list_res = output_info->cell_partc_list_res;
        int* part_num = output_info->cell_partc_num_res;
        int cell_num = output_info->cell_num;
        int max_num_part = output_info->col_num;

        int* cell_list = (int *) malloc(sizeof(int) * cell_num * max_num_part);
        int* cell_list_count = (int *) malloc(sizeof(int) * cell_num);
        int j;
        count = 0;
        for (j = 0; j < cell_num; j++) {
            int k;
            for (k = 0; k < max_num_part; k++) {
                cell_list[count++] = cell_partc_list_res[j][k];
            }
            cell_list_count[j] = part_num[j]; 
        }
        double residual_count = 0;
        double residual_list = 0;
        for (j = 0; j < cell_num; j++) {
            residual_count += (cell_list_count[j] - part_num[j]);
            int k;
            for (k = 0; k < max_num_part; k++) {
                residual_list += (cell_list[j * max_num_part + k] - cell_partc_list_res[j][k]);
            }
        }

        printf("residual list %f count %f \n", residual_list, residual_count);
        int* dev_cell_list;
        int* dev_cell_list_count;
        int size_cell_list = sizeof(int) * cell_num * max_num_part;
        cudaMalloc((void **) &dev_cell_list, size_cell_list);
        cudaMalloc((void **) &dev_cell_list_count, sizeof(int) * cell_num);
        cudaMemcpy(dev_cell_list, cell_list, size_cell_list, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_cell_list_count, cell_list_count, sizeof(int) * cell_num, cudaMemcpyHostToDevice);
 
        curandState *d_state;
        cudaMalloc((void**) &d_state, N);
        int ntimestep = 5000;
        double m = 1.0;
        int blockSize = (int) floor(N/1 + 1);
        printf("go to call iteration\n");
        cout << "verify f " << verify_f(f) << endl;
        ofstream outputfile;
        outputfile.open("dump_gpu.md", ios::out);

        for (i = 0; i <= 5000; i++) {
            init<<<blockSize, 1>>>(dev_r, dev_v, dev_f);
            cudaDeviceSynchronize();
            f_clear<<<blockSize, 1>>>(dev_f);
            cudaDeviceSynchronize();
            iteration<<<blockSize, 1>>>(d_state, dev_r, dev_f, dev_v, dev_cell_list, dev_cell_list_count, max_num_part);
            cudaDeviceSynchronize();
            post_int<<<blockSize, 1>>>(dev_v, dev_f);
            cudaMemcpy(cpu_r, dev_r, size, cudaMemcpyDeviceToHost);
            cudaMemcpy(cpu_v, dev_v, size, cudaMemcpyDeviceToHost);
            cudaMemcpy(cpu_f, dev_f, size, cudaMemcpyDeviceToHost);           
            double res = 0;
            for (int j = 0; j < 3 * N; j++) {
                res += cpu_f[j];
            }
            //cout << " cpu_f is " << res << endl; 
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < 3; k++) {
                    r[j][k] = cpu_r[j * 3 + k];
                    v[j][k] = cpu_v[j * 3 + k];
                    f[j][k] = cpu_f[j * 3 + k];
                }
            }
            //rebuild cell list
            int partc_num_def = 4000;
            output_info =  func_partc_incell_stat(r, partc_num_def, len_cell, len_x, len_y, len_z);
            cell_partc_list_res = output_info->cell_partc_list_res;
            part_num = output_info->cell_partc_num_res;
            cell_num = output_info->cell_num;
            max_num_part = output_info->col_num;
            free(cell_list);
            free(cell_list_count);
            cell_list = (int *) malloc(sizeof(int) * cell_num * max_num_part);
            cell_list_count = (int *) malloc(sizeof(int) * cell_num);

            count = 0;
            for (int j = 0; j < cell_num; j++) {
                int k;
                for (k = 0; k < max_num_part; k++) {
                    cell_list[count++] = cell_partc_list_res[j][k];
                }
                cell_list_count[j] = part_num[j]; 
            }
            cudaFree(dev_cell_list);
            cudaFree(dev_cell_list_count);
            int size_cell_list = sizeof(int) * cell_num * max_num_part;
            cudaMalloc((void **) &dev_cell_list, size_cell_list);
            cudaMalloc((void **) &dev_cell_list_count, sizeof(int) * cell_num);

            cudaMemcpy(dev_cell_list, cell_list, size_cell_list, cudaMemcpyHostToDevice);
            cudaMemcpy(dev_cell_list_count, cell_list_count, sizeof(int) * cell_num, cudaMemcpyHostToDevice);
 
            if(i % 100 == 0) {
                double ke = computeKE(v);
                cout << i << " temp is " << ke * 2 / (3 * 4000 * 1) << endl;
                writeDump(outputfile, r, f, i);
                //cout << "verify force on gpu at step " << i << " is " << verify_f(f) << endl;

            }
        }
        printf("iteration done\n");
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
        cudaFree(d_state);
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

    int type_of_device = atoi(argv[1]); // 0 - CPU; 1 - GPU

    const char* select = argv[2];
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
        gpu_dpd(r, f, v, select, random, output_info,len_cell, len_x, len_y, len_z);  
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





