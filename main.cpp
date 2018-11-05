#include	<stdio.h>
#include 	<stdlib.h>
#include 	<string.h>
#include	<math.h>
#include	<iostream>
#include	<fstream>
#include    "sim_run.h"
#include    "calc_neighbor_list.h"
#include    "random_mars.h"

using namespace std;

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

int main()
{

    //unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine generator (seed1);
    //std::normal_distribution<double> distribution (0.0,(double) 1/sqrt(dt));
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
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed1);
    std::normal_distribution<double> distribution (0.0,(double) 1/sqrt(dt));

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
    ofstream outputfile;
    outputfile.open("dump.md", ios::out);

    //double rand_num = distribution(generator);
    //init force compute
    //compute_force(r, v, f, random, N, output_info);
    double rand_num = 0;
    compute_force_std(r, v, f, generator, distribution, random, N);
    writeDump(outputfile, r, v, 0);

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
        //vector<vector<int>> cell_list;
        //buildNeighborList(neighborlist, r);
        //force computation
        //update f
        clear_force(f, N);
        //rand_num = distribution(generator);
        //compute_force(r, v, f, random, N, output_info);
        compute_force_std(r, v, f, generator, distribution, random, N);

        //full integration
        for(int j = 0; j < N; j++) {
            v[j][0] += 0.5 * f[j][0] * dt;
            v[j][1] += 0.5 * f[j][1] * dt;
            v[j][2] += 0.5 * f[j][2] * dt;
        }
         if(i % 100 == 0) {
            double ke = computeKE(v);
            cout << i << " temp is " << ke * 2 / (3 * 4000 * 1) << endl;
            //writeDump(outputfile, r, v, i);
        }
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





