#include	<stdio.h>
#include 	<stdlib.h>
#include 	<string.h>
#include	<math.h>
#include	<iostream>
#include	<fstream>
#include    "sim_run.h"
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

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);

	FILE*		file_ptr;
	
	char		str_input[5];
	double		result[3];
	int			count;
	int			end;
	int			flag;
	
	//Main:
	count = 0;
	if((file_ptr = fopen("4000_new.txt","r")) == NULL){
		printf("Cannt open the file!");
		exit(1);
	}
	
	next_func(file_ptr);
	
	//while (!feof(file_ptr)){
		//count = count+1;
		//flag = load_func(file_ptr, result);
		//if (flag){
		//	printf("Line:%d, %.2f, %.2f, %.2f\n", count, result[0], result[1], result[2]);
		//	
		//}
		
	//}

	fclose(file_ptr);
	//scanf("%d",&end);

    int N = 4000;
    double** r = new double* [N];
    double** v = new double* [N];
    double** f = new double* [N];
    int ntimestep = 0;


    for(int i = 0; i < N; i++) {
        r[i] = new double[3];
        v[i] = new double[3];
        f[i] = new double[3];
    }

    for (int i = 0; i < N; i++) {
        r[i][0] = result[0];
        r[i][1] = result[1];
        r[i][2] = result[2];
    }

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < 3; j++) {
            v[i][j] = 0;
            f[i][j] = 0;
        }
    }

    double rand_num = distribution(generator);
    //init force compute
    compute_force(r, v, f, rand_num, N);
    int ntimestep = 1000;
    for (int i = 0; i < ntimestep; i++) {

        //half integration
        for(int j = 0; j < N; j++) {
            v[j][0] += 0.5 * f[j][0] / m * dt;
            v[j][1] += 0.5 * f[j][1] / m * dt;
            v[j][2] += 0.5 * f[j][2] / m * dt;

            r[j][0] += v[j][0] * dt;
            r[j][1] += v[j][1] * dt;
            r[j][2] += v[j][2] * dt;
        }
 
        vector<vector<int>> cell_list;
        buildNeighborList(neighborlist, r);
        //force computation
        //update f
        clear_force(f);
        computeForce(r, v, f, distribution(generator), cell_list);

        //full integration
        for(int j = 0; j < N; j++) {
            v[j][0] += 0.5 * f[j][0] / m * dt;
            v[j][1] += 0.5 * f[j][1] / m * dt;
            v[j][2] += 0.5 * f[j][2] / m * dt;
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
    cout << "position is read, f is computed" << endl;
 
}





