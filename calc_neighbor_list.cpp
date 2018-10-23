#include "calc_neighbor_list.h"

/***********************************************************************************************/

//Function:		output_info = func_partc_incell_stat(input_cube, input_partc)

//				Type			Size					Name				Note
//input_cube:	cells information,
//				int				1						cell_num,
//				double			1						len_x, len_y, len_z,
//				double**		[cell_num*3]			cell_pos,			[pos_x, pos_y, pos_z]
//input_partc:	positions in cube,
//				int				1						partc_num,
//				double**		[partc_num*3]			partc_pos,			[pos_x, pos_y, pos_z]
//output_info:	each cell's partc_list
//				int				1						col_num
//				int*			[1*cell_num]			cell_partc_num
//				int**			[cell_num*col_num]		cell_partc_list

/***********************************************************************************************/


//const int		CONST_MAX_RAND = 1000;
//
//
//typedef struct OUTPUT_struct{
//	int				cell_num;
//	int				col_num; // return max #partices in a cell
//	int*			cell_partc_num_res; // return num of particles in each cell
//	int**			cell_partc_list_res; // return cell list
//	
//	int**			cell_id2loc;
//}OUTPUT_struct;


// function to calculate cell list 
OUTPUT_struct* func_partc_incell_stat(double** partc_pos, int partc_num, double len_cell, double len_x, double len_y, double len_z){
	
	int		num_cx, num_cy, num_cz;
	int*	partc_cell_id;
	int		i;
	int		p_cx, p_cy, p_cz;
	int		cell_num;
	int*	count_cell_id;
	int		max_partc_cell = 0;
	int**	cell_partc_list;
	int*	temp_cell_id;
	int		curr_cell_id;
	int		temp;
	int		cell_id;
	int**	cell_id2loc;
	OUTPUT_struct* result;
	
	// 1. initialization
	/* get #cells in x, y, z dimensions: num_cx, num_cy, num_cz */
	num_cx = (int) ceil(len_x/len_cell);
	num_cy = (int) ceil(len_y/len_cell);
	num_cz = (int) ceil(len_z/len_cell);

	// 2. loop over all N particles to get the list of cell that each particle belongs to
	partc_cell_id = (int*)malloc(partc_num*sizeof(int));
	
	for (i = 0; i<partc_num; i++){
		// get particle i's position according to cell:
		p_cx = (int) floor(partc_pos[i][0]/len_cell); 
		p_cy = (int) floor(partc_pos[i][1]/len_cell); 
		p_cz = (int) floor(partc_pos[i][2]/len_cell); 

		// get particle i's cell ID using cell index: 
		cell_id = p_cz + num_cz*p_cy + num_cy*num_cz*p_cx;
		
		partc_cell_id[i] = cell_id;
	}
	
	cell_num = num_cx*num_cy*num_cz;
	
	//cell_id mapping rules
	cell_id2loc = (int**)malloc(cell_num*sizeof(int*));
	for(i = 0; i<cell_num; i++){
		cell_id = i;
		cell_id2loc[i] = (int*)malloc(3*sizeof(int));
		
		
		cell_id2loc[i][0] = cell_id/(num_cy*num_cz);
		temp = cell_id%(num_cy*num_cz);
		cell_id2loc[i][1] = temp/num_cz;
		temp = temp%num_cz;
		cell_id2loc[i][2] = temp;
		
	}
	
	
	// 3. loop over the partc_cell_id array to get: #particles in each cell, max #particles in a cell
	count_cell_id = (int*)malloc(cell_num*sizeof(int));
	
	// initialize
	for (i=0; i<cell_num; i++){
		count_cell_id[i] = 0;
	}
	//
	for (i=0; i<partc_num; i++){
		count_cell_id[partc_cell_id[i]]++;
	}

	for (i=0; i<cell_num; i++){
		if (count_cell_id[i] > max_partc_cell){
			max_partc_cell = count_cell_id[i];
		}
	}

	// 4. create output cell list using: partc_cell_id[partc_num], count_cell_id[cell_num], max_partc_cell
	cell_partc_list = (int**)malloc(cell_num*sizeof(int*));
	for (i=0; i<cell_num; i++){
		cell_partc_list[i] = (int*)malloc(max_partc_cell*sizeof(int));
		
	}
	
	temp_cell_id = (int*)malloc(cell_num*sizeof(int));
	memcpy(temp_cell_id, count_cell_id, cell_num*sizeof(int));
	
	// get particles in each cell
	for (i=0; i<partc_num; i++){
		curr_cell_id = partc_cell_id[i]; // this means that ith particle belongs to curr_cell_id'th cell
		cell_partc_list[curr_cell_id][temp_cell_id[curr_cell_id]-1] = i;
		temp_cell_id[curr_cell_id]--;
	}
	
	// pack into output struct and return 
	result = (OUTPUT_struct*)malloc(sizeof(OUTPUT_struct));
	result->cell_num = cell_num;
	result->col_num = max_partc_cell;
	result->cell_partc_num_res = count_cell_id;
	result->cell_partc_list_res = cell_partc_list;
	result->cell_id2loc = cell_id2loc;
	
	
	free(temp_cell_id);
	free(partc_cell_id);
	return result;

}


void func_print_output_info(OUTPUT_struct*	output_info){
	
	int			loop, loop1;
	
	printf("Total cell_num is %d.\n\n", output_info->cell_num);
	
	printf("Particle information of each cell:\n");
	
	for(loop=0; loop<output_info->cell_num; loop++){
		printf("Cell[%d] loc.[%d][%d][%d] ------ [%d] num. of particles\n",
			loop,
			output_info->cell_id2loc[loop][0],
			output_info->cell_id2loc[loop][1],
			output_info->cell_id2loc[loop][2],
			output_info->cell_partc_num_res[loop]
			);
		
		if(output_info->cell_partc_num_res[loop]){
			printf("including: {");
			for(loop1=0; loop1<output_info->cell_partc_num_res[loop]; loop1++){
				
				printf("[%d]", output_info->cell_partc_list_res[loop][loop1]);
				
			}
			printf("}\n");
		}
		printf("\n");
		
	}
	
}


//int main(){
//	
//	int		loop;
//	unsigned int seed;
//	
//	//input_cube
//	double	len_cell;
//	double	len_x, len_y, len_z;
//	
//	//input_partc
//	int		partc_num;
//	double**	partc_pos;
//	
//	int		end_flag;
//	
//	//output
//	OUTPUT_struct*	output_info;
//	
//	seed = 10;
//	srand(seed);
//	
//	//INPUT information:
//	len_cell = 1.0;
//	len_x = 2.0;
//	len_y = 2.0;
//	len_z = 5.0;
//	
//	partc_num = 10;
//	partc_pos = (double**)malloc(partc_num*sizeof(double*));
//	for(loop=0; loop<partc_num; loop++){
//		partc_pos[loop] = (double*)malloc(3*sizeof(double));
//		partc_pos[loop][0] = (double)(rand()%CONST_MAX_RAND)/(CONST_MAX_RAND-1)*len_x;
//		partc_pos[loop][1] = (double)(rand()%CONST_MAX_RAND)/(CONST_MAX_RAND-1)*len_y;
//		partc_pos[loop][2] = (double)(rand()%CONST_MAX_RAND)/(CONST_MAX_RAND-1)*len_z;
//		
//	}
//	
//	
//	
//	output_info =  func_partc_incell_stat(partc_pos, partc_num, len_cell, len_x, len_y, len_z);
//	
//	
//	//PRINT basic information
//	//if (1){
//		
//		
//	//}
//	
//	//PRINT output information
//	if (1){
//		func_print_output_info(output_info);
//	}
//	
//	scanf("%d",&end_flag);
//	
//	
//	
//}
//



