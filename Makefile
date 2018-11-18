dpd.cpu: main.o sim_run.o calc_neighbor_list.o random_mars.o
	nvcc -o dpd.cpu main.o sim_run.o calc_neighbor_list.o random_mars.o
main.o: main.cu 
	nvcc --resource-usage --std=c++11 -c main.cu
sim_run.o: sim_run.cpp sim_run.h
	g++  --std=c++11 -c sim_run.cpp
calc_neighbor_list.o: calc_neighbor_list.cpp calc_neighbor_list.h
	g++  --std=c++11 -c calc_neighbor_list.cpp
random_mars.o: random_mars.cpp random_mars.h
	g++  --std=c++11 -c random_mars.cpp
clean: 
	rm dpd.cpu main.o sim_run.o calc_neighbor_list.o random_mars.o
