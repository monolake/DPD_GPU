dpd.cpu : main.o sim_run.o calc_neighbor_list.o random_mars.o
	g++ -std=c++11 -o dpd.cpu main.o sim_run.o calc_neighbor_list.o random_mars.o
main.o : main.cpp
	g++ -g -std=c++11 -c main.cpp
sim_run.o : sim_run.cpp
	g++ -g -std=c++11 -c sim_run.cpp
calc_neighbor_list.o : calc_neighbor_list.cpp
	g++ -g -std=c++11 -c calc_neighbor_list.cpp
random_mars.o:	random_mars.cpp
	g++ -g -std=c++11 -c random_mars.cpp
clean: 
	rm dpd.cpu main.o sim_run.o calc_neighbor_list.o random_mars.o
