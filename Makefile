dpd.cpu : load_main.o sim_run.o
	g++ -std=c++11 -o dpd.cpu load_main.o sim_run.o
load_main.o : load_main.cpp
	g++ -std=c++11 -c load_main.cpp
sim_run.o : sim_run.cpp
	g++ -std=c++11 -c sim_run.cpp
clean: 
	rm dpd.cpu load_main.o sim_run.o
