dpd.cpu : load_main.o
	g++ -o dpd.cpu load_main.o
load_main.o : load_main.cpp
	g++ -c load_main.cpp
clean: 
	rm dpd.cpu load_main.o
