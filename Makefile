# makefile syntax for a rule:
# target: prerequisites
# <TAB> recipe

CODE = ./code/
BOOST_ROOT = /import/vldb01/2/scratch/mazhuo/boost/boost_1_77_0

# flags:
# -c 			- compile source code to object file
# -Ofast 		- highest level of optimization, longer compilation time
# -fopenmp		- 多核编程框架, link OpenMP runtime library
# -I/path		- path to header files
# -o 			- output executable file
# std=c++14		- use C++14 standard

all: wcsd
wcsd: main.o WGraph.o
	g++ -Ofast  -fopenmp  main.o WGraph.o -o wcsd -std=c++14
	rm *.o
main.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)main.cpp -I$(BOOST_ROOT) -std=c++14
WGraph.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)WGraph.cpp -I$(BOOST_ROOT) -std=c++14
