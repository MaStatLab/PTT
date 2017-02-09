%CC= g++ -Wall -g
CC = g++ -O3 


all: markov_apt_density_C.so # optree_test test

markov_apt_density_C.so: markov_apt_kD_density_C.cpp
	R CMD SHLIB markov_apt_kD_density_C.cpp


clean:
	rm *.o *.so
