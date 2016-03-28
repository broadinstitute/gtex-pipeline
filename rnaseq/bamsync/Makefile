
CCFLAGS := -O3 -Wall -Wextra\
	-I/usr/local/include/bamtools

LDFLAGS := -L /usr/local/lib/bamtools

bamsync.o: bamsync.cpp
	g++ -std=c++11 $(CCFLAGS) bamsync.cpp $(LDFLAGS) -o bamsync -lboost_program_options -lbamtools
