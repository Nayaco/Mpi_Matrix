CXX = gcc
CXXFLAG = -O0 -g -Wall
MPICXX = mpicc
MPICXXFLAG = -O0 -g -Wall

main: main.o matMPI.o mat.o
	$(MPICXX) $(MPICXXFLAG) main.o matMPI.o mat.o -o main
	rm main.o matMPI.o mat.o

main.o: main.c
	$(MPICXX) $(MPICXXFLAG) -c main.c -o main.o

matMPI.o: matMPI.c
	$(MPICXX) $(MPICXXFLAG) -c matMPI.c -o matMPI.o
	#$(CXX) $(CXXFLAG) mat.o tmp.o -o matMPI.o
	#rm tmp.o
	
mat.o: mat.c
	$(CXX) $(CXXFLAG) -c mat.c -o mat.o

gener: gener.c
	$(CXX) -O2 -Wall gener.c -o gener

clean:
	rm main gener



	