CXX = gcc
CXXFLAG = -O2 -Wall
MPICXX = mpicc
MPICXXFLAG = -O2 -Wall

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

tester: tester.o mat.o
	$(CXX) $(CXXFLAG) tester.o mat.o -o tester

tester.o: tester.c
	$(CXX) $(CXXFLAG) -c tester.c -o tester.o


clean:
	rm main gener



	