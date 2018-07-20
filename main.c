#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "matMPi.h"
#include "mat.h"

int wise(matPtr T){
    for(int i = 0; i < T->r * T->c; i++){
        for(int j = 0; j < 1000; j++){
            T->data[i] += 2.333;
            T->data[i] /= 2.333;
            T->data[i] *= 1.414;
        }
    }
    return 1;
}

int main(int argc, char const *argv[])
{
    MPI_Init(NULL, NULL);
    connPtr conn = CreateConn(MPI_COMM_WORLD);
    Addfun(conn, wise);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    matPtr A;
    FILE *file;
    FILE *out1;
    if(rank == 0){
        file = fopen("test.in", "r");
        out1 = fopen("test.out1", "w");
        int n, m;
        fscanf(file, "%d %d", &n, &m);
        A = CreateMat(n, m);
        GetMat(A,file);
    }else{
        A = CreateMat(0, 0);
    }
    int s = clock();
    MWiseMpi(A , 0, conn, SYNC);
    printf("%d : %.4f\n", rank, (clock()-s )/1000000.0);
    if(rank == 0){
        PrintMat(A, out1);
        fclose(file);
        fclose(out1);
    }
    FreeConn(conn);
    MPI_Finalize();
    return 0;
}
