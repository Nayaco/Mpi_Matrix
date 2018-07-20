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
    FILE *out2;
    FILE *out3;
    FILE *out4;
    int s;
    matPtr C;
    if(rank == 0){
        file = fopen("test.in", "r");
        out1 = fopen("test.out1", "w");
        out2 = fopen("test.out2", "w");
        out3 = fopen("test.out3", "w");
        out4 = fopen("test.out4", "w");
        int n, m;
        fscanf(file, "%d %d", &n, &m);
        A = CreateMat(n, m);
        GetMat(A,file);
        s = clock();
    }else{
        A = CreateMat(0, 0);
    }
    
    C = AddMpi(A , A, conn, SYNC);
    if(rank == 0){
        printf("%.4f\n", (clock()-s )/1000000.0);
        //PrintMat(C, out1);
        FreeMat(C);
        s = clock();
    }

    C = MultMpi(A , A, conn, SYNC);
    if(rank == 0){
        printf("%.4f\n", (clock()-s )/1000000.0);
        //PrintMat(C, out2);
        FreeMat(C);
        s = clock();
    }
    
    C = StrussenMpi(A , A, conn, SYNC);
    if(rank == 0){
        printf("%.4f\n", (clock()-s )/1000000.0);
        //PrintMat(C, out3);
        FreeMat(C);
        s = clock();
    }

    MWiseMpi(A, 0, conn, SYNC);
    if(rank == 0){
        printf("%.4f\n", (clock()-s )/1000000.0);
        //PrintMat(A, out4);
    }

    if(rank == 0){
        fclose(file);
        fclose(out1);
        fclose(out2);
        fclose(out3);
        fclose(out4);
    }

    FreeConn(conn);
    MPI_Finalize();
    return 0;
}
