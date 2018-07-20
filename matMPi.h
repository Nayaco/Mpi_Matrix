#ifndef  MATMPI_H
#define MATMPI_H
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#include "mat.h"

#define NOSYNC          0
#define SYNC            1

#define MIN_LEN_P       131072

struct conn{
    MPI_Comm _comm;
    mwise_callback *_wiseList;
    int list_len;
    int list_max;
};
typedef struct conn *connPtr;

connPtr CreateConn(MPI_Comm comm);

void FreeConn(connPtr conn);

void Addfun(connPtr conn, mwise_callback func);

void synclizeMat(void* T, int size, MPI_Comm comm);

matPtr AddMpi(matPtr A, matPtr B, connPtr conn, int sync);

matPtr MultMpi(matPtr A, matPtr B, connPtr conn, int sync);

int MWiseMpi(matPtr A, uint64_t func, connPtr conn, int sync);

#endif