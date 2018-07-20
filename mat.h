#ifndef MAT_H
#define MAT_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef int (*mwise_callback)(matPtr); 

typedef struct MAT{
    int r, c;
    double *data;
}rawmat, *matPtr;

matPtr CreateMat(int r, int c);

void FreeMat(matPtr T);

int GetMat(matPtr T, FILE *stream);

void PrintMat(const matPtr T, FILE *stream);

matPtr Trans(const matPtr T);

void Tranself(matPtr T);

matPtr ADD(const matPtr A,const matPtr B);

matPtr MULT(const matPtr A, const matPtr B);

int MWISE(matPtr T, mwise_callback func);

int SPLIT_R(const matPtr T, int r, matPtr Dest_1, matPtr Dest_2);

int SPLIT_C(const matPtr T, int c, matPtr Dest_1, matPtr Dest_2);

matPtr Strussen(const matPtr A, const matPtr B);

#endif