#include "mat.h"

matPtr CreateMat(int r, int c){
    matPtr T = malloc(sizeof(rawmat));
    T->c = c;
    T->r = r;
    T->data = calloc(r * c, sizeof(double));
    return T;
}

void FreeMat(matPtr T){
    free(T->data);
    free(T);
}

int GetMat(matPtr T, FILE *stream){
    int len = T->r * T->c;
    double *_data = T->data;
    for(int i = 0; i < len; i++){
        if(fscanf(stream, "%lf", _data + i) == -1)return 1;
    }
    return 0;
}

void PrintMat(const matPtr T, FILE *stream){
    int len = T->r * T->c;
    double *_data = T->data;
    for(int i = 0; i < len; i++){
        if(i % T->c == 0 && i)fprintf(stream, "\n");
        fprintf(stream, "%.0lf ", _data[i]);
    }
}

void Tranself(matPtr T){
    int tr = T->r;
    int tc = T->c;
    double *buffer = malloc(tr * tc * sizeof(double));
    for(int i = 0; i < T->r; i++){
        int tx = tc * i;
        for(int j = 0; j < T->c; j++)
            buffer[i + tr * j] = (T->data)[tx + j]; 
    }
    T->c = tr;
    T->r = tc;
    memcpy(T->data, buffer, tr * tc * sizeof(double));
    free(buffer);
    return;
}

matPtr Trans(const matPtr T){
    matPtr C = CreateMat(T->c, T->r);
    int tr = T->r;
    int tc = T->c;
    for(int i = 0; i < T->r; i++){
        int tx = tc * i;
        for(int j = 0; j < T->c; j++)
            (C->data)[i + tr * j] = (T->data)[tx + j]; 
    }
    return C;
}

matPtr ADD(const matPtr A, const matPtr B){
    if(A->r != B->r || A->c != B->c)return NULL;
    int len = A->r * A->c;
    matPtr T = CreateMat(A->r, A->c);
    double *_data = T->data;
    double *_dataA = A->data;
    double *_dataB = B->data;
    for(int i = 0; i < len; i++)_data[i] = _dataA[i] + _dataB[i];
    return T;
}

matPtr MULT(const matPtr A, const matPtr B){
    if(A->c != B->r )return NULL;
    matPtr T = CreateMat(A->r, B->c);
    int len = A->c;
    int ar = A->r;
    int bc = B->c;
    double *_data = T->data;
    double *_dataA = A->data;
    double *_dataB = B->data;

    for(int i = 0; i < ar; i++){
        int ax = i * len;
        int x = i * bc;
        for(int k = 0; k < len; k++){
            int temp = _dataA[ax + k];
            int bx = k * bc; 
            for(int j = 0; j < bc; j++){
                _data[x + j] += temp * _dataB[bx + j]; 
            }
        }
    }
    return T;
}

int MWISE(matPtr T, int (*func)(matPtr)){
    return func(T);
}

int SLICE_R(const matPtr T, int r, matPtr Dest_1, matPtr Dest_2){
    if(Dest_1 == NULL || Dest_2 == NULL)return 1;
    int size_1 = T->c * r;
    int size_2 = T->c * (T->r - r);
    Dest_1->data = realloc(Dest_1->data, size_1 * sizeof(double));
    Dest_2->data = realloc(Dest_2->data, size_2 * sizeof(double));
    Dest_1->r = r;
    Dest_1->c = T->c;
    Dest_2->r = T->r - r;
    Dest_2->c = T->c;
    memcpy(Dest_1->data, T->data, size_1 * sizeof(double));
    memcpy(Dest_2->data, T->data + size_1, size_2 * sizeof(double));
    return 0;
}

int SLICE_C(const matPtr T, int c, matPtr Dest_1, matPtr Dest_2){
    if(Dest_1 == NULL || Dest_2 == NULL)return 1;
    int size_1 = T->r * c;
    int size_2 = T->r * (T->c - c);
    int size = T->r * T->c;
    int Tc = T->c;
    int Tr = T->r;
    Dest_1->data = realloc(Dest_1->data, size_1 * sizeof(double));
    Dest_2->data = realloc(Dest_2->data, size_2 * sizeof(double));
    Dest_1->r = T->r;
    Dest_1->c = c;
    Dest_2->r = T->r;
    Dest_2->c = T->c - c;
    for(int i =  0; i < Tr; i++){
        int x_1 = i * c;
        int x_2 = i * (Tc - c);
        int x = i * Tc;
        memcpy(Dest_1->data + x_1, T->data + x, c * sizeof(double));
        memcpy(Dest_2->data + x_2, T->data + x +c, (Tc - c) * sizeof(double));
    }
    return 0;
}
                                                                         