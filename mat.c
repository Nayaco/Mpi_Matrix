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

matPtr Strussen(const matPtr A, const matPtr B){
    if(A->c * A->r < 2500 && B->c * B->r < 2500)return MULT(A, B);
    matPtr C = CreateMat(A->r, B->c); 
    double *dataA = A->data;
    double *dataB = B->data;
    double *dataC = C->data;
    matPtr S[10];
    matPtr P[7];
    int ra = A->r / 2;
    int ca = A->c / 2;
    int rb = B->r / 2;
    int cb = B->c / 2;
    matPtr A11 = CreateMat(ra, ca);
    matPtr A12 = CreateMat(ra, ca);
    matPtr A21 = CreateMat(ra, ca);
    matPtr A22 = CreateMat(ra, ca);
    matPtr B11 = CreateMat(rb, cb);
    matPtr B12 = CreateMat(rb, cb);
    matPtr B21 = CreateMat(rb, cb);
    matPtr B22 = CreateMat(rb, cb);
    S[1] = CreateMat(ra, ca);
    S[2] = CreateMat(ra, ca);
    S[4] = CreateMat(ra, ca);
    S[6] = CreateMat(ra, ca);
    S[8] = CreateMat(ra, ca);
    S[0] = CreateMat(rb, cb);
    S[3] = CreateMat(rb, cb);
    S[5] = CreateMat(rb, cb);
    S[7] = CreateMat(rb, cb);
    S[9] = CreateMat(rb, cb);
    
    for(int i = 0; i < ra; i++){
        int ax1 = i * A->c;         double *a11 = dataA + ax1;
        int ax2 = ax1 + ca;         double *a12 = dataA + ax2;
        int ax3 = (i + ra) * A->c;  double *a21 = dataA + ax3;
        int ax4 = ax3 + ca;         double *a22 = dataA + ax4;
        int axx = i * ca;
        for(int j = 0; j < ca; j++){
            A11->data[axx + j] = a11[j];
            A12->data[axx + j] = a12[j];
            A21->data[axx + j] = a21[j];
            A22->data[axx + j] = a22[j];
            S[1]->data[axx + j] = a11[j] + a12[j];
            S[2]->data[axx + j] = a21[j] + a22[j];
            S[4]->data[axx + j] = a11[j] + a22[j];
            S[6]->data[axx + j] = a12[j] - a22[j];
            S[8]->data[axx + j] = a11[j] - a21[j];
        }
    }

    for(int i = 0; i < rb; i++){
        int bx1 = i * B->c;         double *b11 = dataB + bx1;
        int bx2 = bx1 + cb;         double *b12 = dataB + bx2;
        int bx3 = (i + rb) * B->c;  double *b21 = dataB + bx3;
        int bx4 = bx3 + cb;         double *b22 = dataB + bx4;
        int bxx = i * cb;
        for(int j = 0; j < cb; j++){
            B11->data[bxx + j] = b11[j];
            B12->data[bxx + j] = b12[j];
            B21->data[bxx + j] = b21[j];
            B22->data[bxx + j] = b22[j];
            S[0]->data[bxx + j] = b12[j] - b22[j];
            S[3]->data[bxx + j] = b21[j] - b11[j];
            S[5]->data[bxx + j] = b11[j] + b22[j];
            S[7]->data[bxx + j] = b21[j] + b22[j];
            S[9]->data[bxx + j] = b11[j] + b12[j];
        }
    }

    P[0] = Strussen(A11, S[0]);
    P[1] = Strussen(S[1], B22);
    P[2] = Strussen(S[2], B11);
    P[3] = Strussen(A22, S[3]);
    P[4] = Strussen(S[4], S[5]);
    P[5] = Strussen(S[6], S[7]);
    P[6] = Strussen(S[8], S[9]);

    for(int i = 0; i < ra; i++){
        int cx1 = i * C->c;         double *c11 = dataC + cx1;
        int cx2 = cx1 + cb;         double *c12 = dataC + cx2;
        int cx3 = (i + ra) * C->c;  double *c21 = dataC + cx3;
        int cx4 = cx3 + cb;         double *c22 = dataC + cx4;
        int bxx = i * cb;
        double *P1 = P[0]->data + bxx;
        double *P2 = P[1]->data + bxx;
        double *P3 = P[2]->data + bxx;
        double *P4 = P[3]->data + bxx;
        double *P5 = P[4]->data + bxx;
        double *P6 = P[5]->data + bxx;
        double *P7 = P[6]->data + bxx;
        for(int j = 0; j < cb; j++){
            c11[j] = P5[j] + P4[j] - P2[j] + P6[j];
            c12[j] = P1[j] + P2[j];
            c21[j] = P3[j] + P4[j];
            c22[j] = P5[j] + P1[j] - P3[j] - P7[j];
        }
    }

    if((A->c) & 1){
        int bbx = B->c * (B->r - 1);
        for(int i = 0; i < ra * 2; i++){
            int x = i * C->c;
            int aax = i * A->c + A->c - 1;
            for(int j = 0; j < cb * 2; j++){
                C->data[x + j] += A->data[aax] * B->data[bbx + j]; 
            }
        }
    }
    if((A->r) & 1){
        int aax = A->c * (A->r - 1);
        int ccx = C->c * (C->r - 1);
        for(int i = 0; i < cb * 2; i++)C->data[ccx + i] = 0;
        for(int i = 0; i < A->c; i++){
            int temp = A->data[aax + i];
            int bbx = i * A->c;
            for(int j = 0; j < cb * 2; j++){
                C->data[ccx + j] += temp * B->data[bbx + j]; 
            }
        }
    }
    if((B->c) & 1){
        for(int i = 1; i <= ra * 2; i++)C->data[i * C->c - 1] = 0;
        for(int i = 0; i < ra * 2; i++){
            int ccx = (i + 1) * C->c - 1;
            int aax = A->c * i;
            for(int j = 0; j < A->c; j++){
                C->data[ccx] += A->data[aax + j] * B->data[(j + 1) * B->c - 1]; 
            }
        }
    }
    if((B->c) & 1 && (A->r) & 1){
        int aax = A->c * (A->r - 1);
        int last = C->r * C->c - 1;
        C->data[last] = 0;
        for(int i = 0; i < A->c; i++)
            C->data[last] += A->data[aax + i] * B->data[(i + 1) * B->c - 1];
    }
    for(int i = 0; i < 10; i++)FreeMat(S[i]);
    for(int i = 0; i < 7; i++)FreeMat(P[i]);
    FreeMat(A11);
    FreeMat(A12);
    FreeMat(A21);
    FreeMat(A22);
    FreeMat(B11);
    FreeMat(B12);
    FreeMat(B21);
    FreeMat(B22);
    return C;
}
                                                        