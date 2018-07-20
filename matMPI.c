#include "matMPi.h"

connPtr CreateConn(MPI_Comm comm){
    connPtr conn = calloc(1, sizeof(struct conn));
    conn->_comm = comm;
    conn->_wiseList = calloc(0, sizeof(mwise_callback));
    conn->list_len = conn->list_max = 0;
    return conn;
}

void FreeConn(connPtr conn){
    free(conn->_wiseList);
    free(conn);
}

void Addfun(connPtr conn, mwise_callback func){
    if(conn->list_max == conn->list_len){
        conn->list_max += 10;
        conn->_wiseList = realloc(conn->_wiseList, conn->list_max * sizeof(mwise_callback));
    }
    conn->_wiseList[conn->list_len] = func;
    conn->list_len++;
}

void synclizeMat(void* T, int size, MPI_Comm comm){
    MPI_Bcast((int*)T, size, MPI_INT, 0, comm);
}

//Plus//////////////////////////////////////////////////////////////////
matPtr AddMpi(matPtr A, matPtr B, connPtr conn, int sync){
    MPI_Comm comm = conn->_comm;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int np;
    MPI_Comm_size(comm, &np);
    
    matPtr T;
    if(rank == 0)T = A;
        else T = CreateMat(0, 0);
    
    if(sync == SYNC)synclizeMat(T, 2, comm);
        else if(rank > 0)T[0] = *A, T[1] = *B;

    matPtr C;
    double *sendbuffer = NULL;
    double *buffer = NULL;
    double *bufferb = NULL;
    double *bufferc = NULL;
    int *_size = NULL, *_disp = NULL;
    int len = (T->c * T->r) / np;

    if(rank == 0){ 
        C = CreateMat(T->r, T->c);
        int local_len = (T->c * T->r) - (np - 1) * len;
        _size = malloc(np * sizeof(int));
        _disp = malloc(np * sizeof(int));
        sendbuffer = malloc(T->r * T->c * 2 * sizeof(double));
        buffer = malloc(local_len * 2 * sizeof(double));
        bufferc = malloc(local_len * sizeof(double));

        _size[0] = local_len * 2;
        _disp[0] = 0;
        memcpy(sendbuffer, A->data, local_len * sizeof(double));
        memcpy(sendbuffer + local_len, B->data, local_len * sizeof(double));
        
        for(int i = 1; i < np; i++){
            _size[i] = len * 2;
            _disp[i] = _disp[i - 1] + _size[i - 1];
            memcpy(sendbuffer + _disp[i], A->data + len * (i - 1) + local_len, len * sizeof(double));
            memcpy(sendbuffer + _disp[i] + len, B->data + len * (i - 1) + local_len, len * sizeof(double));
        }
        
        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, _size[0], MPI_DOUBLE, 0, comm);
        bufferb = buffer + local_len; 
        for(int i = 0; i < local_len; i++)
            bufferc[i] = buffer[i] + bufferb[i];

        for(int i = 0; i < np; i++)
            _size[i] /= 2, _disp[i] /= 2;

        MPI_Gatherv(bufferc, _size[0], MPI_DOUBLE, C->data, _size, _disp, MPI_DOUBLE, 0, comm);
        
        free(_size);
        free(_disp);
        free(sendbuffer);
        free(buffer);
        free(bufferc);
        return C;
    }else{

        C = CreateMat(0, 0);
        buffer = malloc(len * 2 * sizeof(double));
        bufferc = malloc(len * sizeof(double));
        
        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, len * 2, MPI_DOUBLE, 0, comm);
        
        bufferb = buffer + len;
        for(int i = 0; i < len; i++)
            bufferc[i] = buffer[i] + bufferb[i];
        
        MPI_Gatherv(bufferc, len, MPI_DOUBLE, C->data, _size, _disp, MPI_DOUBLE, 0, comm);
        
        FreeMat(T);
        free(buffer);
        free(bufferc);
        return NULL;
    }
}

//Mult//////////////////////////////////////////////////////////////////
matPtr MultMpi(matPtr A, matPtr B, connPtr conn, int sync){
    MPI_Comm comm = conn->_comm;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int np;
    MPI_Comm_size(comm, &np);
    
    rawmat T[2];
    if(rank == 0)T[0] = *A, T[1] = *B;
    
    if(sync == SYNC)synclizeMat(T, 2*sizeof(rawmat) / sizeof(int), comm);
        else if(rank > 0)T[0] = *A, T[1] = *B;

    if(T[0].c != T[1].r)return NULL;
    
    matPtr temp;
    rawmat local_a;
    rawmat local_b;
    matPtr S;
    matPtr C;
    double *buffer = NULL;
    double *sendbuffer = NULL;

    int *_size = NULL, *_disp = NULL;
    int lena = (T[0].c / np) * T[0].r;
    int lenb = (T[1].r / np) * T[1].c;

    if(rank == 0){
        temp = Trans(A);
        C = CreateMat(T[0].r, T[1].c);

        int local_lena = (T[0].r * T[0].c) - (np - 1) * lena;
        int local_lenb = (T[1].r * T[1].c) - (np - 1) * lenb;

        buffer = malloc((local_lena + local_lenb) * sizeof(double));
        sendbuffer = malloc((T[0].c * T[0].r + T[1].c * T[1].r) * sizeof(double));

        _size = malloc(np * sizeof(int));
        _disp = malloc(np * sizeof(int));
        _size[0] = local_lena + local_lenb;
        _disp[0] = 0;
        memcpy(sendbuffer, temp->data, local_lena * sizeof(double));
        memcpy(sendbuffer + local_lena, B->data, local_lenb * sizeof(double));
        for(int i = 1; i < np; i++){
            _size[i] = lena + lenb;
            _disp[i] = _disp[i - 1] + _size[i - 1];
            memcpy(sendbuffer + _disp[i], temp->data + local_lena + lena * (i - 1), lena * sizeof(double));
            memcpy(sendbuffer + _disp[i] + lena, B->data + local_lenb + lenb * (i - 1), lenb * sizeof(double));
        }

        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, _size[0], MPI_DOUBLE, 0, comm);
        
        local_a.r = T[0].c - (T[0].c / np) * (np - 1);
        local_a.c = T[0].r;
        local_b.r = T[1].r - (T[1].r / np) * (np - 1);
        local_b.c = T[1].c;

        local_a.data = buffer;
        local_b.data = buffer + local_lena;

        Tranself(&local_a);
        S = MULT(&local_a, &local_b);
        MPI_Reduce(
            S->data,
            C->data,
            T[0].r * T[1].c,
            MPI_DOUBLE,
            MPI_SUM,
            0,
            comm);
        
        free(_size);
        free(_disp);
        FreeMat(S);
        FreeMat(temp);
        free(buffer);
        free(sendbuffer);
        return C;
    }else{
        C = CreateMat(0, 0);
        buffer = malloc((lena + lenb) * sizeof(double));

        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, lena + lenb, MPI_DOUBLE, 0, comm);
        
        local_a.r = T[0].c / np;
        local_a.c = T[0].r;
        local_b.r = T[1].r / np;
        local_b.c = T[1].c;

        local_a.data = buffer;
        local_b.data = buffer + lena;

        Tranself(&local_a);       
        S = MULT(&local_a, &local_b);
        
        MPI_Reduce(
            S->data,
            C->data,
            T[0].r * T[1].c,
            MPI_DOUBLE,
            MPI_SUM,
            0,
            comm);
           
        FreeMat(S);
        free(buffer);
        return NULL;
    }
}

//Wise//////////////////////////////////////////////////////////////////
int MWiseMpi(matPtr A, uint64_t func, connPtr conn, int sync){
    MPI_Comm comm = conn->_comm;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int np;
    MPI_Comm_size(comm, &np);
    
    rawmat T;

    if(rank == 0)T = *A, T.data = (double*)func;

    if(sync == SYNC)synclizeMat(&T, sizeof(rawmat) / sizeof(int), comm);
        else if(rank > 0)T = *A;
    mwise_callback tfunc = conn->_wiseList[(uint64_t)T.data];
    
    rawmat local_mat;
    double *buffer = NULL;
    double *sendbuffer = NULL;

    int *_size = NULL, *_disp = NULL;
    int len = (T.r / np) * T.c;

    if(rank == 0){
        int local_len = (T.r * T.c) - (np - 1) * len;

        buffer = malloc((local_len) * sizeof(double));
        sendbuffer = A->data;

        _size = malloc(np * sizeof(int));
        _disp = malloc(np * sizeof(int));
        _size[0] = local_len;
        _disp[0] = 0;

        for(int i = 1; i < np; i++){
            _size[i] = len;
            _disp[i] = _disp[i - 1] + _size[i - 1];
        }

        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, local_len, MPI_DOUBLE, 0, comm);
        
        T.r = local_len / T.c;
        T.data = buffer;
        int ans = MWISE(&T, tfunc);

        MPI_Gatherv(T.data, local_len, MPI_DOUBLE, A->data, _size, _disp, MPI_DOUBLE, 0, comm);        

        free(_size);
        free(_disp);
        free(buffer);
        return ans;
    }else{
        buffer = malloc((len) * sizeof(double));
        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, len, MPI_DOUBLE, 0, comm);
        
        T.r = len / T.c;    
        T.data = buffer;
        int ans = MWISE(&T, tfunc);
        MPI_Gatherv(T.data, len, MPI_DOUBLE, A->data, _size, _disp, MPI_DOUBLE, 0, comm);        
        
        free(buffer);
        return ans;
    }
}

//strussen/////////////////////////////////////////////////
matPtr StrussenMpi(matPtr A, matPtr B, connPtr conn, int sync){
    MPI_Comm comm = conn->_comm;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int np;
    MPI_Comm_size(comm, &np);
    
    rawmat T[2];
    if(rank == 0)T[0] = *A, T[1] = *B;
    
    if(sync == SYNC)synclizeMat(T, 2*sizeof(rawmat) / sizeof(int), comm);
        else if(rank > 0)T[0] = *A, T[1] = *B;

    if(T[0].c != T[1].r)return NULL;
    
    matPtr temp;
    rawmat local_a;
    rawmat local_b;
    matPtr S;
    matPtr C;
    double *buffer = NULL;
    double *sendbuffer = NULL;

    int *_size = NULL, *_disp = NULL;
    int lena = (T[0].c / np) * T[0].r;
    int lenb = (T[1].r / np) * T[1].c;

    if(rank == 0){
        temp = Trans(A);
        C = CreateMat(T[0].r, T[1].c);

        int local_lena = (T[0].r * T[0].c) - (np - 1) * lena;
        int local_lenb = (T[1].r * T[1].c) - (np - 1) * lenb;

        buffer = malloc((local_lena + local_lenb) * sizeof(double));
        sendbuffer = malloc((T[0].c * T[0].r + T[1].c * T[1].r) * sizeof(double));

        _size = malloc(np * sizeof(int));
        _disp = malloc(np * sizeof(int));
        _size[0] = local_lena + local_lenb;
        _disp[0] = 0;
        memcpy(sendbuffer, temp->data, local_lena * sizeof(double));
        memcpy(sendbuffer + local_lena, B->data, local_lenb * sizeof(double));
        for(int i = 1; i < np; i++){
            _size[i] = lena + lenb;
            _disp[i] = _disp[i - 1] + _size[i - 1];
            memcpy(sendbuffer + _disp[i], temp->data + local_lena + lena * (i - 1), lena * sizeof(double));
            memcpy(sendbuffer + _disp[i] + lena, B->data + local_lenb + lenb * (i - 1), lenb * sizeof(double));
        }

        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, _size[0], MPI_DOUBLE, 0, comm);
        
        local_a.r = T[0].c - (T[0].c / np) * (np - 1);
        local_a.c = T[0].r;
        local_b.r = T[1].r - (T[1].r / np) * (np - 1);
        local_b.c = T[1].c;

        local_a.data = buffer;
        local_b.data = buffer + local_lena;

        Tranself(&local_a);
        S = MULT(&local_a, &local_b);
        MPI_Reduce(
            S->data,
            C->data,
            T[0].r * T[1].c,
            MPI_DOUBLE,
            MPI_SUM,
            0,
            comm);
        
        free(_size);
        free(_disp);
        FreeMat(S);
        FreeMat(temp);
        free(buffer);
        free(sendbuffer);
        return C;
    }else{
        C = CreateMat(0, 0);
        buffer = malloc((lena + lenb) * sizeof(double));

        MPI_Scatterv(sendbuffer, _size, _disp, MPI_DOUBLE, buffer, lena + lenb, MPI_DOUBLE, 0, comm);
        
        local_a.r = T[0].c / np;
        local_a.c = T[0].r;
        local_b.r = T[1].r / np;
        local_b.c = T[1].c;

        local_a.data = buffer;
        local_b.data = buffer + lena;

        Tranself(&local_a);       
        S = Strussen(&local_a, &local_b);
        
        MPI_Reduce(
            S->data,
            C->data,
            T[0].r * T[1].c,
            MPI_DOUBLE,
            MPI_SUM,
            0,
            comm);
           
        FreeMat(S);
        free(buffer);
        return NULL;
    }
}