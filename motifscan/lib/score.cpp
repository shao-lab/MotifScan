#include <iostream>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

const size_t N_BASE = 4;

typedef struct mat
{
    int n;
    double *a_arr;
    double *c_arr;
    double *g_arr;
    double *t_arr;
} MAT;

typedef struct motif_res
{
    int tarnum;
    double ratio;
    int *tarsite;
    double *tarratio;
} MOTIF_RES;

typedef struct peakNode
{
    int id;
    double value;
} PEAK_NODE;

typedef struct peakList
{
    int n;
    PEAK_NODE * peak;
} PEAK_LIST;

MAT *copyMAT(MAT *s)
{
    MAT *t = (MAT*)malloc(sizeof(MAT));

    t->n = s->n;
    t->a_arr = (double *)malloc(sizeof(double)*t->n);
    t->c_arr = (double *)malloc(sizeof(double)*t->n);
    t->g_arr = (double *)malloc(sizeof(double)*t->n);
    t->t_arr = (double *)malloc(sizeof(double)*t->n);
    memcpy(t->a_arr,s->a_arr,sizeof(double)*t->n);
    memcpy(t->c_arr,s->a_arr,sizeof(double)*t->n);
    memcpy(t->g_arr,s->a_arr,sizeof(double)*t->n);
    memcpy(t->t_arr,s->a_arr,sizeof(double)*t->n);
    return t;
}

void freeMAT(MAT *m)
{
    free(m->a_arr);
    free(m->c_arr);
    free(m->g_arr);
    free(m->t_arr);
    free(m);
    m = NULL;
}

extern "C" void freeMOTIF_RES(MOTIF_RES *mr)
{
    free(mr->tarsite);
    free(mr->tarratio);
    free(mr);
    mr = NULL;
}

MAT *expandMAT(double *B, size_t s)
{
    MAT *M = (MAT*)malloc(sizeof(MAT));
    M->n = s;

    M->a_arr = (double *)malloc(sizeof(double)*s);
    M->c_arr = (double *)malloc(sizeof(double)*s);
    M->g_arr = (double *)malloc(sizeof(double)*s);
    M->t_arr = (double *)malloc(sizeof(double)*s);
    int i = 0;
    for(i=0;i<s;i++)
    {
       M->a_arr[i] = B[0];
       M->c_arr[i] = B[1];
       M->g_arr[i] = B[2];
       M->t_arr[i] = B[3];
    }
  
    return M;
}

extern "C" MAT* logMAT(MAT *S)
{  
    int n = S->n;
    MAT *M = (MAT*)malloc(sizeof(MAT));
    M->a_arr = (double*)malloc(sizeof(double)*n);
    M->c_arr = (double*)malloc(sizeof(double)*n);
    M->g_arr = (double*)malloc(sizeof(double)*n);
    M->t_arr = (double*)malloc(sizeof(double)*n);
    M->n = n;
    double *a = (double*)malloc(sizeof(double)*n);
    double *c = (double*)malloc(sizeof(double)*n);
    double *g = (double*)malloc(sizeof(double)*n);
    double *t = (double*)malloc(sizeof(double)*n);
   
    int i = 0;
    for(i=0;i<n;i++)
    {
        a[i] = log(S->a_arr[i]);
        c[i] = log(S->c_arr[i]);
        g[i] = log(S->g_arr[i]);
        t[i] = log(S->t_arr[i]);
    }
    memcpy(M->a_arr, a, sizeof(double)*n);
    memcpy(M->c_arr, c, sizeof(double)*n);
    memcpy(M->g_arr, g, sizeof(double)*n);
    memcpy(M->t_arr, t, sizeof(double)*n);
    free(a);
    free(c);
    free(g);
    free(t);
    a = NULL;
    c = NULL;
    g = NULL;
    t = NULL;
    return M;
}

MAT* revMAT(MAT *S)
{
    int n = S->n;
    MAT *M = (MAT*)malloc(sizeof(MAT));
    M->a_arr = (double*)malloc(sizeof(double)*n);
    M->c_arr = (double*)malloc(sizeof(double)*n);
    M->g_arr = (double*)malloc(sizeof(double)*n);
    M->t_arr = (double*)malloc(sizeof(double)*n);
    M->n = n;
    double *a = (double*)malloc(sizeof(double)*n);
    double *c = (double*)malloc(sizeof(double)*n);
    double *g = (double*)malloc(sizeof(double)*n);
    double *t = (double*)malloc(sizeof(double)*n);
    int i = 0;
    for(i=0;i<n;i++)
    {
        a[i] = S->a_arr[n-i-1];
        c[i] = S->c_arr[n-i-1];
        g[i] = S->g_arr[n-i-1];
        t[i] = S->t_arr[n-i-1];
    }
    memcpy(M->a_arr, a, sizeof(double)*n);
    memcpy(M->c_arr, c, sizeof(double)*n);
    memcpy(M->g_arr, g, sizeof(double)*n);
    memcpy(M->t_arr, t, sizeof(double)*n);
    free(a);
    free(c);
    free(g);
    free(t);
    a = NULL;
    c = NULL;
    g = NULL;
    t = NULL;
    return M;
}

double * diff(double *A, int a, double *B, int b)
{  
   // assert(a==b);
    double *r = (double *)malloc(sizeof(double)*a);
    int i = 0;
    for(i=0;i<a;i++)
    {
        r[i] = A[i] - B[i];
    }
    return r;
}

double *sliding_score(double *w,  int n_w,  double *s,  int n_s)
{   
    int n_r = n_s-n_w+1;
    double *r = (double*)malloc(sizeof(double)*n_r);

    double sum = 0;
    int i = 0;
    for(i=0;i<n_r;i++)
    {  
        double z = 0;
        int j,k;
        for(j=i, k=0; k<n_w; j++,k++)
        {
            z += w[k]*s[j];
        };
        r[i] = z;
    }
    return r;
}

extern "C" double * sliding_score_mat(MAT *W, MAT *S)
{
    int n_w = W->n;
    int n_s = S->n;
    double * a_r_arr = sliding_score(W->a_arr,n_w,S->a_arr,n_s);
    double * c_r_arr = sliding_score(W->c_arr,n_w,S->c_arr,n_s);
    double * g_r_arr = sliding_score(W->g_arr,n_w,S->g_arr,n_s);
    double * t_r_arr = sliding_score(W->t_arr,n_w,S->t_arr,n_s);
    double * score_arr = (double*)malloc(sizeof(double)*(n_s-n_w+1));
    int i = 0;
    for(i=0;i<n_s-n_w+1;i++)
    {   
        score_arr[i] = a_r_arr[i]+c_r_arr[i]+g_r_arr[i]+t_r_arr[i];
    }
    free(a_r_arr);
    free(c_r_arr);
    free(g_r_arr);
    free(t_r_arr);
    return score_arr;
}

extern "C" double * sliding_score_mat_rev(MAT *W, MAT *S)
{
    int n_w = W->n;
    int n_s = S->n;
    MAT *W_r = revMAT(W);
    double * a_r_arr = sliding_score(W_r->t_arr,n_w,S->a_arr,n_s);
    double * c_r_arr = sliding_score(W_r->g_arr,n_w,S->c_arr,n_s);
    double * g_r_arr = sliding_score(W_r->c_arr,n_w,S->g_arr,n_s);
    double * t_r_arr = sliding_score(W_r->a_arr,n_w,S->t_arr,n_s);
    double *score_arr = (double*)malloc(sizeof(double)*(n_s-n_w+1));
    int i = 0;
    for(i=0;i<n_s-n_w+1;i++)
    {   
        score_arr[i] = a_r_arr[i]+c_r_arr[i]+g_r_arr[i]+t_r_arr[i];
    }
    free(a_r_arr);
    free(c_r_arr);
    free(g_r_arr);
    free(t_r_arr);
    freeMAT(W_r);
    return score_arr;
}

void printMAT(MAT* m)
{
    int len = m->n;
    int i = 0;
    printf("a: ");
    for(i=0;i<len;i++)
    {
        printf("%.3f ",m->a_arr[i]);
    }
    printf("\n");
    printf("c: ");
    for(i=0;i<len;i++)
    {
        printf("%.3f ",m->c_arr[i]);
    }
    printf("\n");
    printf("g: ");
    for(i=0;i<len;i++)
    {
        printf("%.3f ",m->g_arr[i]);
    }
    printf("\n");
    printf("t: ");
    for(i=0;i<len;i++)
    {
        printf("%.3f ",m->t_arr[i]);
    }
    printf("\n");

}

extern "C" MOTIF_RES * motif_scan_core(MAT *S, MAT *M, double* B, double max_score, double score_cutoff)
{
    //data prep.
    int motif_len = M->n;
    int seq_len = S->n;
    int ratio_len = seq_len - motif_len + 1;
    MAT *logM = logMAT(M);
    MAT *B_m = expandMAT(B,motif_len);
    MAT *logB = logMAT(B_m);
    double *by = sliding_score_mat(logB, S);
    // scan plus strand
    double *my = sliding_score_mat(logM, S);
    double *ratio_f = diff(my,ratio_len,by,ratio_len);

    // scan minus strand
    double *my_r = sliding_score_mat_rev(logM, S);
    double *ratio_r = diff(my_r,ratio_len,by,ratio_len);
    freeMAT(logM);
    freeMAT(B_m);
    freeMAT(logB);

    free(my);
    free(by);
    free(my_r);

    int tarnum = 0;
    double ratio = 0;
    double * allratio = (double*)malloc(sizeof(double)*ratio_len);
    int * idx = (int*)malloc(sizeof(int)*ratio_len);
    int i = 0;
    for(i = 0; i<ratio_len; i++)
    {   
        if(ratio_f[i]> ratio_r[i])
        {
            allratio[i] = ratio_f[i]/max_score;
        }else
        {
            allratio[i] = ratio_r[i]/max_score;
        }

        //filter out the target sites
        if(allratio[i]>=score_cutoff)
        {
            idx[i] = 1;
            tarnum++;
        }else
        {
            idx[i] = 0;
        }

        //extract the global ratio
        if(allratio[i]>ratio)
        {
            ratio = allratio[i];
        }
    }

    free(ratio_f);
    free(ratio_r);

    int * tarsite = (int*)malloc(sizeof(int)*tarnum);
    double * tarratio = (double*)malloc(sizeof(double)*tarnum);
    int j = 0;
    for(i=0,j=0; i<tarnum,j<ratio_len; j++)
    {
        if (idx[j] == 1)
        {
            tarsite[i] = j;
            tarratio[i] = allratio[j];
            i++;
        }
    }

    MOTIF_RES * mr = (MOTIF_RES *)malloc(sizeof(MOTIF_RES));
    mr->tarsite = (int *)malloc(sizeof(int)*tarnum);
    mr->tarratio = (double *)malloc(sizeof(double)*tarnum);
    mr->tarnum = tarnum;
    mr->ratio = ratio;
    memcpy(mr->tarsite,tarsite,sizeof(int)*tarnum);
    memcpy(mr->tarratio,tarratio,sizeof(double)*tarnum);
    
    free(allratio);
    free(tarratio);
    free(tarsite);
    free(idx);

    return mr;
}

extern "C" double motif_scan_core_simulation(MAT *S, MAT *M, double* B, double max_score)
{
    //data prep.
    int motif_len = M->n;
    int seq_len = S->n;
    int ratio_len = seq_len - motif_len + 1;
    MAT *logM = logMAT(M);
    MAT *B_m = expandMAT(B,motif_len);
    MAT *logB = logMAT(B_m);
    double *by = sliding_score_mat(logB, S);
    // scan plus strand
    double *my = sliding_score_mat(logM, S);
    double *ratio_f = diff(my,ratio_len,by,ratio_len);
    // scan minus strand
    double *my_r = sliding_score_mat_rev(logM, S);
    double *ratio_r = diff(my_r,ratio_len,by,ratio_len);
    freeMAT(logM);
    freeMAT(B_m);
    freeMAT(logB);
   
    free(my);
    free(by);
    free(my_r);

    double ratio = 0;
    if(ratio_f[0]> ratio_r[0])
    {
        ratio = ratio_f[0]/max_score;
    }else
    {
        ratio = ratio_r[0]/max_score;
    }
    
    free(ratio_f);
    free(ratio_r);

    return ratio;
}