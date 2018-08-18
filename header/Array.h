#ifndef __Array_h__
#define __Array_h__

#include "../header/Argument.h"

#define __iterate__(v,st,ed) for((v)=(st);(v)<(ed);++(v))

// Array Index

#define Index5D(a,s1,s2,s3,s4,i0,i1,i2,i3,i4) ((a)[((((i0)*(s1)+(i1))*(s2)+(i2))*(s3)+(i3))*(s4)+(i4)])
#define Index4D(a,s1,s2,s3,i0,i1,i2,i3)       ((a)[ (((i0)*(s1)+(i1))*(s2)+(i2))*(s3)+(i3)])
#define Index3D(a,s1,s2,i0,i1,i2)             ((a)[  ((i0)*(s1)+(i1))*(s2)+(i2)])
#define Index2D(a,s1,i0,i1)                   ((a)[   (i0)*(s1)+(i1)])

// array of pointer to 1Darray

static inline void aop2a5DI(int *****p, int *a, int s0, int s1, int s2, int s3, int s4) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3) __iterate__(i4,0,s4)
    {
        Index5D(a,s1,s2,s3,s4,i0,i1,i2,i3,i4) = p[i0][i1][i2][i3][i4];
    }
}
static inline void aop2a4DI(int ****p, int *a, int s0, int s1, int s2, int s3) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3)
    {
        Index4D(a,s1,s2,s3,i0,i1,i2,i3) = p[i0][i1][i2][i3];
    }
}
static inline void aop2a3DI(int ***p, int *a, int s0, int s1, int s2) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2)
    {
        Index3D(a,s1,s2,i0,i1,i2) = p[i0][i1][i2];
    }
}
static inline void aop2a2DI(int **p, int *a, int s0, int s1) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1)
    {
        Index2D(a,s1,i0,i1) = p[i0][i1];
    }
}

static inline void aop2a5DF(Real *****p, Real *a, int s0, int s1, int s2, int s3, int s4) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3) __iterate__(i4,0,s4)
    {
        Index5D(a,s1,s2,s3,s4,i0,i1,i2,i3,i4) = p[i0][i1][i2][i3][i4];
    }
}
static inline void aop2a4DF(Real ****p, Real *a, int s0, int s1, int s2, int s3) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3)
    {
        Index4D(a,s1,s2,s3,i0,i1,i2,i3) = p[i0][i1][i2][i3];
    }
}
static inline void aop2a3DF(Real ***p, Real *a, int s0, int s1, int s2) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2)
    {
        Index3D(a,s1,s2,i0,i1,i2) = p[i0][i1][i2];
    }
}
static inline void aop2a2DF(Real **p, Real *a, int s0, int s1) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1)
    {
        Index2D(a,s1,i0,i1) = p[i0][i1];
    }
}

// 1Darray to pointer of array


static inline void a2aop5DI(int *a, int *****p, int s0, int s1, int s2, int s3, int s4) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3) __iterate__(i4,0,s4)
    {
        p[i0][i1][i2][i3][i4] = Index5D(a,s1,s2,s3,s4,i0,i1,i2,i3,i4);
    }
}
static inline void a2aop4DI(int *a, int ****p, int s0, int s1, int s2, int s3) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3)
    {
        p[i0][i1][i2][i3] = Index4D(a,s1,s2,s3,i0,i1,i2,i3);
    }
}
static inline void a2aop3DI(int *a, int ***p, int s0, int s1, int s2) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2)
    {
        p[i0][i1][i2] = Index3D(a,s1,s2,i0,i1,i2);
    }
}
static inline void a2aop2DI(int *a, int **p, int s0, int s1) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1)
    {
        p[i0][i1] = Index2D(a,s1,i0,i1);
    }
}

static inline void a2aop5DF(Real *a, Real *****p, int s0, int s1, int s2, int s3, int s4) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3) __iterate__(i4,0,s4)
    {
        p[i0][i1][i2][i3][i4] = Index5D(a,s1,s2,s3,s4,i0,i1,i2,i3,i4);
    }
}
static inline void a2aop4DF(Real *a, Real ****p, int s0, int s1, int s2, int s3) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2) \
    __iterate__(i3,0,s3)
    {
        p[i0][i1][i2][i3] = Index4D(a,s1,s2,s3,i0,i1,i2,i3);
    }
}
static inline void a2aop3DF(Real *a, Real ***p, int s0, int s1, int s2) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1) __iterate__(i2,0,s2)
    {
        p[i0][i1][i2] = Index3D(a,s1,s2,i0,i1,i2);
    }
}
static inline void a2aop2DF(Real *a, Real **p, int s0, int s1) {
    int i0,i1,i2,i3,i4;
    __iterate__(i0,0,s0) __iterate__(i1,0,s1)
    {
        p[i0][i1] = Index2D(a,s1,i0,i1);
    }
}

#endif
