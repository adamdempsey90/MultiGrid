#include "stdlib.h"
#include "stdio.h"
void solve(double *u, double *f) {
    /* Exact solution for poisson problem on a 3x3 grid where,
     * all points are boundary points except for a single 
     * interior point.
     * If level > 0 then we make no modifications to u
     */
    int i,j,indx;
    double h, h2;
    int n = 3;
    printf("SOLVE\n");

    h = .5;
    h2 = h*h;

    for(i=0;i<n*n;i++) u[i] = 0;

    indx = 1 + n;
    u[indx] = -h2*f[indx]*.25;
    printf("RHS\n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            printf("%lg\t",f[j+i*n]);
        }
        printf("\n");
    }

    printf("SOL\n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            printf("%lg\t",u[j+i*n]);
        }
        printf("\n");
    }
    return;
}
