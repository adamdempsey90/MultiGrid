#include "stdlib.h"
#include "stdio.h"
void solve(double *u, double *f) {
    /* Exact solution for poisson problem on a 3x3 grid where,
     * all points are boundary points except for a single 
     * interior point.
     * If level > 0 then we make no modifications to u
     */
    int i,indx;
    double h, h2;
    int n = 3;

    h = .5;
    h2 = h*h;

    for(i=0;i<n*n;i++) u[i] = 0;

    indx = 1 + n;
    u[indx] = -h2*f[indx]*.25;
    return;
}
