#include "stdlib.h"
extern void matfill(double *, double *, int);
void solve(double *u, double *f, const int level) {
    /* Exact solution for poisson problem on a 3x3 grid where,
     * all points are boundary points except for a single 
     * interior point.
     * If level > 0 then we make no modifications to u
     */
    int i,j,indx,n;
    double h,h2;
    n = 1 + (2 << level);

    if (n != 3) {
        return;
    }

    h = 1./(n-1);
    h2 = h*h;

    i = 1; j = 1;

    indx = j + n*i;
    matfill(u,NULL,n);
    u[indx] = .25*(u[indx+n] + u[indx-n] + u[indx+1] + u[indx-1] - h2*f[indx]); 
    return;
}
