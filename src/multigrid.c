#include "stdlib.h"

extern void matfill(double *,  double *,  int );
extern void restrict2D(double *,  double *,  int );
extern void solve(double *, double *, int);

void multigrid(double *u, double *f, const int lmin, const int lmax, const int numiter) {
/* Multigrid V-cycle for the Poisson equation 
 * u and f are given on the lmax level and are
 * coarsened to the lmin level
 */

    int i,j,k;
    int n;
    int ntop = i + (2 << lmax);
    int nbot = i + (2 << lmin);
    int nlevels = lmax - lmin + 1;


/* Allocate the storage for the grids */
    int *ngrid; 
    double **sol, **rhs;
    ngrid = (int *)malloc(sizeof(int) * nlevels);
    sol = (double **)malloc(sizeof(double *) * nlevels);
    rhs = (double **)malloc(sizeof(double *) * nlevels);
    for(i=lmin; i<=lmax;i++) {
        n = 1 + (2 << i);
        ngrid[i-lmin] = n;
        sol[i-lmin] = (double *)malloc(sizeof(double)*n*n);
        rhs[i-lmin] = (double *)malloc(sizeof(double)*n*n);
    }

/* Copy initial values to finest grid */

    matfill(sol[lmax-lmin], u, ngrid[lmax-lmin]);
    matfill(rhs[lmax-lmin], f, ngrid[lmax-lmin]);

/* Restrict to coarser meshes */
    for(i=lmax-1; i>=lmin; i--) {
        restrict2D(sol[i-lmin], sol[i+1-lmin], i);
        restrict2D(sol[i-lmin], sol[i+1-lmin], i);
    }

/* Solve on lowest level */

    solve(sol[0],rhs[0], 0);

    

/* Free grid storage */
    for(i=lmin; i<=lmax;i++) {
        free(sol[i-lmin]); 
        free(rhs[i-lmin]); 
    }
    free(sol);
    free(rhs);
    free(ngrid);

    return;

}
