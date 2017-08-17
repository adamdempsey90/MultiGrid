#include "stdlib.h"
#include "stdio.h"


extern void matfill(double *,  double *,  int );
extern void restrict2D(double *,  double *,  int );
extern void solve(double *, double *);
extern void gauss_seidel(double *, double *, int, int);
extern void prolongate2D(double *, double *, int );
extern void residual(double *, double *, double *, int);
extern void vcycle(double **, double *,int, int,int);



void multigrid(double *u, double *f, const int lmax, const int numcycles, const int numpre, const int numpost) {
/* Multigrid V-cycle for the Poisson equation 
 * u and f are given on the lmax level and are
 * coarsened to the lmin level
 */

    int i,l;
    int n;
    int nlevels = lmax + 1;


/* Allocate the storage for the grids */
    int *ngrid; 
    double **sol, **rhs, **res;
    ngrid = (int *)malloc(sizeof(int) * nlevels);
    sol = (double **)malloc(sizeof(double *) * nlevels);
    rhs = (double **)malloc(sizeof(double *) * nlevels);
    res = (double **)malloc(sizeof(double *) * nlevels);
    for(i=0; i<=lmax;i++) {
        n = 1 + (2 << i);
        ngrid[i] = n*n;
        sol[i] = (double *)malloc(sizeof(double)*n*n);
        rhs[i] = (double *)malloc(sizeof(double)*n*n);
        res[i] = (double *)malloc(sizeof(double)*n*n);
    }


/* Copy initial values to finest grid */

    matfill(sol[lmax], u, ngrid[lmax]);
    matfill(rhs[lmax], f, ngrid[lmax]);
    matfill(res[lmax], NULL, ngrid[lmax]);



/* Restrict to coarser meshes */
    for(i=lmax; i>0; i--) {
        restrict2D(sol[i-1], sol[i], i-1);
        restrict2D(rhs[i-1], rhs[i], i-1);
        matfill(res[i], NULL, ngrid[i]);
    }

/* Solve on lowest level */

    solve(sol[0],rhs[0]);
/* Start FMG */

    for(l=0;l<lmax;l++) {

        /* Move solution up a level */
        prolongate2D(sol[l+1], sol[l], l);
    
        /* Vcycles */
        for(i=0; i < numcycles; i++) {
            vcycle(sol,rhs[l+1],l+1,numpre,numpost);
        }
    }

    matfill(u,sol[lmax],ngrid[lmax]);

            

    

/* Free grid storage */
    for(i=0; i<=lmax;i++) {
        free(sol[i]); 
        free(rhs[i]); 
        free(res[i]); 
    }

    free(sol);
    free(rhs);
    free(res);
    free(ngrid);

    return;

}


