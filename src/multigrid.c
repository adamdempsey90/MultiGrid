#include "stdlib.h"

extern void matfill(double *,  double *,  int );
extern void restrict2D(double *,  double *,  int );
extern void solve(double *, double *, int);
extern void gauss_seidel(double *, double *, int, int);
extern void prolongate2D(double *, double *, int );
extern void residual(double *, double *, double *, int);

void multigrid(double *u, double *f, const int lmin, const int lmax, const int numcycles, const int numpre, const int numpost) {
/* Multigrid V-cycle for the Poisson equation 
 * u and f are given on the lmax level and are
 * coarsened to the lmin level
 */

    int i,j,k,l;
    int n;
    int ntop = i + (2 << lmax);
    int nbot = i + (2 << lmin);
    int nlevels = lmax - lmin + 1;


/* Allocate the storage for the grids */
    int *ngrid; 
    double **sol, **rhs, **res;
    ngrid = (int *)malloc(sizeof(int) * nlevels);
    sol = (double **)malloc(sizeof(double *) * nlevels);
    rhs = (double **)malloc(sizeof(double *) * nlevels);
    res = (double **)malloc(sizeof(double *) * nlevels);
    for(i=lmin; i<=lmax;i++) {
        n = 1 + (2 << i);
        ngrid[i-lmin] = n;
        sol[i-lmin] = (double *)malloc(sizeof(double)*n*n);
        rhs[i-lmin] = (double *)malloc(sizeof(double)*n*n);
        res[i-lmin] = (double *)malloc(sizeof(double)*n*n);
    }

/* Copy initial values to finest grid */

    matfill(sol[lmax-lmin], u, ngrid[lmax-lmin]);
    matfill(rhs[lmax-lmin], f, ngrid[lmax-lmin]);
    matfill(res[lmax-lmin], NULL, ngrid[lmax-lmin]);

/* Restrict to coarser meshes */
    for(i=lmax-1; i>=lmin; i--) {
        restrict2D(sol[i-lmin], sol[i+1-lmin], i);
        restrict2D(sol[i-lmin], sol[i+1-lmin], i);
        matfill(res[i-lmin], NULL, ngrid[i]);
    }

/* Solve on lowest level */

    solve(sol[0],rhs[0], 0);

/* Start FMG */
    for(l=lmin;l<lmax;l++) {

        /* Move solution up a level */
        prolongate2D(sol[l+1-lmin], sol[l-lmin], l);
    
        /* Vcycles */
        for(i=0; i < numcycles; i++) {
            
            for(j=l; j > lmin; j--) {
                /* Going down. 
                 * Pre-smoothing
                 */
                gauss_seidel(sol[j-lmin],rhs[j-lmin], j,numpre);

                residual(res[j-lmin], sol[j-lmin], rhs[j-lmin], j);

                restrict2D(rhs[j-lmin-1],res[j-lmin],j-1);

                matfill(sol[j-lmin-1], NULL, 1 + (2 << (j-1)));
            }
            /* Solve at lowest level */
            solve(sol[0],rhs[0], 0);

            /* Going up */
            for(j=1; j<l; j++) {
                prolongate2D(sol[j-lmin],sol[j-1-lmin],j);
                for(k=0;k<1 + (2 << j);k++) {
                    sol[j-lmin][k] += res[j-lmin][k];
                }
                gauss_seidel(sol[j-lmin],rhs[j-lmin], j,numpost);
            }
        }
    }


            

    

/* Free grid storage */
    for(i=lmin; i<=lmax;i++) {
        free(sol[i-lmin]); 
        free(rhs[i-lmin]); 
        free(res[i-lmin]); 
    }
    free(sol);
    free(rhs);
    free(res);
    free(ngrid);

    return;

}
