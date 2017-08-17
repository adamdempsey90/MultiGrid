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

    int ii,jj;
    int i,j,k,l;
    int n,n_curr;
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

    printf("Allocated\n");

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

    printf("Coarse\n");
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++) {
            printf("%lg\t",rhs[0][j+3*i]);
        }
        printf("\n");
    }

    printf("2\n");
/* Solve on lowest level */

    solve(sol[0],rhs[0]);
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++) {
            printf("%lg\t",rhs[0][j+i*3]);
        }
        printf("\n");
    }

    printf("3\n");
/* Start FMG */
 //   vcycle(sol,rhs,res,lmax,numpre,numpost);
    for(l=0;l<lmax;l++) {

        printf("Starting Vcycle %d\n",l);
        /* Move solution up a level */
        prolongate2D(sol[l+1], sol[l], l);
        printf("3.1\n");
    
        /* Vcycles */
        for(i=0; i < numcycles; i++) {
            vcycle(sol,rhs[l+1],l+1,numpre,numpost);
        }
        printf("Finished Vcycle %d\n",l);
    }
    printf("Ended FMG\n");

    matfill(u,sol[lmax],ngrid[lmax]);

            

    

/* Free grid storage */
    for(i=0; i<=lmax;i++) {
        printf("free sol %d\n",i);
        free(sol[i]); 
        printf("free rhs %d\n",i);
        free(rhs[i]); 
        printf("free res %d\n",i);
        free(res[i]); 
    }

    printf("5\n");
    printf("free sol\n");
    free(sol);
    printf("free rhs\n");
    free(rhs);
    printf("free res\n");
    free(res);
    printf("free ngrid\n");
    free(ngrid);

    return;

}


