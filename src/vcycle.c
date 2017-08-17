#include <stdlib.h>
#include <stdio.h>
extern void gauss_seidel(double *,double *,int,int);
extern void residual(double *, double *,double *, int);
extern void restrict2D(double *,double *,int);
extern void addint(double *, double *, double *,int);
extern void solve(double *, double *);
extern void matfill(double *, double *, int);

void vcycle(double **u, double *f , int j, const int numpre, const int numpost) {
    /* Recursive function to do one V-cycle loop 
     * u contains the initial solution on the finest grid
     * f contains the rhs on the finest grid
     * there are 2^j + 1 points on the finest grid
     */
    printf("Vcycle %d\n",j);
    if (j == 0) {
    /* Exact solve at bottom of V */
        solve(u[j],f);
        return;
    }

    int nf = 1 + (2 << j);
    int nc = 1 + (2 << (j-1));

    double *tmp = (double *)malloc(sizeof(double)*nf*nf);
    double *res = (double *)malloc(sizeof(double)*nc*nc);

    matfill(tmp,NULL,nf*nf);
    matfill(res,NULL,nc*nc);


    /* 1. Downstroke 
     * 1.1 Pre-smoothing 
     */
    gauss_seidel(u[j],f, j,numpre);
    /* 1.2 Compute residual on this level */
    residual(tmp, u[j], f, j);
    /* 1.3 Restrict that residual down another level*/
    restrict2D(res,tmp,j-1);
    /* 1.4 Repeat until we hit the bottom */
    matfill(u[j-1],NULL,nc*nc);
    vcycle(u,res, j-1,numpre,numpost);

    /* 2. Upstroke
     * 2.1 Prolongate solution up a level
     */

    printf("V %d\tAddint %d -> %d\n",j,j-1,j);
    addint(u[j],u[j-1],tmp,j-1);

    /* Apply any post smoothing */
    gauss_seidel(u[j],f, j,numpost);

    free(tmp);
    free(res);


    return;
}
