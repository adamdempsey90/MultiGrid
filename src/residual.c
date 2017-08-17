
void residual(double *res, const double *u, const double *f, const int level) {
    /* Computes the negative residual for the poisson problem. */ 
    int i,j,indx,n;
    double h,invh2;

    n = (2 << level);

    h = 1./n;
    invh2 = 1./(h*h);
    n+=1;

    for(i=1; i<n-1 ; i++) {
        for(j=1; j <n-1;j++) {
            indx=  j+n*i;

            res[indx] = -invh2*( u[indx+1] + u[indx-1] 
                                +u[indx+n] + u[indx-n]
                                - 4. * u[indx]) + f[indx];
        }
    }

    /* Boundary points have zero residual */

    for(i=0;i<n;i++) {
        res[n*i] = 0;
        res[n*i + n-1] = 0;
        res[i] = 0;
        res[i + n*(n-1)] = 0;
    }

    return;
}
