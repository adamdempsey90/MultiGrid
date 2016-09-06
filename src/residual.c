
void residual(double *res, const double *u, const double *f, const int level) {
    /* Computes the negative residual for the poisson problem. */ 
    int i,j,indx,n;
    double h,invh2;

    n = 1 + (2 << level);

    h = 1./(n-1);
    invh2 = 1./(h*h);

    for(i=1; i<n-1 ; i++) {
        for(j=1; j <n-1;j ++) {
            indx=  j+indx*i;
            res[indx] = -invh2*( u[indx+1] + u[indx-1] 
                                +u[indx+n] + u[indx-n]
                                - 4 * u[indx]) + f[indx];
        }
    }

    /* Boundary points have zero residual */

    for(i=0;i<n;i++) {
        indx = n*i;
        res[indx] = 0;
        res[indx + n-1] = 0;
    }
    for(j=1;j < n-1;j++) {
        indx = n*(n-1);
        res[j] = 0;
        res[j + indx] = 0;
    }

    return;
}
