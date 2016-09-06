

void gauss_seidel(double *u, const double *f, const int level, const int numiter) {
    /* Red-black Gauss-Seidal relaxation, for poisson operator.
     * u contains initial guess
     */
    
    int i,j,k,indx;
    int n = 2 << level ;
    double h = 1./n;
    double h2 = h*h;
    n += 1;

    for(k=0;k<numiter;k++) {

    /* Red Points
     * Exclude boundary values 
     */
        for(i=1; i<n-1;i++) {
            for(j=1 + (i&1); j< n-1;j+=2) {    // i&1 = 0 if i is even
                indx = j + n*i;
                u[indx] = .25*(u[indx+1] + u[indx-1]
                            +  u[indx + n] + u[indx-n]
                            + h2*f[indx]);
            }
        }
    /* Black Points */
        for(i=1;i<n-1;i++) {
            for(j=2-(i&1);j<n-1;j+=2) {
                indx= j + n*i;
                u[indx] = .25*(u[indx+1] + u[indx-1]
                            +  u[indx + n] + u[indx-n]
                            + h2*f[indx]);
            }
        }

    }

    return;
}
