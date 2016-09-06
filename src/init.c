
void init(double *u, double *rhs, const int n) {
    int i,j,indx;


    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            indx = j + i*n;
            u[indx] = 0.0;
            rhs[indx] = 1.0;
        }
    }

    return;
}
