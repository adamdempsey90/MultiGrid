
void prolongate2D(double *uf, const double *uc, const int level) {
    /* Prolongate from coarse to fine mesh.
     * The coarse mesh has 2^level + 1 points.
     * The fine mesh has 2^(level+1) + 1 points.
     * Uses bilinear interpolation.
     */
    int nc, nf;
    int ic, jc, iif, jf;
    int indx;

    nc = 2 << level;
    nf = nc << 1;
    nc += 1;
    nf += 1;

/* Copy values at matching nodes */
    for(ic=0; ic < nc; ic++) {
        for(jc = 0; jc < nc; jc++) {
            uf[2*jc + 2*ic*nf] = uc[jc + nc*ic];
        }
    }

/* Odd # columns horizontal interpolation */

    for(iif = 0; iif < nf; iif += 2) {
        for(jf=1; jf < nf; jf +=2) {
            indx = jf + iif*nf;
            uf[indx] = .5*( uf[indx - nf] + uf[indx + nf]);
        }
    }

/* Even # colums vertical interpolation */

    for(iif=1; iif < nf; iif +=2) {
        for(jf=0;jf<nf;jf+=2) {
            indx = jf +iif*nf;
            uf[indx] = .5*(uf[indx-1] + uf[indx+1]);
        }
    }

/* In between points */
    for(iif=1; iif < nf; iif += 2) {
        for(jf=1; jf<nf; jf+=2) {
            indx = jf + iif*nf;
            uf[indx] = .5*( uf[indx-1] 
                            + uf[indx+1] 
                            + uf[indx + nf]
                            + uf[indx - nf])
                       + .25*(uf[indx -1 - nf]
                               + uf[indx +1 + nf]
                               + uf[indx +1 - nf]
                               + uf[indx -1 + nf]);
        }
    }

    return;
}
