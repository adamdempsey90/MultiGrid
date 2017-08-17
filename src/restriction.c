
void restrict2D(double *uc, const double *uf, const int level) {
    /* Restrict from fine grid to coarse grid.
     * Uses half weighting.
     * Fine grid has 2^level + 1 points.
     * Coarse grid has 2^(level-1) + 1 points
     */

    int nc,nf;
    int ic, jc, iif, jf;
    int i,j;
    int indx;

    nf = 1 + (2 << (level+1));
    nc = 1 + (2 << level);


/* Interior points */
    for(i=1;i<nc-1;i++) {
        for(j=1;j<nc-1;j++) {
            indx = 2*j + 2*i*nf;
#ifdef FULL
            uc[j + nc*i] = .25 * uf[indx]
                            + .125*( uf[indx+1]
                                    + uf[indx -1]
                                    + uf[indx+nf]
                                    + uf[indx-nf])
                            + .0625*( uf[indx+1+nf]
                                    + uf[indx+1-nf]
                                    + uf[indx-1+nf]
                                    + uf[indx-1-nf]);
#else
            uc[j + nc*i] = .5 * uf[indx]
                            + .125*( uf[indx+1]
                                    + uf[indx -1]
                                    + uf[indx+nf]
                                    + uf[indx-nf]);
#endif
        }
    }

/*
    for(iif = 2; iif < nf-2; iif += 2) {
        for(jf=2; jf < nf-2; jf += 2) {
            ic = iif >> 1;
            jc = jf >> 1;
            indx = jf + nf*iif;

          //  uc[jc + nc*ic] = .5 * uf[indx]
          //                  + .125*( uf[indx+1]
          //                          + uf[indx -1]
          //                          + uf[indx+nf]
          //                          + uf[indx-nf]);
            uc[jc + nc*ic] = .25 * uf[indx]
                            + .125*( uf[indx+1]
                                    + uf[indx -1]
                                    + uf[indx+nf]
                                    + uf[indx-nf])
                            + .0625*( uf[indx+1+nf]
                                    + uf[indx+1-nf]
                                    + uf[indx-1+nf]
                                    + uf[indx-1-nf]);
        }
    }
*/
/* Boundary Points */

    for(iif = 0; iif < nf; iif += nf-1) {
        for(jf=0; jf < nf; jf+= 2) {
            ic = iif >> 1;
            jc = jf >> 1;
            uc[jc + nc*ic] = uf[jf + iif*nf];
        }
    }

    for(iif=1; iif < nf-1; iif += 2) {
        for(jf=0; jf < nf; jf += nf-1) {
            ic = iif >> 1;
            jc = jf >> 1;
            uc[jc + nc*ic] = uf[jf + nf*iif];
        }
    }

    return;
}
