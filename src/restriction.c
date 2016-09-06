
void restrict2D(double *uc, const double *uf, const int level) {
    /* Restrict from fine grid to coarse grid.
     * Uses half weighting.
     * Fine grid has 2^level + 1 points.
     * Coarse grid has 2^(level-1) + 1 points
     */

    int nc,nf;
    int ic, jc, iif, jf;
    int indx;

    nf = 2 << level;
    nc = nf >> 1;
    nf += 1;
    nc += 1;

/* Interior points */
    for(iif = 2; iif < nf-2; iif += 2) {
        for(jf=2; jf < nf-2; jf += 2) {
            ic = iif >> 1;
            jc = jf >> 1;
            indx = jf + nf*iif;

            uc[jc + nc*ic] = .5 * uf[indx]
                            + .125*( uf[indx+1]
                                    + uf[indx -1]
                                    + uf[indx+nf]
                                    + uf[indx-nf]);
        }
    }

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
