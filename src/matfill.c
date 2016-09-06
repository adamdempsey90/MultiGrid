#include "stdlib.h"
void matfill(double *u, const double *fill, const int n) {
    /* Fill u with the values in fill array.
     * If fill is NULL, then fill u with zeros
     */
    int i;

    if (fill == NULL) {
        for(i=0;i<n;i++) u[i] = 0.0;
    }
    else {
        for(i=0;i<n;i++) u[i] = fill[i];
    }

    return;
}
