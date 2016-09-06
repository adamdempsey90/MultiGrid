#include "stdio.h"
void output(const double *u, const double *rhs, const int n, const char *fname) {

    FILE *f = fopen(fname, "w");

    fwrite(u, sizeof(double),n, f);
    fwrite(rhs, sizeof(double),n, f);

    fclose(f);
    return;
}
