#include "stdio.h"
#include "stdlib.h"


extern void multigrid(double *, double *, int, int, int, int, int);
extern void init(double *, double *, int);
extern void output(double *, double *, int, char *);

int main(int argc, char *argv[]) {

    if (argc < 4) {
        printf("Need at least 4 arguments!\n");
        printf("./fmg lmin lmax numcycles outputdir\n");

    }
    int lmin, lmax,numcycles;
    char fname[256];
    double *u, *f;

    lmin = atoi(argv[1]);
    lmax = atoi(argv[2]);
    numcycles = atoi(argv[3]);
    sprintf(fname, "%soutput.dat",argv[4]);

    int n = 1 + (2 << lmax);

    u = (double *)malloc(sizeof(double)*n*n);
    f = (double *)malloc(sizeof(double)*n*n);

    init(u,f,n);
    multigrid(u,f,lmin,lmax,numcycles,1,1); 
    output(u,f,n,fname);
    free(u);
    free(f);
    return 1;
}
