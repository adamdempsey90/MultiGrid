#include "stdio.h"
#include "stdlib.h"
#define TRUE 1
#define FALSE 0

#define PRINT if (display) printf

extern void restrict2D(double *, double *, int );
extern void prolongate2D(double *, double *, int );

int main(int argc, char *argv[]) {
    int i,j,k,n,indx,ngrids,lf;
    double **uc, **uc_ans;
    int *nvals;
    
    lf = atoi(argv[1]);

    int display = (lf < 5);
    ngrids = lf+1;


    uc = (double **)malloc(sizeof(double *)*ngrids);
    uc_ans = (double **)malloc(sizeof(double *)*ngrids);
    nvals = (int *)malloc(sizeof(int)*ngrids);


    for(k=0;k<ngrids;k++) {
        n = 1 + (2<<k);
        nvals[k] = n; 
        uc[k] = (double *)malloc(sizeof(double)*n*n);
        uc_ans[k] = (double *)malloc(sizeof(double)*n*n); 
        for(j=0;j<n*n;j++) uc[k][j] = 0.;
    }

    
    for(k=0;k<ngrids;k++) {
        n = nvals[k];
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                if ( (i==0) || (i==n-1) || (j==0) || (j==n-1)) {
                    uc_ans[k][j+i*n] = 0.;
                }
                else {
                    uc_ans[k][j+i*n] = 1.;
                }
            }
        }
    }
                

    k = ngrids-1;
    n = nvals[k];
    PRINT("Start with\n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            indx = j + i*n;
            uc[k][indx] = uc_ans[k][indx];
            PRINT("%lg\t",uc_ans[k][indx]);
        }
        PRINT("\n");
    }
    

    k = ngrids-1;
    int pass = TRUE;
    
    while (k > 0) {

        PRINT("Going down to %d\n",k-1);
        restrict2D(uc[k-1], uc[k], k);
        k--;
        n = nvals[k];
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                indx = j + i*n;
                PRINT("%lg\t",uc[k][indx]);
            }
            PRINT("\n");
        }
        /*
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                indx = j + i*n;
                if (uc[k][indx] != uc_ans[k][indx]) {
                    pass = FALSE;
                    printf("test_restriction failed!\n");
                    printf("i\tj\tvalue\tcorrect\n");
                    printf("%d\t%d\t%lg\t%lg\n", i,j,uc[k][indx],uc_ans[k][indx]);
                    return 0;
                }
            }
        }
        */
    }

    printf("Lowest 2 Grids\n");
    for(k=1;k>=0;k--) {
        n = nvals[k];
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                indx = j + i*n;
                printf("%lg\t",uc[k][indx]);
            }
            printf("\n");
        }
        printf("\n");
    }



    printf("test_restriction Passed!\n");


    k = 0;
    while (k < ngrids-1) {

        PRINT("Going up to %d\n",k+1);
        prolongate2D(uc[k+1], uc[k], k);
        k++;
        n = nvals[k];
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                indx = j + i*n;
                PRINT("%lg\t",uc[k][indx]);
            }
            PRINT("\n");
        }
        /*
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                indx = j + i*n;
                if (uc[k][indx] != uc_ans[k][indx]) {
                    pass = FALSE;
                    printf("test_restriction failed!\n");
                    printf("i\tj\tvalue\tcorrect\n");
                    printf("%d\t%d\t%lg\t%lg\n", i,j,uc[k][indx],uc_ans[k][indx]);
                    return 0;
                }
            }
        }
        */
    }
    printf("test prolongation passed!\n");
    for(i=0;i<ngrids;i++) {
        free(uc[i]);
        free(uc_ans[i]);
    }
    free(uc);
    free(uc_ans);
    free(nvals);
    return 1;
}

