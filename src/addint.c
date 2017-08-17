
extern void prolongate2D(double *,const double *,const int);
void addint(double *uf, const double *uc, double *res, const int level) {
    /* Prolongates solution from course mesh to fine mesh
     * and adds that to fine mesh solution
     * level is the course level, where there are 2^level + 1 points
     * res is used for temp storage
     */
    int i;
    int n = 1 + (2 << (level+1));
    prolongate2D(res,uc,level);
    for(i=0;i<n*n;i++) {
        uf[i] += res[i];
    }

    return;
}
