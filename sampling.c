#include <string.h>
#include <stdlib.h>
#include "DOT2D.h"
#include "generateGeo.h"

void resetShape(struct matlab_variables* p_vars)
{
    int j;
    for(j = 0; j < p_vars->n_shape; ++j)
      free(p_vars->shape[j]);
    free(p_vars->shape);
    p_vars->shape = generateGeo(p_vars->n_shape, 0.1, 100);
    for(j = 0; j < p_vars->n_shape; ++j)
    {
        p_vars->shape[j][0] *= p_vars->R;
        p_vars->shape[j][1] *= p_vars->R;
        p_vars->shape[j][2] *= p_vars->R;
        p_vars->shape[j][3] *= p_vars->R;
    }
}

int main(int argc, char* argv[])
{
    if(argc < 4) {
        printf("Usage: %s <N_shape> <N_samples> <filename> <seed> \n", argv[0]);
        return 1;
    }
    int Nshape = atoi(argv[1]);
    int Nsample = atoi(argv[2]);
    int seed = atoi(argv[4]);
    char* filename_shape; 
    char* filename_measure;
    filename_shape = (char*)malloc((strlen(argv[3]) + 6) * sizeof(char));
    filename_measure = (char*)malloc((strlen(argv[3]) + 6) * sizeof(char));
    sprintf(filename_shape, "%s_s.txt", argv[3]);
    sprintf(filename_measure, "%s_m.txt", argv[3]);
    printf("filenames: %s\t%s\n", filename_shape, filename_measure);

    int i, j, k;
    struct matlab_variables vars;
    clock_t t1, t2;
    double** measurement;
    vars = defaultSetup( );
    vars.slevel = 5;
    srand(time(0) + 100000 * seed);

    for(j = 0; j < vars.n_shape; ++j)
      free(vars.shape[j]);
    free(vars.shape);
    vars.n_shape = Nshape;
    vars.n_us = vars.n_shape;
    free(vars.val_us);
    vars.val_us = malloc((vars.n_us+1)*sizeof(double));
    vars.val_us[0] = 1.;
    for(i = 0; i < vars.n_us; ++i)
      vars.val_us[i+1] = 2.;

    vars.shape = generateGeo(vars.n_shape, 0.1, 100);

    FILE* f_s = fopen(filename_shape, "w");
    FILE* f_m = fopen(filename_measure, "w");
    for(k = 0; k < Nsample; ++k)
    {
        resetShape(&vars);
        t1 = clock();
        measurement = DOT2D(vars);
		t2=clock();
		printf("sample %d, runtime: %f\n", k+1, (t2-t1) / (double)CLOCKS_PER_SEC);
        fprintf(f_s, "%d\t", vars.n_shape);
        for(i = 0; i < vars.n_shape; ++i) {
            for(j = 0; j < 5; ++j)
              fprintf(f_s, "%.16e\t", vars.shape[i][j]);
        }
        fprintf(f_s, "\n");
        for(i=0;i<vars.Nsor;i++)
        {   for(j=0;j<vars.Ndet;j++)
            {
                fprintf(f_m, "%.16e\t", measurement[i][j]);
            }
        }
        fprintf(f_m, "\n");
        for(i=0; i < vars.Nsor; ++i)
          free(measurement[i]);
        free(measurement);
    }
    fclose(f_s);
    fclose(f_m);

    free(vars.val_ua);
    free(vars.val_us);
    for(j = 0; j < vars.n_shape; ++j)
      free(vars.shape[j]);
    free(vars.shape);
}
