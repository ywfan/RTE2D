struct matlab_variables defaultSetup( )
{
    struct matlab_variables vars;
    vars.R = 20;
    vars.Bsor = 1; 
    vars.Asor = 1;
    vars.Nsor = 16; 
    vars.Bdet = 1; 
    vars.Adet = 1; 
    vars.Ndet = 16;
    vars.index_i = 1.5;
    vars.index_o = 1.0;
    vars.alevel = 5;
    vars.alevel0 = 0;
    vars.slevel = 4;
    vars.slevel0 = 0;
    vars.tol = 1.e-6;
    vars.whichmg = 5;
    vars.FMG = 1;
    vars.n1 = 3;
    vars.n2 = 3;
    vars.n3 = 1;
    vars.n_max = 10;

    vars.n_ua = 0;
    vars.n_us = 2;
    vars.val_ua = malloc((vars.n_ua+1)*sizeof(double));
    vars.val_ua[0] = 0.01;
    vars.val_us = malloc((vars.n_us+1)*sizeof(double));
    vars.val_us[0] = 1.;
    vars.val_us[1] = 2.;
    vars.val_us[2] = 2.;

    FILE* f = fopen("shape.txt", "r");
    fscanf(f, "%d", &(vars.n_shape));
    vars.shape = malloc((vars.n_shape)*sizeof(double *));
    int i, j;
    for (j = 0; j < vars.n_shape; ++j) 
    {
       vars.shape[j] = malloc(5*sizeof(double));
       for(i = 0; i < 5; ++i)
       {
           fscanf(f, "%le", &(vars.shape[j][i]));
       }
    }
    fclose(f);
    return vars;
}
