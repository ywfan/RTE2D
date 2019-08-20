#include "utilityfunctions.h"
#include "RTE_2D.h"

void RTE_2D(struct variables_fp fp);

double** DOT2D(struct matlab_variables vars)
{
    struct variables_fp fp;
    struct detector sor,det;
    int j,k,tri,nt,ns,is,id,alevel,slevel,level;
    double source_corr,temp,*a,s_int2[3],**measurement; //,*output;

    // 0. initialization
    // load, malloc and compute "amesh", "smesh" and "fp" ("noflevel", "ua", "us", "RHS", "d", "flux" and "b", see "initlalization")
    fp = initialization_fp(vars);

    alevel=fp.alevel;slevel=fp.slevel;level=fp.level;
    a=fp.smesh[slevel].a;

    sor=det_fp(fp, vars, 1); // ,prhs[4]);
    det=det_fp(fp, vars, 0); // ,prhs[4]);

    ns=fp.amesh[alevel].ns;nt=fp.smesh[slevel].nt;

    measurement=malloc(sor.n*sizeof(double *));
    for(is=0;is<sor.n;is++)
    {measurement[is]=malloc(det.n*sizeof(double));}


    for(is=0;is<sor.n;is++)
    {
        tri=sor.t[is];
        // 1. assemble light source into "RHS"
        source_corr=a[tri]/12;
        for(j=0;j<ns;j++)
        {   
            temp=0;
            for(k=0;k<3;k++)
            {   
                s_int2[k]=sor.i[is][j][k]/source_corr;
                temp+=s_int2[k];
            }
            temp=temp/4;
            for(k=0;k<3;k++)
            {
                fp.RHS[level][j][tri][k]=s_int2[k]-temp;
            }
        }

        // 2. solving RTE
        RTE_2D(fp);

        // measurement
        for(id=0;id<det.n;id++)
        {
            measurement[is][id]=0;
            tri=det.t[id];
            for(j=0;j<ns;j++)
            {   for(k=0;k<3;k++)
                {
                    measurement[is][id]+=det.i[id][j][k]*fp.flux[level][j][tri][k];
                }
            }
        }
    }
    // 4. free variables
    // for(is=0;is<sor.n;is++)
    // {free(measurement[is]);}
    // free(measurement);

    for(j=0;j<det.n;j++)
    {   for(k=0;k<ns;k++)
        {free(det.i[j][k]);}
        free(det.i[j]);
    }
    free(det.i);
    free(det.t);

    for(j=0;j<sor.n;j++)
    {   for(k=0;k<ns;k++)
        {free(sor.i[j][k]);}
        free(sor.i[j]);
    }
    free(sor.i);
    free(sor.t);

    free_space_fp(fp);

    return measurement;
}
