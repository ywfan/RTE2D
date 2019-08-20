struct angularmesh{         // angular mesh
	int ns;                 // number of angular nodes (data from MATLAB)
	double **a;             // angular coordinates: a[ns][3] with a[i][3]:=(cos(theta(i)), sin(theta(i)), theta(i)). (data from MATLAB)
	double **w;             // angular weights: w[ns][ns] (data from MATLAB)
};

struct spatialmesh{         // spatial mesh
    int nt;                 // number of triangles (data from MATLAB)
    int np;                 // number of nodes (data from MATLAB)
    int ne;                 // number of edges (data from MATLAB)
    int **so;               // sweep ordering: so[ns][nt] (data from MATLAB)
    double **p;             // nodal coordinates: p[np][2] (data from MATLAB)
    int **t;                // nodes contained in one triangle: t[nt][3] (data from MATLAB)
    int **e;                // nodes contained in one boundary edge: e[ne][4] (data from MATLAB and see "boundary")
    double **c;             // center of triangles: c[nt][2] (see "initialization")
    double **ec;            // center of edges: ec[ne][2] (see "initialization")
    double *a;              // area of triangles: a[nt] (see "initialization")
    int **p2;               // triangles adjacent to one node: p2[np][p2[np]+1] (see "initialization")
    int **smap;             // spatial position mapping between coarse and fine mesh: smap[nt_c][smap[nt_c]+1] (see "spatialmapping")
    double ****cf;          // spatial coarse-to-fine nodal-value mapping: cf[nt_c][3][smap[nt_c]][3] (see "spatialmapping2")
    double ****fc;          // spatial fine-to-coarse nodal-value mapping: fc[nt_c][3][smap[nt_c]][3] (see "spatialmapping2")
    int *ori;               // boundary information: ori[ne] (see "boundary")
    int **e2;               // boundary information: e2[ne][2] (see "boundary")
    int **so2;              // boundary information: so2[nt][3] (see "boundary")
    double **n;             // boundary information: n[ne][2] (see "boundary")
    int ***bd;              // edge upwind or outgoing flux: bd[ns][nt][3*3] (see "edgeterm")
    double ***bd2;          // edge upwind or outgoing flux: bd[ns][nt][3]   (see "edgeterm")
};

struct boundarycoupling{    // boundary reflection or refraction coupling (see "bc_reflection")
    double **ratio;
    int **nr;
    int ***r;              // reflection contribution to incoming flux from outgoing flux: ri[ne][ns][3]    --- interpolated direction
    double ***r2;          // reflection contribution to incoming flux from outgoing flux: ri2[ne][ns][3]   --- interpolated weight
};

struct variables_fp{
    int level;
    int alevel;
    int alevel0;
    int slevel;
    int slevel0;
    int whichmg;
    int FMG;
    int vacuum;
    int n1;
    int n2;
    int n3;
    int n_max;
    double index_i;
    double index_o;
    double tol;
    double **ua;
    double **us;
    int **noflevel;
    double ****flux;
    double ****RHS;
    double ****d;
    struct boundarycoupling *b;
	struct angularmesh *amesh;
	struct spatialmesh *smesh;
};

//////////////// fyw ///////////////////
struct matlab_variables {
    double R;
    int Bsor; 
    int Asor;
    int Nsor; 
    int Bdet; 
    int Adet; 
    int Ndet;
    double index_i;
    double index_o;
    int alevel;
    int alevel0;
    int slevel;
    int slevel0;
    double tol;
    int whichmg;
    int FMG;
    int n1;
    int n2;
    int n3;
    int n_max;

    int n_ua;
    int n_us;
    double* val_ua;
    double* val_us;

    int n_shape;
    double** shape;
};
//////////////////////////////////////

struct detector{
    int B;
    int A;
    int n;
    int *t;
    double ***i;
};

// 1. initialization.h
// struct variables_fp initialization_fp(const mxArray *para_amesh,const mxArray *para_smesh,const mxArray *para_fp,const mxArray *para_uaus);
// struct variables_fp initialization_fp( );
struct variables_fp initialization_fp(struct matlab_variables vars);
// 2. spatialmapping.h
void spatialmapping(struct spatialmesh cmesh,struct spatialmesh fmesh,int **smap);
void spatialmapping2(struct spatialmesh cmesh,struct spatialmesh fmesh,int **smap,double ****cf,double ****fc,int tempsize);
// 3. boundary.h
void boundary(int ne,int nt,int **t,int **p2,double **p,int **e,int **e2,int **so2,double **n,int *ori);
// 4. edgeterm.h
void edgeterm(int M,double theta[3],double **p,int **p2,int **t,int **bd,double **bd2,int **so2);
// 5. bc_reflection.h
void bc_reflection(int ns,double **theta,struct spatialmesh smesh,double index_i,double index_o,struct boundarycoupling b);
void intepolation_a(double theta_m,double dtheta,int ns,int *b,double *b2,double constant);
double reflection(double theta_i,double ni,double no);

// 7. mgcycle.h
double mgcycle(struct angularmesh *amesh,struct spatialmesh *smesh,
struct boundarycoupling *b,double ****RHS,double **ua,double **us,double ****flux,double ****d,
int n1,int n2,int alevel,int alevel0,int slevel,int slevel0,int Ns,int vacuum,int whichmg);
double residual(int nt,int ns,double ***d,double *a);
// 8. relaxation.h
void relaxation(int Ns,struct angularmesh amesh,struct spatialmesh smesh,double ***RHS,double *ua,double *us,double ***flux,
struct boundarycoupling bb,int vacuum);
// 9. matrixsolver3by3.h
void matrixsolver(double A[3][3],double B[3],double sol[3]);
// 10. defect.h
void defect(int Ns,struct angularmesh amesh,struct spatialmesh smesh,double ***RHS,double *ua,double *us,double ***flux,
struct boundarycoupling bb,double ***res,int vacuum);
// 11. interpolation.h
void ftoc_a(int nt,int ns_c,double ***flux,double ***cflux);
void ctof_a(int nt,int ns_c,double ***flux,double ***dcflux);
void ftoc_s(int nt_c,int ns,double ***flux,double ***cflux,int **smap,double ****fc);
void ctof_s(int nt_c,int ns,double ***flux,double ***dcflux,int **smap,double ****cf);
void ftoc(int nt_c,int ns_c,double ***flux,double ***cflux,int **smap,double ****fc);
void ctof(int nt_c,int ns_c,double ***flux,double ***dcflux,int **smap,double ****cf);
void ftoc_s3(int nt_c,double *f,double *cf,int **smap,double *area);
void ctof_s3(int nt_c,double *f,double *cf,int **smap);

// 13. free_space_fp.h
void free_space_fp(struct variables_fp fp);

//struct detector det_fp(struct variables_fp fp,int SOR,const mxArray *para_sd);
struct detector det_fp(struct variables_fp fp, struct matlab_variables vars, int SOR);

#include "edgeterm_fp.h"
#include "spatialmapping_fp.h"
#include "boundary_fp.h"
#include "bc_reflection2_fp.h"
#include "relaxation_fp.h"
#include "defect_fp.h"
#include "matrixsolver3by3_fp.h"
#include "interpolation_fp.h"
#include "initialization_fp.h"
#include "mgcycle_fp.h"
#include "free_space_fp.h"
#include "detector_fp.h"
#include "defaultSetup.h"

void RTE_2D(struct variables_fp fp)
{
    int alevel,slevel,level,alevel0,slevel0,whichmg,vacuum,n1,n2,FMG,n3,n_max;
    double tol;
	struct boundarycoupling *b;
    struct angularmesh *amesh;
    struct spatialmesh *smesh;
	double ****flux,****RHS,****d;
	double **ua,**us;
	int **noflevel;

	int i,j,k,m,n,NS,nt1,nt2,ns1,ns2,nf=0;
	double res=0,res0=1,rho=1.0;
	clock_t t1,t2;

    alevel0=fp.alevel0;slevel0=fp.slevel0;
    whichmg=fp.whichmg;vacuum=fp.vacuum;n1=fp.n1;n2=fp.n2;FMG=fp.FMG;n3=fp.n3;n_max=fp.n_max;
    tol=fp.tol;ua=fp.ua;us=fp.us;amesh=fp.amesh;smesh=fp.smesh;

    alevel=fp.alevel;slevel=fp.slevel;level=fp.level;
    b=fp.b;flux=fp.flux;RHS=fp.RHS;d=fp.d;noflevel=fp.noflevel;

	NS=amesh[fp.alevel].ns;

    // Step 0: initial residual
    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
    res=residual(smesh[slevel].nt,amesh[alevel].ns,d[level],smesh[slevel].a);
    tol=tol*res;
//    printf(" Initial residual is %f and the tolerance is %f\n",res,tol);

    // Step 1: good initial guess by full multigrid method
	for(n=0;n<=level;n++)
	{   for (i=0;i<amesh[noflevel[n][1]].ns;i++)
        {	for(j=0;j<smesh[noflevel[n][0]].nt;j++)
            {   for(k=0;k<3;k++)
                {flux[n][i][j][k]=0;}
            }
        }
	}

    if (FMG==1)
    {   nt2=smesh[noflevel[level][0]].nt;ns2=amesh[noflevel[level][1]].ns;
        for(n=level-1;n>=0;n--)
        {   nt1=smesh[noflevel[n][0]].nt;ns1=amesh[noflevel[n][1]].ns;
            if (nt1==nt2)
            {ftoc_a(nt1,ns1,RHS[n+1],RHS[n]);}
            else
            {   if(ns1==ns2)
                {ftoc_s(nt1,ns1,RHS[n+1],RHS[n],smesh[noflevel[n][0]+1].smap,smesh[noflevel[n][0]+1].fc);}
                else
                {ftoc(nt1,ns1,RHS[n+1],RHS[n],smesh[noflevel[n][0]+1].smap,smesh[noflevel[n][0]+1].fc);}
            }
            nt2=nt1;ns2=ns1;
        }

        nt1=smesh[noflevel[0][0]].nt;ns1=amesh[noflevel[0][1]].ns;

        for(n=0;n<level;n++)
        {
            if(whichmg==MG4_a)
            {   if (mod(level-n,2)==0)
                {   for(i=0;i<n3;i++)
                    {mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,noflevel[n][1],alevel0,noflevel[n][0],slevel0,NS,vacuum,MG4_a);}
                }
                else
                {   for(i=0;i<n3;i++)
                    {mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,noflevel[n][1],alevel0,noflevel[n][0],slevel0,NS,vacuum,MG4_s);}
                }
            }
            else
            {   if(whichmg==MG4_s)
                {   if (mod(level-n,2)==0)
                    {   for(i=0;i<n3;i++)
                        {mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,noflevel[n][1],alevel0,noflevel[n][0],slevel0,NS,vacuum,MG4_s);}
                    }
                    else
                    {   for(i=0;i<n3;i++)
                        {mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,noflevel[n][1],alevel0,noflevel[n][0],slevel0,NS,vacuum,MG4_a);}
                    }
                }
                else
                {   for(i=0;i<n3;i++)
                    {mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,noflevel[n][1],alevel0,noflevel[n][0],slevel0,NS,vacuum,whichmg);}
                }
            }

            nt2=smesh[noflevel[n+1][0]].nt;ns2=amesh[noflevel[n+1][1]].ns;
            if (nt1==nt2)
            {ctof_a(nt1,ns1,flux[n+1],flux[n]);}
            else
            {   if(ns1==ns2)
                {ctof_s(nt1,ns1,flux[n+1],flux[n],smesh[noflevel[n][0]+1].smap,smesh[noflevel[n][0]+1].cf);}
                else
                {ctof(nt1,ns1,flux[n+1],flux[n],smesh[noflevel[n][0]+1].smap,smesh[noflevel[n][0]+1].cf);}
            }
            nt1=nt2;ns1=ns2;
            for(m=0;m<=n;m++)
            {   for (i=0;i<amesh[noflevel[m][1]].ns;i++)
                {	for(j=0;j<smesh[noflevel[m][0]].nt;j++)
                    {   for(k=0;k<3;k++)
                        {flux[m][i][j][k]=0;}
                    }
                }
            }
        }
    }
	// Step 2: multigrid solver on the finest mesh
	n=0;
//    printf(" Iter\t Res\t Rho \t T (in seconds)\n");
	while(n<n_max)
	{
		t1=clock();
		n++;
        res=mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel,slevel0,NS,vacuum,whichmg);
        for(m=0;m<level;m++)
        {   for(i=0;i<amesh[noflevel[m][1]].ns;i++)
            {   for(j=0;j<smesh[noflevel[m][0]].nt;j++)
                {   for(k=0;k<3;k++)
                    {flux[m][i][j][k]=0;}
                }
            }
        }

		t2=clock();
		// printf("runtime of iteration %d\t %f. res = %f\n", n, (t2-t1) / (double)CLOCKS_PER_SEC, res);
		if (n>0)
		{
		    rho*=res/res0;
//		    mexEvalString("pause(.001);");
//			printf(" %d\t%le\t%f\t",n,res,res/res0);
            if(res<tol)
            {
                rho=pow(rho,1.0/(n-1));
                nf=n;
                //n=n_max;
                break;
            }
		}
		else
		{
//		    mexEvalString("pause(.001);");
//		    printf(" %d\t%le\t%f\t",n,res,res/res0);
        }
		res0=res;
	}

    // Step 3: postsolver
	for(n=0;n<=level;n++)
	{   for (i=0;i<amesh[noflevel[n][1]].ns;i++)
        {	for(j=0;j<smesh[noflevel[n][0]].nt;j++)
            {   for(k=0;k<3;k++)
                {   RHS[n][i][j][k]=0;
                    d[n][i][j][k]=0;
                }
            }
        }
	}
}
