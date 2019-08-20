
//struct variables_fp initialization_fp(const mxArray *para_amesh,const mxArray *para_smesh,const mxArray *para_fp,const mxArray *para_uaus)
struct variables_fp initialization_fp(struct matlab_variables vars)
{
   	int tempsize=20,tempsize2=20;
   	int **t,**e,**p2,**smap;
   	double **p;
   	int i,j,k,m,n,nt,np,ne,ns,da,ds;
   	double x1,x2,x3,y1,y2,y3;

    int alevel,alevel0,slevel,slevel0,level=-1,whichmg,vacuum,whichslevel;
    double index_i,index_o;
    struct angularmesh *amesh;
    struct spatialmesh *smesh;
	struct boundarycoupling *b=NULL;
	double ****flux,****RHS,****d;
	double **ua,**us; //,*utemp;
	int **noflevel;
    struct variables_fp fp;
    double tempd;

    // 1.1. load parameters from "para_fp"
    // load para_fp from "parafp.txt"
    index_i = vars.index_i;
    index_o = vars.index_o;
    alevel = vars.alevel - 1;
    alevel0 = vars.alevel0;
    slevel = vars.slevel - 1;
    slevel0 = vars.slevel0;
    fp.tol = vars.tol;
    whichmg = vars.whichmg;
    fp.FMG = vars.FMG;
    fp.n1 = vars.n1;
    fp.n2 = vars.n2;
    fp.n3 = vars.n3;
    fp.n_max = vars.n_max;
   	FILE *f;
    // f = fopen("parafp.txt", "r");
    // fscanf(f, "%le", &index_i);
    // fscanf(f, "%le", &index_o);
    // fscanf(f, "%d", &alevel);
    // alevel -= 1;
    // alevel0 = 0;
    // fscanf(f, "%d", &slevel);
    // slevel -= 1;
    // slevel0 = 0;
    // fscanf(f, "%le", &(fp.tol));
    // fscanf(f, "%d", &whichmg);
    // fscanf(f, "%d", &(fp.FMG));
    // fscanf(f, "%d", &(fp.n1));
    // fscanf(f, "%d", &(fp.n2));
    // fscanf(f, "%d", &(fp.n3));
    // fscanf(f, "%d", &(fp.n_max));
    // fscanf(f, "%le", &(fp.R));
    // fclose(f);

	if(abs2(index_i-index_o)/index_i<0.01) // refraction index mismatch at the boundary
	{vacuum=1;}
	else
	{vacuum=0;}

    level=computelevel(alevel,alevel0,slevel,slevel0,whichmg);// compute "level" (the indicator of mesh levels in multigrid)

	// 1.2. load para.amesh
	amesh=malloc((alevel+1)*sizeof(struct angularmesh));
    f = fopen("amesh.txt", "r");
    for (i=0;i<=alevel;i++)
    {
        fscanf(f, "%le", &tempd); amesh[i].ns=(int)tempd;
        amesh[i].a=malloc(amesh[i].ns*sizeof(double *));
        for(j=0;j<amesh[i].ns;j++)
        {   amesh[i].a[j]=malloc(3*sizeof(double));
            for(k=0;k<3;k++)
            {
                fscanf(f, "%le", &(amesh[i].a[j][k]));
            }
        }
        amesh[i].w=malloc(amesh[i].ns*sizeof(double *));
        for(j=0;j<amesh[i].ns;j++)
        {	amesh[i].w[j]=malloc(amesh[i].ns*sizeof(double));
            for(k=0;k<amesh[i].ns;k++)
            {fscanf(f, "%le", &(amesh[i].w[j][k]));}
        }
    }
    fclose(f);

	// 1.3. load smesh.txt
	//      Notice the index difference in c programming: array indexes from 0 instead of 1,
	//      we subtract "1" from every integer-valued index here as for "so", "t" and "e" as follow.
    f=fopen("smesh.txt", "r");
	smesh=malloc((slevel+1)*sizeof(struct spatialmesh));
	for (i=0;i<=slevel;i++)
	{
		fscanf(f, "%le", &tempd);smesh[i].nt=(int)tempd;
		fscanf(f, "%le", &tempd);smesh[i].np=(int)tempd;
        fscanf(f, "%le", &tempd);smesh[i].ne=(int)tempd;
        smesh[i].so=malloc(amesh[alevel].ns*sizeof(int *));
        for(j=0;j<amesh[alevel].ns;j++)
        {   smesh[i].so[j]=malloc(smesh[i].nt*sizeof(int));
            for(k=0;k<smesh[i].nt;k++)
            {fscanf(f, "%le", &tempd);smesh[i].so[j][k]=(int)tempd-1;}
        }
        smesh[i].p=malloc(smesh[i].np*sizeof(double *));
        for(j=0;j<smesh[i].np;j++)
        {   smesh[i].p[j]=malloc(2*sizeof(double));
            for(k=0;k<2;k++)
            {fscanf(f, "%le", &tempd);smesh[i].p[j][k]=tempd;}
        }
        smesh[i].t=malloc(smesh[i].nt*sizeof(int *));
        for(j=0;j<smesh[i].nt;j++)
        {   smesh[i].t[j]=malloc(3*sizeof(int));
            for(k=0;k<3;k++)
            {fscanf(f, "%le", &tempd);smesh[i].t[j][k]=(int)tempd-1;}
        }
        smesh[i].e=malloc(smesh[i].ne*sizeof(int *));
        for(j=0;j<smesh[i].ne;j++)
        {   smesh[i].e[j]=malloc(4*sizeof(int));
            smesh[i].e[j][0]=-1;smesh[i].e[j][3]=-1;
            for(k=1;k<3;k++)
            {fscanf(f, "%le", &tempd);smesh[i].e[j][k]=(int)tempd-1;}
        }
	}
	fclose(f);

	// 2.1. compute "c", "ec", "a" and "p2"
    //      p2[np][p2[np][0]+1]: triangles adjacent to one node
    //      For the 2nd index, the 1st element is the total number of triangles adjacent to this node,
    //      the corresponding triangles are saved from the 2nd.
    //      Example:    assume triangles [5 16 28 67] are adjacent to the node 10, then
    //                  p2[10][0]=4, p2[10][1]=5, p2[10][2]=16, p2[10][3]=28 and p2[10][4]=67.
	for(i=0;i<=slevel;i++)
	{   
        p=smesh[i].p;t=smesh[i].t;nt=smesh[i].nt;np=smesh[i].np;e=smesh[i].e;ne=smesh[i].ne;
        smesh[i].c=malloc(nt*sizeof(double *));
        for(j=0;j<nt;j++)
        {   smesh[i].c[j]=malloc(2*sizeof(double));
            for(k=0;k<2;k++)
            {smesh[i].c[j][k]=(p[t[j][0]][k]+p[t[j][1]][k]+p[t[j][2]][k])/3;}// center of triangle
        }

        smesh[i].ec=malloc(ne*sizeof(double *));
        for(j=0;j<ne;j++)
        {
            smesh[i].ec[j]=malloc(2*sizeof(double));
            for(k=0;k<2;k++)
            {smesh[i].ec[j][k]=(p[e[j][1]][k]+p[e[j][2]][k])/2;}// center of edge
        }

        smesh[i].a=malloc(nt*sizeof(double));
	    for(j=0;j<nt;j++)
        {
            x1=p[t[j][0]][0];y1=p[t[j][0]][1];
            x2=p[t[j][1]][0];y2=p[t[j][1]][1];
            x3=p[t[j][2]][0];y3=p[t[j][2]][1];
            smesh[i].a[j]=AREA(x1,y1,x2,y2,x3,y3);//area of triangle
        }

        p2=malloc(np*sizeof(int *));
        for(j=0;j<np;j++)
        {   p2[j]=malloc(tempsize*sizeof(int));
            // tempsize is the initial length of the second index of "p2", and it may cause problem if it is too small.
            p2[j][0]=0;
            for(k=1;k<tempsize;k++)
            {p2[j][k]=-1;}
        }

        for(j=0;j<nt;j++)
        {   for(k=0;k<3;k++)
            {
                p2[t[j][k]][0]+=1;
                p2[t[j][k]][p2[t[j][k]][0]]=j;
            }
        }

        for(j=0;j<np;j++)
        {   if(p2[j][0]+1>tempsize)
            {   printf("WARNING: tempsize for p2 is too small!!\n");
                goto stop;
            }
        }

	    smesh[i].p2=malloc(np*sizeof(int *));
        for(j=0;j<np;j++)
        {   smesh[i].p2[j]=malloc((p2[j][0]+1)*sizeof(int));
            for(k=0;k<=p2[j][0];k++)
            {smesh[i].p2[j][k]=p2[j][k];}
        }

        for(j=0;j<np;j++)
        {free(p2[j]);}
        free(p2);
    }
	// 2.2. compute "smap", "cf" and "fc"
    //      For the data structure of "smap", see "spatialmapping";
    //      For the data structure of "cf" and "fc", see "spatialmapping2".
    //      Note: those arrays are not defined on the coarsest spatial mesh (i=0),
    //      since they are always saved on the fine level instead of the coarse level.
	for (i=1;i<=slevel;i++)
	{
        smap=malloc(smesh[i-1].nt*sizeof(int *));
        for(j=0;j<smesh[i-1].nt;j++)
        {   smap[j]=malloc(tempsize2*sizeof(int));
            // tempsize2 is the initial length of the second index of "smap", and it may cause problem if it is too small.
            smap[j][0]=0;
            for(k=1;k<tempsize2;k++)
            {smap[j][k]=-1;}
        }
        spatialmapping(smesh[i-1],smesh[i],smap);

        for(j=0;j<smesh[i-1].nt;j++)
        {   if(smap[j][0]>tempsize2-1)
            {   printf("WARNING: tempsize2 for smap is too small!!\n");
                goto stop;
            }
        }

	    smesh[i].smap=malloc(smesh[i-1].nt*sizeof(int *));
        for(j=0;j<smesh[i-1].nt;j++)
        {   smesh[i].smap[j]=malloc((smap[j][0]+1)*sizeof(int));
            for(k=0;k<=smap[j][0];k++)
            {smesh[i].smap[j][k]=smap[j][k];}
        }

        for(j=0;j<smesh[i-1].nt;j++)
        {free(smap[j]);}
        free(smap);

        smesh[i].cf=malloc(smesh[i-1].nt*sizeof(double ***));
        for(j=0;j<smesh[i-1].nt;j++)
        {   smesh[i].cf[j]=malloc(3*sizeof(double **));
            for(k=0;k<3;k++)
            {   smesh[i].cf[j][k]=malloc(smesh[i].smap[j][0]*sizeof(double *));
                for(m=0;m<smesh[i].smap[j][0];m++)
                {smesh[i].cf[j][k][m]=malloc(3*sizeof(double));}
            }
        }
        smesh[i].fc=malloc(smesh[i-1].nt*sizeof(double ***));
        for(j=0;j<smesh[i-1].nt;j++)
        {   smesh[i].fc[j]=malloc(3*sizeof(double **));
            for(k=0;k<3;k++)
            {   smesh[i].fc[j][k]=malloc(smesh[i].smap[j][0]*sizeof(double *));
                for(m=0;m<smesh[i].smap[j][0];m++)
                {smesh[i].fc[j][k][m]=malloc(3*sizeof(double));}
            }
        }
        spatialmapping2(smesh[i-1],smesh[i],smesh[i].smap,smesh[i].cf,smesh[i].fc,tempsize2);
	}

    // 2.3. compute "e", "e2", "so2", "n" and "ori"
    //      For the data structure of "eo", "e2","so2", "n" and "ori", see "boundary".
    for(i=0;i<=slevel;i++)
    {
        smesh[i].e2=malloc(smesh[i].ne*sizeof(int *));
        for(j=0;j<smesh[i].ne;j++)
        {smesh[i].e2[j]=malloc(3*sizeof(int));}
        smesh[i].so2=malloc(smesh[i].nt*sizeof(int *));
        for(j=0;j<smesh[i].nt;j++)
        {smesh[i].so2[j]=malloc(3*sizeof(int));}
        smesh[i].n=malloc(smesh[i].ne*sizeof(double *));
        for(j=0;j<smesh[i].ne;j++)
        {smesh[i].n[j]=malloc(2*sizeof(double));}
        smesh[i].ori=malloc(smesh[i].ne*sizeof(int));

        boundary(smesh[i].ne,smesh[i].nt,smesh[i].t,smesh[i].p2,smesh[i].p,smesh[i].e,smesh[i].e2,smesh[i].so2,smesh[i].n,smesh[i].ori);
    }

    // 2.4. compute "bd" and "bd2"
    //      For the data structure of "bd" and "bd2", see "edgeterm".
    for(i=0;i<=slevel;i++)
    {   smesh[i].bd=malloc(amesh[alevel].ns*sizeof(int **));
        for(j=0;j<amesh[alevel].ns;j++)
        {   smesh[i].bd[j]=malloc(smesh[i].nt*sizeof(int *));
            for(k=0;k<smesh[i].nt;k++)
            {   smesh[i].bd[j][k]=malloc(9*sizeof(int));
                for(m=0;m<9;m++)
                {smesh[i].bd[j][k][m]=-1;}
            }
        }
    }
    for(i=0;i<=slevel;i++)
    {   smesh[i].bd2=malloc(amesh[alevel].ns*sizeof(double **));
        for(j=0;j<amesh[alevel].ns;j++)
        {   smesh[i].bd2[j]=malloc(smesh[i].nt*sizeof(double *));
            for(k=0;k<smesh[i].nt;k++)
            {smesh[i].bd2[j][k]=malloc(3*sizeof(double));}
        }
    }
    for(i=0;i<=slevel;i++)
    {   for(j=0;j<amesh[alevel].ns;j++)
        {edgeterm(smesh[i].nt,amesh[alevel].a[j],smesh[i].p,smesh[i].p2,smesh[i].t,smesh[i].bd[j],smesh[i].bd2[j],smesh[i].so2);}
    }

	// 3.1. compute "noflevel"
	//      level: the mesh layers for multgrid; the bigger value represents finer mesh.
	//      noflevel[level][2]: the corresponding angular mesh level and spatial mesh level for each multigrid mesh level
    //      Example:    assume noflevel[i][0]=3 and noflevel[i][1]=2, then
    //                  the spatial mesh level is "3" and the angular mesh level is "2" on the multgrid mesh level "i".
    da=alevel-alevel0;ds=slevel-slevel0;
    noflevel=malloc((level+1)*sizeof(int *));
	switch(whichmg)
	{
        case AMG:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            for(i=0;i<=da;i++)
            {
                noflevel[i][0]=slevel;
                noflevel[i][1]=i+alevel0;
            }
        break;
        case SMG:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            for(i=0;i<=ds;i++)
            {
                noflevel[i][0]=i+slevel0;
                noflevel[i][1]=alevel;
            }
        break;
        case MG1:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            if(ds>da)
            {   for(i=0;i<ds-da;i++)
                {
                    noflevel[i][0]=i+slevel0;
                    noflevel[i][1]=alevel0;
                }
                for(i=0;i<=da;i++)
                {
                    noflevel[i+ds-da][0]=i+ds-da+slevel0;
                    noflevel[i+ds-da][1]=i+alevel0;
                }
            }
            else
            {   for(i=0;i<da-ds;i++)
                {
                    noflevel[i][0]=slevel0;
                    noflevel[i][1]=i+alevel0;
                }
                for(i=0;i<=ds;i++)
                {
                    noflevel[i+da-ds][0]=i+slevel0;
                    noflevel[i+da-ds][1]=i+da-ds+alevel0;
                }
            }
        break;
        case MG2:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            for(i=slevel0;i<=slevel;i++)
            {
                noflevel[i-slevel0][0]=i;
                noflevel[i-slevel0][1]=alevel0;
            }
            for(i=alevel0+1;i<=alevel;i++)
            {
                noflevel[i-alevel0+slevel-slevel0][0]=slevel;
                noflevel[i-alevel0+slevel-slevel0][1]=i;
            }
        break;
        case MG3:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            for(i=alevel0;i<=alevel;i++)
            {
                noflevel[i-alevel0][0]=slevel0;
                noflevel[i-alevel0][1]=i;
            }
            for(i=slevel0+1;i<=slevel;i++)
            {
                noflevel[i-slevel0+alevel-alevel0][0]=i;
                noflevel[i-slevel0+alevel-alevel0][1]=alevel;
            }
        break;
        case MG4_a:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            if(ds>=da)
            {   for(i=0;i<=ds-da;i++)
                {
                    noflevel[i][0]=i+slevel0;
                    noflevel[i][1]=alevel0;
                }
                for(i=1;i<=da;i++)
                {
                    noflevel[ds-da+2*i-1][0]=slevel0+ds-da+i;
                    noflevel[ds-da+2*i-1][1]=alevel0+i-1;
                    noflevel[ds-da+2*i][0]=slevel0+ds-da+i;
                    noflevel[ds-da+2*i][1]=alevel0+i;
                }
            }
            else
            {   for(i=0;i<=da-ds;i++)
                {
                    noflevel[i][0]=slevel0;
                    noflevel[i][1]=i+alevel0;
                }
                for(i=1;i<=ds;i++)
                {
                    noflevel[da-ds+2*i-1][0]=slevel0+i;
                    noflevel[da-ds+2*i-1][1]=alevel0+da-ds+i-1;
                    noflevel[da-ds+2*i][0]=slevel0+i;
                    noflevel[da-ds+2*i][1]=alevel0+da-ds+i;
                }
            }
        break;
        case MG4_s:
            for(i=0;i<=level;i++)
            {noflevel[i]=malloc(2*sizeof(int));}
            if(ds>=da)
            {   for(i=0;i<=ds-da;i++)
                {
                    noflevel[i][0]=i+slevel0;
                    noflevel[i][1]=alevel0;
                }
                for(i=1;i<=da;i++)
                {
                    noflevel[ds-da+2*i-1][0]=slevel0+ds-da+i-1;
                    noflevel[ds-da+2*i-1][1]=alevel0+i;
                    noflevel[ds-da+2*i][0]=slevel0+ds-da+i;
                    noflevel[ds-da+2*i][1]=alevel0+i;
                }
            }
            else
            {   for(i=0;i<=da-ds;i++)
                {
                    noflevel[i][0]=slevel0;
                    noflevel[i][1]=i+alevel0;
                }
                for(i=1;i<=ds;i++)
                {
                    noflevel[da-ds+2*i-1][0]=slevel0+i-1;
                    noflevel[da-ds+2*i-1][1]=alevel0+da-ds+i;
                    noflevel[da-ds+2*i][0]=slevel0+i;
                    noflevel[da-ds+2*i][1]=alevel0+da-ds+i;
                }
            }
        break;
	}

    // 3.2. malloc "ua", "us", "flux", "RHS", "d" and "q"
    //      ua[level][nt]:absorption coefficient
    //      us[level][nt]:scattering coefficient
    //      flux[level][ns][nt][3]: photon flux
    //      RHS[level][ns][nt][3]: source term
    //      d[level][ns][nt][3]: residual
	{
    whichslevel = slevel;
    // ua
	ua=malloc((slevel+1)*sizeof(double *));
	for(i=0;i<=slevel;i++)
	{ua[i]=malloc(smesh[i].nt*sizeof(double));}

    //f = fopen("ua.txt", "r");
    //for(j=0;j<smesh[whichslevel].nt;j++)
    //{
    //    fscanf(f, "%le", &(ua[whichslevel][j]));
    //}
    //fclose(f);

    // us
	us=malloc((slevel+1)*sizeof(double *));
	for(i=0;i<=slevel;i++)
	{us[i]=malloc(smesh[i].nt*sizeof(double));}

    //f = fopen("us.txt", "r");
    //for(j=0;j<smesh[whichslevel].nt;j++)
    //{
    //    fscanf(f, "%le", &(us[whichslevel][j]));
    //}
    //fclose(f);
    

    //f = fopen("uaus.txt", "r");
    //int n_ua, n_us;
    //double *val_ua, *val_us;
    //fscanf(f, "%d", &n_ua);
    //val_ua = malloc((n_ua+1)*sizeof(double));
    //for (i = 0; i < n_ua + 1; ++i) 
    //{
    //   fscanf(f, "%le", &(val_ua[i])); 
    //}
    //fscanf(f, "%d", &n_us);
    //val_us = malloc((n_us+1)*sizeof(double));
    //for (i = 0; i < n_us + 1; ++i) 
    //{
    //   fscanf(f, "%le", &(val_us[i])); 
    //}
    //fclose(f);
    //f = fopen("shape.txt", "r");
    //int n_shape;
    //double **shape;
    //fscanf(f, "%d", &n_shape);
    //shape = malloc((n_shape)*sizeof(double *));
    //for (j = 0; j < n_shape; ++j) 
    //{
    //   shape[j] = malloc(5*sizeof(double));
    //   for(i = 0; i < 5; ++i)
    //   {
    //       fscanf(f, "%lf", &(shape[j][i]));
    //   }
    //}
    //fclose(f);

    double** shape = vars.shape;
    double x1, x2, y1, y2;
    for(j = 0; j < smesh[whichslevel].nt; j++)
    {
        //printf("%f\t%f\n", smesh[whichslevel].c[j][0], smesh[whichslevel].c[j][1]);
        ua[whichslevel][j] = vars.val_ua[0];
        us[whichslevel][j] = vars.val_us[0];
        for(i = 0; i < vars.n_ua; ++i)
        {
            x1 = smesh[whichslevel].c[j][0] - shape[i][0];
            x2 = smesh[whichslevel].c[j][1] - shape[i][1];
            y1 = x1 * cos(shape[i][4]) + x2 * sin(shape[i][4]);
            y2 = -x1 * sin(shape[i][4]) + x2 * cos(shape[i][4]);
            if(y1*y1/(shape[i][2]*shape[i][2]) + y2*y2/(shape[i][3]*shape[i][3]) < 1.)
              ua[whichslevel][j] = vars.val_ua[i+1];
        }
        for(i = vars.n_ua; i < vars.n_ua + vars.n_us; ++i)
        {
            x1 = smesh[whichslevel].c[j][0] - shape[i][0];
            x2 = smesh[whichslevel].c[j][1] - shape[i][1];
            y1 = x1 * cos(shape[i][4]) + x2 * sin(shape[i][4]);
            y2 = -x1 * sin(shape[i][4]) + x2 * cos(shape[i][4]);
            if(y1*y1/(shape[i][2]*shape[i][2]) + y2*y2/(shape[i][3]*shape[i][3]) < 1.)
            {
              us[whichslevel][j] = vars.val_us[i-vars.n_ua+1];
            }
        }
    }

	for(i=whichslevel+1;i<=slevel;i++)
	{ctof_s3(smesh[i-1].nt,ua[i],ua[i-1],smesh[i].smap);}
	for(i=whichslevel-1;i>=0;i--)
	{ftoc_s3(smesh[i].nt,ua[i+1],ua[i],smesh[i+1].smap,smesh[i+1].a);}
   
	for(i=whichslevel+1;i<=slevel;i++)
	{ctof_s3(smesh[i-1].nt,us[i],us[i-1],smesh[i].smap);}
	for(i=whichslevel-1;i>=0;i--)
	{ftoc_s3(smesh[i].nt,us[i+1],us[i],smesh[i+1].smap,smesh[i+1].a);}
    //RHS
	RHS=malloc((level+1)*sizeof(double ***));
	for(n=0;n<=level;n++)
	{   nt=smesh[noflevel[n][0]].nt;
        ns=amesh[noflevel[n][1]].ns;
	    RHS[n]=malloc(ns*sizeof(double **));
        for (i=0;i<ns;i++)
        {	RHS[n][i]=malloc(nt*sizeof(double *));
            for(j=0;j<nt;j++)
            {   RHS[n][i][j]=malloc(3*sizeof(double));
                for(k=0;k<3;k++)
                {   RHS[n][i][j][k]=0;
                }
            }
        }
	}
	d=malloc((level+1)*sizeof(double ***));
	for(n=0;n<=level;n++)
	{   nt=smesh[noflevel[n][0]].nt;
        ns=amesh[noflevel[n][1]].ns;
	    d[n]=malloc(ns*sizeof(double **));
        for (i=0;i<ns;i++)
        {	d[n][i]=malloc(nt*sizeof(double *));
            for(j=0;j<nt;j++)
            {   d[n][i][j]=malloc(3*sizeof(double));
                for(k=0;k<3;k++)
                {   d[n][i][j][k]=0;
                }
            }
        }
	}
	flux=malloc((level+1)*sizeof(double ***));
	for(n=0;n<=level;n++)
	{   nt=smesh[noflevel[n][0]].nt;
        ns=amesh[noflevel[n][1]].ns;
	    flux[n]=malloc(ns*sizeof(double **));
        for (i=0;i<ns;i++)
        {	flux[n][i]=malloc(nt*sizeof(double *));
            for(j=0;j<nt;j++)
            {   flux[n][i][j]=malloc(3*sizeof(double));
                for(k=0;k<3;k++)
                {   flux[n][i][j][k]=0;
                }
            }
        }
	}
	}

	// 3.3. compute "b"
    //      For the data structure of "b", see "boundarycoupling".
    b=malloc((level+1)*sizeof(struct boundarycoupling));
	if(vacuum==0)// we need "b" only in the presence of refraction index mismatch at the domain boundary.
	{
        for(i=0;i<=level;i++)
	    {
	        ne=smesh[noflevel[i][0]].ne;
	        ns=amesh[noflevel[i][1]].ns;

            b[i].ratio=malloc(ne*sizeof(double *));
            for(j=0;j<ne;j++)
            {b[i].ratio[j]=malloc(ns*sizeof(double));}

            b[i].nr=malloc(ne*sizeof(int *));
            for(j=0;j<ne;j++)
            {b[i].nr[j]=malloc(ns*sizeof(int));}
            b[i].r=malloc(ne*sizeof(int **));
            for(j=0;j<ne;j++)
            {   b[i].r[j]=malloc(ns*sizeof(int *));
                for(k=0;k<ns;k++)
                {b[i].r[j][k]=malloc(BSIZE*sizeof(int));}
            }
            b[i].r2=malloc(ne*sizeof(double **));
            for(j=0;j<ne;j++)
            {   b[i].r2[j]=malloc(ns*sizeof(double *));
                for(k=0;k<ns;k++)
                {b[i].r2[j][k]=malloc(BSIZE*sizeof(double));}
            }
            bc_reflection(ns,amesh[noflevel[i][1]].a,smesh[noflevel[i][0]],index_i,index_o,b[i]);
	    }
	}

    // 4. form "fp"
    fp.alevel=alevel;fp.alevel0=alevel0;fp.slevel=slevel;fp.slevel0=slevel0;
    fp.level=level;fp.whichmg=whichmg;fp.vacuum=vacuum;fp.index_i=index_i;fp.index_o=index_o;
    fp.noflevel=noflevel;fp.b=b;fp.flux=flux;fp.RHS=RHS;fp.d=d;fp.ua=ua;fp.us=us;
    fp.amesh=amesh;fp.smesh=smesh;
    stop:;
    return(fp);
}
