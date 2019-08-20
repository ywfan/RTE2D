
double mgcycle(struct angularmesh *amesh,struct spatialmesh *smesh,
struct boundarycoupling *b,double ****RHS,double **ua,double **us,double ****flux,double ****d,
int n1,int n2,int alevel,int alevel0,int slevel,int slevel0,int NS,int vacuum,int whichmg)

// Purpose: this function contains the multigrid methods with V-cycle.
//     AMG: angular multigrid method.
//     SMG: spatial multigrid method.
//     MG1: simultaneou angular and spatial multigrid method.
//     MG2: angular multigrid first, then spatial multigrid.
//     MG3: spatial multigrid first, then angular multigrid.
//     MG4: alternative angular and spatial multigrid.
//          MG4_a: angular multigrid first; MG4_s: spatial multigrid first.

{
    int i,nt,ns,da,ds,level=-1;
    double res=1e10;
    nt=smesh[slevel].nt;ns=amesh[alevel].ns;
    ds=slevel-slevel0;da=alevel-alevel0;

	switch(whichmg)
	{   case AMG:
            level=da;
        break;
        case SMG:
            level=ds;
        break;
        case MG1:
            level=max2(da,ds);
        break;
        case MG2:
            level=ds+da;
        break;
        case MG3:
            level=ds+da;
        break;
        case MG4_a:
            level=ds+da;
        break;
        case MG4_s:
            level=ds+da;
        break;
	}

    for(i=0;i<n1;i++)
    {relaxation(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],vacuum);}

	switch(whichmg)
	{   case AMG:
        {   if(alevel==alevel0)
            {}
            else
            {
                defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                res=residual(nt,ns,d[level],smesh[slevel].a);
                ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
            }
        }
        break;
        case SMG:
        {   if(slevel==slevel0)
            {}
            else
            {
                defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                res=residual(nt,ns,d[level],smesh[slevel].a);
                ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
            }
        }
        break;
        case MG1:
        {   if(alevel==alevel0)
            {   if(slevel==slevel0)
                {}
                else
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                    ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
                }
            }
            else
            {
                if(slevel==slevel0)
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                    ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
                }
                else
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc(smesh[slevel-1].nt,amesh[alevel-1].ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                    ctof(smesh[slevel-1].nt,amesh[alevel-1].ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
                }
            }
        }
        break;
        case MG2:
        {   if(alevel==alevel0)
            {   if(slevel==slevel0)
                {}
                else
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                    ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
                }
            }
            else
            {
                defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                res=residual(nt,ns,d[level],smesh[slevel].a);
                ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
            }
        }
        break;
        case MG3:
        {   if(slevel==slevel0)
            {   if(alevel==alevel0)
                {}
                else
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                    ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
                }
            }
            else
            {
                defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                res=residual(nt,ns,d[level],smesh[slevel].a);
                ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
            }
        }
        break;
        case MG4_a:
        {   if(alevel==alevel0)
            {   if(slevel==slevel0)
                {}
                else
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                    ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
                }
            }
            else
            {   if(slevel==slevel0)
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                    ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
                }
                else
                {   whichmg=MG4_s;
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                    ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
                }
            }
        }
        break;
        case MG4_s:
        {   if(alevel==alevel0)
            {   if(slevel==slevel0)
                {}
                else
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                    ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
                }
            }
            else
            {   if(slevel==slevel0)
                {
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_a(nt,amesh[alevel-1].ns,d[level],RHS[level-1]);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel-1,alevel0,slevel,slevel0,NS,vacuum,whichmg);
                    ctof_a(nt,amesh[alevel-1].ns,flux[level],flux[level-1]);
                }
                else
                {   whichmg=MG4_a;
                    defect(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],d[level],vacuum);
                    res=residual(nt,ns,d[level],smesh[slevel].a);
                    ftoc_s(smesh[slevel-1].nt,ns,d[level],RHS[level-1],smesh[slevel].smap,smesh[slevel].fc);
                    mgcycle(amesh,smesh,b,RHS,ua,us,flux,d,n1,n2,alevel,alevel0,slevel-1,slevel0,NS,vacuum,whichmg);
                    ctof_s(smesh[slevel-1].nt,ns,flux[level],flux[level-1],smesh[slevel].smap,smesh[slevel].cf);
                }
            }
        }
        break;
	}


    for(i=0;i<n2;i++)
    {relaxation(NS,amesh[alevel],smesh[slevel],RHS[level],ua[slevel],us[slevel],flux[level],b[level],vacuum);}
    return res;
}

double residual(int nt,int ns,double ***d,double *a)
// Purpose: this function is to compute the residual.
{
    double res=0,temp;
    int i,j,k;

    for (i=0;i<nt;i++)
    {	temp=0;
        for(j=0;j<ns;j++)
        {   for(k=0;k<3;k++)
            {temp+=abs2(d[j][i][k]);}
        }
        res+=temp*a[i];
    }
    return res*2*Pi/ns/3;
}

