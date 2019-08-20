void free_space_fp(struct variables_fp fp)

// Purpose: this function is to free variables malloced in initialization.

{
    int i,j,k,m,n;

    int alevel,slevel,level,vacuum;
	struct boundarycoupling *b;
    struct angularmesh *amesh;
    struct spatialmesh *smesh;
	double ****flux,****RHS,****d;
	double **ua,**us;
	int **noflevel;

    alevel=fp.alevel;slevel=fp.slevel;level=fp.level;vacuum=fp.vacuum;
    noflevel=fp.noflevel;b=fp.b;flux=fp.flux;RHS=fp.RHS;d=fp.d;ua=fp.ua;us=fp.us;
    amesh=fp.amesh;smesh=fp.smesh;

    if(vacuum==0)
    {   for(i=0;i<=level;i++)
	    {   for(j=0;j<smesh[noflevel[i][0]].ne;j++)
            {   for(k=0;k<amesh[noflevel[i][1]].ns;k++)
                {   free(b[i].r[j][k]);free(b[i].r2[j][k]);
                }
                free(b[i].ratio[j]);free(b[i].nr[j]);
                free(b[i].r[j]);free(b[i].r2[j]);
            }
            free(b[i].ratio);free(b[i].nr);
            free(b[i].r);free(b[i].r2);
	    }
	}
    free(b);

    for(n=0;n<=level;n++)
    {   for(i=0;i<amesh[noflevel[n][1]].ns;i++)
        {   for(j=0;j<smesh[noflevel[n][0]].nt;j++)
            {free(RHS[n][i][j]);free(d[n][i][j]);free(flux[n][i][j]);}
            free(RHS[n][i]);free(d[n][i]);free(flux[n][i]);
        }
        free(RHS[n]);free(d[n]);free(flux[n]);free(noflevel[n]);
    }
    free(RHS);free(flux);free(d);free(noflevel);

    for(n=0;n<=slevel;n++)
    {free(ua[n]);free(us[n]);}
    free(ua);free(us);

    for(i=0;i<=slevel;i++)
    {   for(j=0;j<amesh[alevel].ns;j++)
        {   for(k=0;k<smesh[i].nt;k++)
            {free(smesh[i].bd[j][k]);free(smesh[i].bd2[j][k]);}
            free(smesh[i].bd[j]);free(smesh[i].bd2[j]);
        }
        free(smesh[i].bd);free(smesh[i].bd2);
    }

    for(i=0;i<=slevel;i++)
    {   for(j=0;j<smesh[i].nt;j++)
        {free(smesh[i].so2[j]);}
        free(smesh[i].so2);
        for(j=0;j<smesh[i].ne;j++)
        {free(smesh[i].n[j]);free(smesh[i].e2[j]);}
        free(smesh[i].n);free(smesh[i].e2);free(smesh[i].ori);
    }

	for (i=1;i<=slevel;i++)
	{	for(j=0;j<smesh[i-1].nt;j++)
		{   for(k=0;k<3;k++)
            {   for(m=0;m<smesh[i].smap[j][0];m++)
                {free(smesh[i].cf[j][k][m]);free(smesh[i].fc[j][k][m]);}
                free(smesh[i].cf[j][k]);free(smesh[i].fc[j][k]);
            }
            free(smesh[i].cf[j]);free(smesh[i].fc[j]);free(smesh[i].smap[j]);
		}
		free(smesh[i].cf);free(smesh[i].fc);free(smesh[i].smap);
	}

	for (i=0;i<=slevel;i++)
	{	for(j=0;j<smesh[i].nt;j++)
		{free(smesh[i].t[j]);free(smesh[i].c[j]);}
        free(smesh[i].t);free(smesh[i].c);
		for(j=0;j<smesh[i].ne;j++)
		{free(smesh[i].e[j]);free(smesh[i].ec[j]);}
        free(smesh[i].e);free(smesh[i].ec);
		for(j=0;j<smesh[i].np;j++)
		{free(smesh[i].p[j]);free(smesh[i].p2[j]);}
		free(smesh[i].p);free(smesh[i].p2);
		for(j=0;j<amesh[alevel].ns;j++)
		{free(smesh[i].so[j]);}
		free(smesh[i].so);free(smesh[i].a);
	}

	for (i=0;i<=alevel;i++)
	{	for(j=0;j<amesh[i].ns;j++)
		{free(amesh[i].w[j]);free(amesh[i].a[j]);}
        free(amesh[i].a);free(amesh[i].w);
    }

    free(smesh);free(amesh);
}
