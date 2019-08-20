void ftoc_a(int nt,int ns_c,double ***flux,double ***cflux)
// Purpose: this function describes angular fine-to-coarse operator by simple restriction.
{
    int i,j,k;
    for (i=0;i<nt;i++)
    {	for(j=0;j<ns_c;j++)
        {   for(k=0;k<3;k++)
            {
                cflux[j][i][k]=flux[2*j][i][k];
            }
        }
    }
}

void ctof_a(int nt,int ns_c,double ***flux,double ***dcflux)
// Purpose: this function describes angular coarse-to-fine operator by linear interpolation.
{
    int i,j,k;
    for (i=0;i<nt;i++)
    {	for(j=0;j<ns_c-1;j++)
        {   for(k=0;k<3;k++)
            {
                flux[2*j][i][k]+=dcflux[j][i][k];
                flux[2*j+1][i][k]+=(dcflux[j][i][k]+dcflux[j+1][i][k])/2;
            }
        }
        j=ns_c-1;
        for(k=0;k<3;k++)
        {
            flux[2*j][i][k]+=dcflux[j][i][k];
            flux[2*j+1][i][k]+=(dcflux[j][i][k]+dcflux[0][i][k])/2;
        }
    }
}

void ftoc_s(int nt_c,int ns,double ***flux,double ***cflux,int **smap,double ****fc)
// Purpose: this function describes spatial fine-to-coarse operator by L2 projection for piecewise linear DG.
{
    int i,j,k,m,n;
    for (i=0;i<ns;i++)
    {	for(j=0;j<nt_c;j++)
        {   for(k=0;k<3;k++)
            {   cflux[i][j][k]=0;
                for(m=1;m<=smap[j][0];m++)
                {   for(n=0;n<3;n++)
                    {
                        cflux[i][j][k]+=flux[i][smap[j][m]][n]*fc[j][k][m-1][n];
                    }
                }
            }
        }
    }
}

void ctof_s(int nt_c,int ns,double ***flux,double ***dcflux,int **smap,double ****cf)
// Purpose: this function describes spatial coarse-to-fine operator by linear interpolation for piecewise linear DG.
{
    int i,j,k,m,n;
    for (i=0;i<ns;i++)
    {	for(j=0;j<nt_c;j++)
        {   for(k=0;k<3;k++)
            {   for(m=1;m<=smap[j][0];m++)
                {   for(n=0;n<3;n++)
                    {
                        flux[i][smap[j][m]][n]+=dcflux[i][j][k]*cf[j][k][m-1][n];
                    }
                }
            }
        }
    }
}

void ftoc_s3(int nt_c,double *f,double *cf,int **smap,double *area)
// Purpose: this function describes spatial fine-to-coarse operator by area-weighted averaging for piecewise constant DG.
{
    int i,j;
    double temp1,temp2;
    for (i=0;i<nt_c;i++)
    {	temp1=0;
        for(j=1;j<=smap[i][0];j++)
        {temp1+=area[smap[i][j]];}
        temp2=0;
        for(j=1;j<=smap[i][0];j++)
        {temp2+=f[smap[i][j]]*area[smap[i][j]];}
        cf[i]=temp2/temp1;
    }
}

void ctof_s3(int nt_c,double *f,double *cf,int **smap)
// Purpose: this function describes spatial coarse-to-fine operator by simple restriction for piecewise linear DG.
{
    int j,m;
    for(j=0;j<nt_c;j++)
    {   for(m=1;m<=smap[j][0];m++)
        {f[smap[j][m]]=cf[j];}
    }
}

void ftoc(int nt_c,int ns_c,double ***flux,double ***cflux,int **smap,double ****fc)
// Purpose: this function describes fine-to-coarse operator simultaneouly in angle and space,
//          i.e., the combination of "ftoc_a" and "ftoc_s".
{
    int i,j,k,m,n;
    for (i=0;i<ns_c;i++)
    {	for(j=0;j<nt_c;j++)
        {   for(k=0;k<3;k++)
            {   cflux[i][j][k]=0;
                for(m=1;m<=smap[j][0];m++)
                {   for(n=0;n<3;n++)
                    {
                        cflux[i][j][k]+=flux[2*i][smap[j][m]][n]*fc[j][k][m-1][n];
                    }
                }
            }
        }
    }
}

void ctof(int nt_c,int ns_c,double ***flux,double ***dcflux,int **smap,double ****cf)
// Purpose: this function describes coarse-to-fine operator simultaneouly in angle and space,
//          i.e., the combination of "ctof_a" and "ctof_s".
{
    int i,j,k,m,n;
    for (i=0;i<ns_c-1;i++)
    {	for(j=0;j<nt_c;j++)
        {   for(k=0;k<3;k++)
            {   for(m=1;m<=smap[j][0];m++)
                {   for(n=0;n<3;n++)
                    {
                        flux[2*i][smap[j][m]][n]+=dcflux[i][j][k]*cf[j][k][m-1][n];
                        flux[2*i+1][smap[j][m]][n]+=0.5*(dcflux[i][j][k]+dcflux[i+1][j][k])*cf[j][k][m-1][n];
                    }
                }
            }
        }
    }
    i=ns_c-1;
    for(j=0;j<nt_c;j++)
    {   for(k=0;k<3;k++)
        {   for(m=1;m<=smap[j][0];m++)
            {   for(n=0;n<3;n++)
                {
                    flux[2*i][smap[j][m]][n]+=dcflux[i][j][k]*cf[j][k][m-1][n];
                    flux[2*i+1][smap[j][m]][n]+=0.5*(dcflux[i][j][k]+dcflux[0][j][k])*cf[j][k][m-1][n];
                }
            }
        }
    }
}
