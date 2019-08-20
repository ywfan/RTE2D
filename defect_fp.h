
void defect(int Ns,struct angularmesh amesh,struct spatialmesh smesh,double ***RHS,double *ua,double *us,double ***flux,
struct boundarycoupling bb,double ***res,int vacuum)

// Purpose: this function is to compute the residual with vacuum or reflection boundary condition.
//          see "relaxation" for more details.

{
    int i,j,k,m,ii,jj,ns,nt,bi,alevel,edge=0;
    double left[3][3],right[3],temp[3],tempd;
    double dettri,cosi,sini,a,b,d,e,f,source_corr;
    double **w,*area,**theta,**p,***r2,***bd2;
    double bv[3][3]={{1.0/3,1.0/6,1.0/6},{1.0/6,1.0/3,1.0/6},{1.0/6,1.0/6,1.0/3}};
    double matrix1[3][3],matrix2[3][3];
    //int **so;
    int **t,**so2,**bnr,***r,**e2,***bd,index[3][2];

    ns=amesh.ns;w=amesh.w;theta=amesh.a;alevel=Ns/ns;
    nt=smesh.nt;
    //so=smesh.so;
    so2=smesh.so2;area=smesh.a;p=smesh.p;t=smesh.t;bd=smesh.bd;bd2=smesh.bd2;e2=smesh.e2;
    bnr=bb.nr;r=bb.r;r2=bb.r2;

    index[0][0]=1;index[0][1]=2;
    index[1][0]=2;index[1][1]=0;
    index[2][0]=0;index[2][1]=1;

    if (vacuum==0)
    {   for(i=0;i<ns;i++)
        {   bi=i*alevel;
            for(j=0;j<nt;j++)
            {   dettri=2*area[j];
                cosi=theta[i][0];sini=theta[i][1];
                a=cosi*(p[t[j][2]][1]-p[t[j][0]][1])+sini*(p[t[j][0]][0]-p[t[j][2]][0]);
                b=cosi*(p[t[j][0]][1]-p[t[j][1]][1])+sini*(p[t[j][1]][0]-p[t[j][0]][0]);
                matrix1[0][0]=(a+b)/6;matrix1[0][1]=matrix1[0][0];matrix1[0][2]=matrix1[0][0];
                matrix1[1][0]=-a/6;matrix1[1][1]=matrix1[1][0];matrix1[1][2]=matrix1[1][0];
                matrix1[2][0]=-b/6;matrix1[2][1]=matrix1[2][0];matrix1[2][2]=matrix1[2][0];

                a=ua[j]+us[j];
                matrix2[0][0]=dettri*a/12;matrix2[0][1]=dettri*a/24;matrix2[0][2]=dettri*a/24;
                matrix2[1][0]=matrix2[0][1];matrix2[1][1]=dettri*a/12;matrix2[1][2]=dettri*a/24;
                matrix2[2][0]=matrix2[0][2];matrix2[2][1]=matrix2[1][2];matrix2[2][2]=dettri*a/12;

                for(ii=0;ii<3;ii++)
                {   for(jj=0;jj<3;jj++)
                    {left[ii][jj]=matrix1[ii][jj]+matrix2[ii][jj];}
                }

                for(ii=0;ii<3;ii++)
                {   temp[ii]=0;
                    for(k=0;k<ns;k++)
                    {temp[ii]+=w[i][k]*flux[k][j][ii];}
                }

                source_corr=area[j]/12;
                a=us[j];
                d=temp[0];e=temp[1];f=temp[2];
                right[0]=dettri*a*(2*d+e+f)/24
                +(2*RHS[i][j][0]+RHS[i][j][1]+RHS[i][j][2])*source_corr;
                right[1]=dettri*a*(d+2*e+f)/24
                +(RHS[i][j][0]+2*RHS[i][j][1]+RHS[i][j][2])*source_corr;
                right[2]=dettri*a*(d+e+2*f)/24
                +(RHS[i][j][0]+RHS[i][j][1]+2*RHS[i][j][2])*source_corr;

                for(k=0;k<3;k++)
                {   if(bd2[bi][j][k]>0)
                    {   for(ii=0;ii<2;ii++)
                        {   for(jj=0;jj<2;jj++)
                            {left[index[k][ii]][index[k][jj]]+=bd2[bi][j][k]*bv[index[k][ii]][index[k][jj]];}
                        }
                    }
                    else if(bd2[bi][j][k]<0)
                    {   if(so2[j][k]>-1)// upwind flux from the internal reflection and boundary source
                        {   edge=so2[j][k];

                            for(ii=0;ii<2;ii++)
                            {   temp[ii]=0;
                                for(m=0;m<bnr[edge][i];m++)
                                {
                                    temp[ii]+=flux[r[edge][i][m]][j][e2[edge][ii]]*r2[edge][i][m];
                                }
                            }

                            for(ii=0;ii<2;ii++)
                            {   tempd=0;
                                for(jj=0;jj<2;jj++)
                                {
                                    tempd+=temp[jj]*bv[e2[edge][ii]][e2[edge][jj]];
                                }
                                right[e2[edge][ii]]+=-bd2[bi][j][k]*tempd;
                            }
                        }
                        else // upwind flux from the adjacent triangle
                        {   for(ii=0;ii<2;ii++)
                            {   tempd=0;
                                for(jj=0;jj<2;jj++)
                                {
                                    tempd+=flux[i][bd[bi][j][3*k]][bd[bi][j][3*k+jj+1]]*bv[index[k][ii]][index[k][jj]];
                                }
                                right[index[k][ii]]+=-bd2[bi][j][k]*tempd;
                            }
                        }
                    }
                }
                // compute "right-left*flux" and then find the corresponding nodal values of "res".
                tempd=0;
                for(ii=0;ii<3;ii++)
                {
                    temp[ii]=(right[ii]-left[ii][0]*flux[i][j][0]-left[ii][1]*flux[i][j][1]-left[ii][2]*flux[i][j][2])/source_corr;
                    tempd+=temp[ii];
                }
                tempd=tempd/4;
                for(ii=0;ii<3;ii++)
                {res[i][j][ii]=temp[ii]-tempd;}
            }
        }
    }
    else
    {   for(i=0;i<ns;i++)
        {   bi=i*alevel;
            for(j=0;j<nt;j++)
            {   dettri=2*area[j];
                cosi=theta[i][0];sini=theta[i][1];
                a=cosi*(p[t[j][2]][1]-p[t[j][0]][1])+sini*(p[t[j][0]][0]-p[t[j][2]][0]);
                b=cosi*(p[t[j][0]][1]-p[t[j][1]][1])+sini*(p[t[j][1]][0]-p[t[j][0]][0]);
                matrix1[0][0]=(a+b)/6;matrix1[0][1]=matrix1[0][0];matrix1[0][2]=matrix1[0][0];
                matrix1[1][0]=-a/6;matrix1[1][1]=matrix1[1][0];matrix1[1][2]=matrix1[1][0];
                matrix1[2][0]=-b/6;matrix1[2][1]=matrix1[2][0];matrix1[2][2]=matrix1[2][0];

                a=ua[j]+us[j];
                matrix2[0][0]=dettri*a/12;matrix2[0][1]=dettri*a/24;matrix2[0][2]=dettri*a/24;
                matrix2[1][0]=matrix2[0][1];matrix2[1][1]=dettri*a/12;matrix2[1][2]=dettri*a/24;
                matrix2[2][0]=matrix2[0][2];matrix2[2][1]=matrix2[1][2];matrix2[2][2]=dettri*a/12;

                for(ii=0;ii<3;ii++)
                {   for(jj=0;jj<3;jj++)
                    {left[ii][jj]=matrix1[ii][jj]+matrix2[ii][jj];}
                }

                for(ii=0;ii<3;ii++)
                {   temp[ii]=0;
                    for(k=0;k<ns;k++)
                    {temp[ii]+=w[i][k]*flux[k][j][ii];}
                }

                source_corr=area[j]/12;
                a=us[j];
                d=temp[0];e=temp[1];f=temp[2];
                right[0]=dettri*a*(2*d+e+f)/24
                +(2*RHS[i][j][0]+RHS[i][j][1]+RHS[i][j][2])*source_corr;
                right[1]=dettri*a*(d+2*e+f)/24
                +(RHS[i][j][0]+2*RHS[i][j][1]+RHS[i][j][2])*source_corr;
                right[2]=dettri*a*(d+e+2*f)/24
                +(RHS[i][j][0]+RHS[i][j][1]+2*RHS[i][j][2])*source_corr;

                for(k=0;k<3;k++)
                {   if(bd2[bi][j][k]>0)
                    {   for(ii=0;ii<2;ii++)
                        {   for(jj=0;jj<2;jj++)
                            {left[index[k][ii]][index[k][jj]]+=bd2[bi][j][k]*bv[index[k][ii]][index[k][jj]];}
                        }
                    }
                    else if(bd2[bi][j][k]<0&&so2[j][k]<0)
                    {   for(ii=0;ii<2;ii++)
                        {   tempd=0;
                            for(jj=0;jj<2;jj++)
                            {
                                tempd+=flux[i][bd[bi][j][3*k]][bd[bi][j][3*k+jj+1]]*bv[index[k][ii]][index[k][jj]];
                            }
                            right[index[k][ii]]+=-bd2[bi][j][k]*tempd;
                        }
                    }
                }
                // compute "right-left*flux" and then find the corresponding nodal values of "res".
                tempd=0;
                for(ii=0;ii<3;ii++)
                {
                    temp[ii]=(right[ii]-left[ii][0]*flux[i][j][0]-left[ii][1]*flux[i][j][1]-left[ii][2]*flux[i][j][2])/source_corr;
                    tempd+=temp[ii];
                }
                tempd=tempd/4;
                for(ii=0;ii<3;ii++)
                {res[i][j][ii]=temp[ii]-tempd;}
            }
        }
    }
}
