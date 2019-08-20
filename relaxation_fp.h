
void relaxation(int Ns,struct angularmesh amesh,struct spatialmesh smesh,double ***RHS,double *ua,double *us,
double ***flux,struct boundarycoupling bb,int vacuum)

// Purpose: this function is improved source-iteration (ISI) with vacuum or reflection boundary condition.

{
    int i,j,k,m,ii,jj,ns,nt,tri,bi,alevel,edge=-1;
    double left[3][3],right[3],temp[3],tempd;
    double dettri,cosi,sini,a,b,d,e,f,source_corr;
    double **w,*area,**theta,**p,***r2,***bd2;
    double bv[3][3]={{1.0/3,1.0/6,1.0/6},{1.0/6,1.0/3,1.0/6},{1.0/6,1.0/6,1.0/3}};
    double matrix1[3][3],matrix2[3][3];
    int **so,**t,**so2,**bnr,***r,***bd,**e2,index[3][2];

    ns=amesh.ns;w=amesh.w;theta=amesh.a;alevel=Ns/ns;
    nt=smesh.nt;so=smesh.so;so2=smesh.so2;area=smesh.a;p=smesh.p;t=smesh.t;bd=smesh.bd;bd2=smesh.bd2;e2=smesh.e2;
    bnr=bb.nr;r=bb.r;r2=bb.r2;

    index[0][0]=1;index[0][1]=2;
    index[1][0]=2;index[1][1]=0;
    index[2][0]=0;index[2][1]=1;

    if (vacuum==0)// reflection boundary condition
    {   for(i=0;i<ns;i++)
        {   bi=i*alevel;
            // "bi" is the angular index of the coarse angle on the fine angular mesh
            // since sweep ordering is saved on the finest angular mesh for each spatial mesh for simplicity.
            for(j=0;j<nt;j++)
            {   tri=so[bi][j];// the current sweep triangle by sweep ordering
                dettri=2*area[tri];
                cosi=theta[i][0];sini=theta[i][1];
                a=cosi*(p[t[tri][2]][1]-p[t[tri][0]][1])+sini*(p[t[tri][0]][0]-p[t[tri][2]][0]);
                b=cosi*(p[t[tri][0]][1]-p[t[tri][1]][1])+sini*(p[t[tri][1]][0]-p[t[tri][0]][0]);

                // matrix1: convection term
                matrix1[0][0]=(a+b)/6;matrix1[0][1]=matrix1[0][0];matrix1[0][2]=matrix1[0][0];
                matrix1[1][0]=-a/6;matrix1[1][1]=matrix1[1][0];matrix1[1][2]=matrix1[1][0];
                matrix1[2][0]=-b/6;matrix1[2][1]=matrix1[2][0];matrix1[2][2]=matrix1[2][0];

                // matrix2: absorption term

                a=ua[tri]+(1-w[i][i])*us[tri];
                matrix2[0][0]=dettri*a/12;matrix2[0][1]=dettri*a/24;matrix2[0][2]=dettri*a/24;
                matrix2[1][0]=matrix2[0][1];matrix2[1][1]=dettri*a/12;matrix2[1][2]=dettri*a/24;
                matrix2[2][0]=matrix2[0][2];matrix2[2][1]=matrix2[1][2];matrix2[2][2]=dettri*a/12;

//                a=ua[tri][0]+(1-w[i][i])*us[tri][0];
//                b=ua[tri][1]+(1-w[i][i])*us[tri][1];
//                c=ua[tri][2]+(1-w[i][i])*us[tri][2];
//                matrix2[0][0]=dettri/60*(3*a+b+c);matrix2[0][1]=dettri/120*(2*a+2*b+c);matrix2[0][2]=dettri/120*(2*a+b+2*c);
//                matrix2[1][0]=matrix2[0][1];matrix2[1][1]=dettri/60*(a+3*b+c);matrix2[1][2]=dettri/120*(a+2*b+2*c);
//                matrix2[2][0]=matrix2[0][2];matrix2[2][1]=matrix2[1][2];matrix2[2][2]=dettri/60*(a+b+3*c);

                // left[][]: convection+absorption
                for(ii=0;ii<3;ii++)
                {   for(jj=0;jj<3;jj++)
                    {left[ii][jj]=matrix1[ii][jj]+matrix2[ii][jj];}
                }

                // temp: scattering contribution to the direction "i" from all directions except "i".
                for(ii=0;ii<3;ii++)
                {   temp[ii]=0;
                    for(k=0;k<ns;k++)
                    {
                        temp[ii]+=w[i][k]*flux[k][tri][ii];
                    }
                    temp[ii]+=-w[i][i]*flux[i][tri][ii];
                }

                source_corr=area[tri]/12;

                a=us[tri];
//                a=us[tri][0];
//                b=us[tri][1];
//                c=us[tri][2];
                d=temp[0];e=temp[1];f=temp[2];
                // right[][]: scattering+light source
                // Note: point (delta) source needs correction.
                right[0]=dettri*a*(2*d+e+f)/24
                +(2*RHS[i][tri][0]+RHS[i][tri][1]+RHS[i][tri][2])*source_corr;
                right[1]=dettri*a*(d+2*e+f)/24
                +(RHS[i][tri][0]+2*RHS[i][tri][1]+RHS[i][tri][2])*source_corr;
                right[2]=dettri*a*(d+e+2*f)/24
                +(RHS[i][tri][0]+RHS[i][tri][1]+2*RHS[i][tri][2])*source_corr;
//                right[0]=dettri/120*(2*a*(3*d+e+f)+b*(2*d+2*e+f)+c*(2*d+e+2*f))
//                +(2*RHS[i][tri][0]+RHS[i][tri][1]+RHS[i][tri][2])*source_corr;
//                right[1]=dettri/120*(a*(2*d+2*e+f)+2*b*(d+3*e+f)+c*(d+2*e+2*f))
//                +(RHS[i][tri][0]+2*RHS[i][tri][1]+RHS[i][tri][2])*source_corr;
//                right[2]=dettri/120*(a*(2*d+e+2*f)+b*(d+2*e+2*f)+2*c*(d+e+3*f))
//                +(RHS[i][tri][0]+RHS[i][tri][1]+2*RHS[i][tri][2])*source_corr;

                // add edge contributions to left, or add upwind fluxes to right from boundary source or the adjacent triangle
                for(k=0;k<3;k++)
                {   // boundary outgoing flux (s dot n>0): add boundary contributions to left
                    if(bd2[bi][tri][k]>0)
                    {   for(ii=0;ii<2;ii++)
                        {   for(jj=0;jj<2;jj++)
                            {left[index[k][ii]][index[k][jj]]+=bd2[bi][tri][k]*bv[index[k][ii]][index[k][jj]];}
                        }
                    }
                    // boundary incoming flux (s dot n<0): add upwind fluxes to right
                    else if(bd2[bi][tri][k]<0)
                    {   if(so2[tri][k]>-1)
                        // "tri" is at the domain boundary: upwind flux is from internal reflection.
                        {   edge=so2[tri][k];

                            for(ii=0;ii<2;ii++)
                            {   temp[ii]=0;
                                for(m=0;m<bnr[edge][i];m++)
                                {
                                    temp[ii]+=flux[r[edge][i][m]][tri][e2[edge][ii]]*r2[edge][i][m];
                                }
                            }

                            for(ii=0;ii<2;ii++)
                            {   tempd=0;
                                for(jj=0;jj<2;jj++)
                                {
                                    tempd+=temp[jj]*bv[e2[edge][ii]][e2[edge][jj]];
                                }
                                right[e2[edge][ii]]+=-bd2[bi][tri][k]*tempd;
                            }
                        }
                        else
                        // "tri" is not at the domain boundary: upwind flux is from the adjacent triangle.
                        {   for(ii=0;ii<2;ii++)
                            {   tempd=0;
                                for(jj=0;jj<2;jj++)
                                {
                                    tempd+=flux[i][bd[bi][tri][3*k]][bd[bi][tri][3*k+jj+1]]*bv[index[k][ii]][index[k][jj]];
                                }
                                right[index[k][ii]]+=-bd2[bi][tri][k]*tempd;
                            }
                        }
                    }
                }
                matrixsolver(left,right,flux[i][tri]);// update the nodal values at "tri"
            }
        }
    }
    else // vacuum boundary condition
    {   for(i=0;i<ns;i++)
        {   bi=i*alevel;
            for(j=0;j<nt;j++)
            {   tri=so[bi][j];
                dettri=2*area[tri];
                cosi=theta[i][0];sini=theta[i][1];
                a=cosi*(p[t[tri][2]][1]-p[t[tri][0]][1])+sini*(p[t[tri][0]][0]-p[t[tri][2]][0]);
                b=cosi*(p[t[tri][0]][1]-p[t[tri][1]][1])+sini*(p[t[tri][1]][0]-p[t[tri][0]][0]);
                matrix1[0][0]=(a+b)/6;matrix1[0][1]=matrix1[0][0];matrix1[0][2]=matrix1[0][0];
                matrix1[1][0]=-a/6;matrix1[1][1]=matrix1[1][0];matrix1[1][2]=matrix1[1][0];
                matrix1[2][0]=-b/6;matrix1[2][1]=matrix1[2][0];matrix1[2][2]=matrix1[2][0];

                a=ua[tri]+(1-w[i][i])*us[tri];
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
                    {
                        temp[ii]+=w[i][k]*flux[k][tri][ii];
                    }
                    temp[ii]+=-w[i][i]*flux[i][tri][ii];
                }

                source_corr=area[tri]/12;

                a=us[tri];
                d=temp[0];e=temp[1];f=temp[2];
                right[0]=dettri*a*(2*d+e+f)/24
                +(2*RHS[i][tri][0]+RHS[i][tri][1]+RHS[i][tri][2])*source_corr;
                right[1]=dettri*a*(d+2*e+f)/24
                +(RHS[i][tri][0]+2*RHS[i][tri][1]+RHS[i][tri][2])*source_corr;
                right[2]=dettri*a*(d+e+2*f)/24
                +(RHS[i][tri][0]+RHS[i][tri][1]+2*RHS[i][tri][2])*source_corr;

                for(k=0;k<3;k++)
                {   if(bd2[bi][tri][k]>0)
                    {   for(ii=0;ii<2;ii++)
                        {   for(jj=0;jj<2;jj++)
                            {left[index[k][ii]][index[k][jj]]+=bd2[bi][tri][k]*bv[index[k][ii]][index[k][jj]];}
                        }
                    }
                    else if(bd2[bi][tri][k]<0&&so2[tri][k]<0)
                    // "tri" is not at the domain boundary
                    {   for(ii=0;ii<2;ii++)
                        {   tempd=0;
                            for(jj=0;jj<2;jj++)
                            {
                                tempd+=flux[i][bd[bi][tri][3*k]][bd[bi][tri][3*k+jj+1]]*bv[index[k][ii]][index[k][jj]];
                            }
                            right[index[k][ii]]+=-bd2[bi][tri][k]*tempd;
                        }
                    }
                }
                matrixsolver(left,right,flux[i][tri]);
            }
        }
    }
}
