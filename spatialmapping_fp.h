
void spatialmapping(struct spatialmesh cmesh,struct spatialmesh fmesh,int **smap)

// Purpose: this function is to find the spatial mapping between coarse and fine mesh, i.e.,
//          the number and the correpsonding fine triangles contained in the coarse triangle.
//
//          smap[nt_c][smap[nt_c][0]+1]:
//              nt_c is the number of triangles in the coarse mesh;
//              smap[nt_c][0] is the number of fine triangles contained in the corresponding coarse triangle;
//              the global nodal order of fine triangles are saved from the 2nd entry to the end.
//          Example:    assume the fine triangles "10 24 38" are contained in the coarse triangle "4", then
//                      smap[4] has 3 elements with smap[4][0]=3, smap[4][1]=10, smap[4][2]=24, and smap[4][3]=38.

{
    double **p,**c_f,**c_c,*a,*distance;
    double area1,area2,area3,x,y,x1,y1,x2,y2,x3,y3;
    int **t,nt_f,nt_c,i,j,tri,flag;
    p=cmesh.p;t=cmesh.t;c_c=cmesh.c;c_f=fmesh.c;a=cmesh.a;
    nt_f=fmesh.nt;nt_c=cmesh.nt;

    if(nt_f==nt_c)
    {   for(i=0;i<nt_c;i++)
        {smap[i][0]=1;smap[i][1]=i;}
    }
    else
    {
        distance=malloc(nt_c*sizeof(double));

        for(i=0;i<nt_c;i++)
        {smap[i][0]=0;}

        for(i=0;i<nt_f;i++)
        {
            flag=0;
            x=c_f[i][0];y=c_f[i][1];
            for(j=0;j<nt_c;j++)
            {distance[j]=(x-c_c[j][0])*(x-c_c[j][0])+(y-c_c[j][1])*(y-c_c[j][1]);}
            for(j=0;j<nt_c;j++)
            {
                tri=findmin(nt_c,distance);// find the triangle with minimal distance between (x,y) and centers of coarse triangles.
                distance[tri]=1e10;
                x1=p[t[tri][0]][0];y1=p[t[tri][0]][1];
                x2=p[t[tri][1]][0];y2=p[t[tri][1]][1];
                x3=p[t[tri][2]][0];y3=p[t[tri][2]][1];

                area1=AREA(x,y,x2,y2,x3,y3);
                area2=AREA(x1,y1,x,y,x3,y3);
                area3=AREA(x1,y1,x2,y2,x,y);
                if(abs2(area1+area2+area3-a[tri])/a[tri]<1e-2)// If true, then the fine triangle "i" is in the coarse triangle "j".
                {
                    flag=1;
                    smap[tri][0]+=1;
                    smap[tri][smap[tri][0]]=i;
                    goto stop;
                }
            }
            if(flag==0)
            {printf("spatial mapping has a problem\n");}
            stop:;
        }
        free(distance);
    }
}

void spatialmapping2(struct spatialmesh cmesh,struct spatialmesh fmesh,int **smap,double ****cf,double ****fc,int tempsize)

// Purpose: this function is to compute "cf" and "fc".
//
//          cf[nt_c][3][smap[nt_c][0]][3]: coarse-to-fine spatial mapping of nodal values computed by linear interpolation
//          fc[nt_c][3][smap[nt_c][0]][3]: fine-to-coarse spatial mapping of nodal values computed by L2 projection
//
//          the 1st indexs the coarse triangles;
//          the 2nd indexes 3 coarse nodes in the coarse triangle;
//          the 3rd indexes the fine triangles contained in the coarse triangle;
//          the 4th indexes 3 fine nodes in the fine triangle.

{
    int **ind_p,np,**temp,nt_c,**tf,flag,i,j,k,m,**t;
    double **temp2,tempd,x1,x2,x3,y1,y2,y3,x,y,**pf,**p;
    double xf[3],yf[3],*a,*af;
    tf=fmesh.t;pf=fmesh.p;af=fmesh.a;
    nt_c=cmesh.nt;p=cmesh.p;t=cmesh.t;a=cmesh.a;

    temp=malloc(3*tempsize*sizeof(int *));
    for(i=0;i<3*tempsize;i++)
    {   temp[i]=malloc(2*sizeof(int));
        for(j=0;j<2;j++)
        {temp[i][j]=-1;}
    }
    temp2=malloc(3*tempsize*sizeof(double *));
    for(i=0;i<3*tempsize;i++)
    {   temp2[i]=malloc(3*sizeof(double));
        for(j=0;j<3;j++)
        {temp2[i][j]=0;}
    }
    ind_p=malloc(tempsize*sizeof(int *));
    for(j=0;j<tempsize;j++)
    {   ind_p[j]=malloc(3*sizeof(int));
        for(k=0;k<3;k++)
        {ind_p[j][k]=-1;}
    }

    for(i=0;i<nt_c;i++)
    {
        // compute "ind_p" and "temp"
        // temp[][]:    find discrete (without repeat) nodes of the fine triangles contained in the coarse triangle
        //              and save it sequentially, i.e., order it in sequence as 1,2,3... and so on.
        //              the 1st indexes fine nodes;
        //              for the 2nd index, 1st entry is the global nodal index on fine mesh, 2nd entry is the repeatence of this node.

        // ind_p[][3]:  the local nodal correspondence in fine triangles contained in the coarse triangle
        //              Example:    if the 2nd node of the 3rd fine triangle is the 4th node as in "temp",
        //                          then ind_p[3][2]=4.

        temp[0][0]=tf[smap[i][1]][0];
        temp[0][1]=1;
        np=1;ind_p[0][0]=0;
        j=0;
        {   for(k=1;k<3;k++)
            {   flag=0;// flag indicates whether it is a new node that has not been found yet.
                for(m=0;m<np;m++)
                {   if(tf[smap[i][j+1]][k]==temp[m][0])// existing node
                    {   ind_p[j][k]=m;
                        temp[m][1]+=1;
                        flag=1;
                        goto stop;
                    }
                }
                stop:;
                if(flag==0)// new node
                {   ind_p[j][k]=np;
                    temp[np][0]=tf[smap[i][j+1]][k];
                    temp[np][1]=1;
                    np+=1;
                }
            }
        }

        for(j=1;j<smap[i][0];j++)
        {   for(k=0;k<3;k++)
            {   flag=0;
                for(m=0;m<np;m++)
                {   if(tf[smap[i][j+1]][k]==temp[m][0])
                    {   ind_p[j][k]=m;
                        temp[m][1]+=1;
                        flag=1;
                        goto stop2;
                    }
                }
                stop2:;
                if(flag==0)
                {   ind_p[j][k]=np;
                    temp[np][0]=tf[smap[i][j+1]][k];
                    temp[np][1]=1;
                    np+=1;
                }
            }
        }

        // compute temp2
        // temp2[][]: the nodal weights of non-repeat fine nodes in the coarse triangle.
        x1=p[t[i][0]][0];y1=p[t[i][0]][1];
        x2=p[t[i][1]][0];y2=p[t[i][1]][1];
        x3=p[t[i][2]][0];y3=p[t[i][2]][1];
        tempd=AREA(x1,y1,x2,y2,x3,y3);
        for(j=0;j<np;j++)
        {
            x=pf[temp[j][0]][0];y=pf[temp[j][0]][1];
            temp2[j][1]=AREA(x1,y1,x,y,x3,y3)/tempd;
            temp2[j][2]=AREA(x1,y1,x2,y2,x,y)/tempd;
            temp2[j][0]=1-temp2[j][1]-temp2[j][2];
        }

        // compute "cf"
//        tempd=0;
        for(j=0;j<smap[i][0];j++)
        {   for(k=0;k<3;k++)
            {   for(m=0;m<3;m++)
                {
                    cf[i][m][j][k]=temp2[ind_p[j][k]][m];
//            tempd+=cf[i][m][j][k];
                }
            }
        }
//        printf("%f\n",tempd);

        // compute "fc"
        for(j=0;j<smap[i][0];j++)
        {
            for(k=0;k<3;k++)
            {   xf[k]=temp2[ind_p[j][k]][1];
                yf[k]=temp2[ind_p[j][k]][2];
            }

            fc[i][1][j][0]=2*xf[0]+xf[1]+xf[2]-1;
            fc[i][2][j][0]=2*yf[0]+yf[1]+yf[2]-1;
            fc[i][0][j][0]=1-fc[i][1][j][0]-fc[i][2][j][0];

            fc[i][1][j][1]=xf[0]+2*xf[1]+xf[2]-1;
            fc[i][2][j][1]=yf[0]+2*yf[1]+yf[2]-1;
            fc[i][0][j][1]=1-fc[i][1][j][1]-fc[i][2][j][1];

            fc[i][1][j][2]=xf[0]+xf[1]+2*xf[2]-1;
            fc[i][2][j][2]=yf[0]+yf[1]+2*yf[2]-1;
            fc[i][0][j][2]=1-fc[i][1][j][2]-fc[i][2][j][2];

            tempd=af[smap[i][j+1]]/a[i];
            for(k=0;k<3;k++)
            {   for(m=0;m<3;m++)
                {fc[i][m][j][k]*=tempd;}
            }
        }
    }
    for(i=0;i<3*tempsize;i++)
    {free(temp[i]);free(temp2[i]);}
    free(temp);free(temp2);
    for(i=0;i<tempsize;i++)
    {free(ind_p[i]);}
    free(ind_p);
}
