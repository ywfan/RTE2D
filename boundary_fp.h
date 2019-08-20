
void boundary(int ne,int nt,int **t,int **p2,double **p,int **e,int **e2,int **so2,double **n,int *ori)

// Purpose: this function is to compute "e", "e2", "so2", "n", "sn" and "ori" with the data structure as follow:
//
//          e[ne][4]:   for the 2nd index,1st is the corresponding triangle containing boundary,
//                      2nd and 3rd are global nodal index of boundary nodes, and 4th
//                      is global index of the other non-boundary node in the triangle.
//
//          e2[ne][3]:  the local order of the boundary nodes in the corresponding triangle with
//                      1st entry for node 1 (e[ne][1]), 2nd entry for node 2 (e[ne][2])
//                      3rd entry for the local order of the boundary in the triangle.
//
//          so2[nt][3]: the corresponding triangle and edge of the boundary
//                      Example:    if boundary "10" is edge "2" of triangle "20",
//                                  then so2[20][2]=10, and "-1" otherwise.
//
//          n[ne][2]:   x and y coordinate of the outgoing normal.
//
//          ori[ne]:    orientation of the boundary edge w.r.t. the triangle, which is useful at ONLY one place in "bc_reflection".

{
    int i,j,k,tri;
    double nx,ny,temp;

    for(i=0;i<nt;i++)
    {   for(j=0;j<3;j++)
        {so2[i][j]=-1;}
    }
    for(i=0;i<ne;i++)
    {
        tri=findtri2(p2[e[i][1]],p2[e[i][2]]);
        e[i][0]=tri;

        // find the local index for boundary nodes in the corresponding triangle
        for(j=0;j<2;j++)
        {   for(k=0;k<3;k++)
            {   if(e[i][j+1]==t[tri][k])
                {e2[i][j]=k;goto stop;}
            }
            stop:;
        }
        // find the non-boundary node and the boundary index in the corresponding triangle
        for(j=0;j<3;j++)
        {   if(t[tri][j]!=e[i][1]&&t[tri][j]!=e[i][2])
            {   e[i][3]=t[tri][j];//non-edge node
                so2[tri][j]=i;//boundary index
                e2[i][2]=j;//local order of the boundary
                goto stop2;
            }
        }
        stop2:;
        // The following is to compute outgoing vector "n" on the boundary
        nx=-(p[e[i][1]][1]-p[e[i][2]][1]);
        ny=p[e[i][1]][0]-p[e[i][2]][0];
        ori[i]=1;
        if(nx*(p[e[i][3]][0]-p[e[i][2]][0])+ny*(p[e[i][3]][1]-p[e[i][2]][1])>0)// the normal is incoming.
        {nx=-nx;ny=-ny;ori[i]=0;}
        temp=sqrt(nx*nx+ny*ny);
        n[i][0]=nx/temp;n[i][1]=ny/temp;
    }
}
