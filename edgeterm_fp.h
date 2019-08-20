
void edgeterm(int nt,double theta[3],double **p,int **p2,int **t,int **bd,double **bd2,int **so2)

// Purpose: this function is to compute edge integrals for each triangle.
//
//          Let "n" be the normal of some edge of the triangle and "s" be the angular direction with
//          the assumption for conforming mesh that one edge is shared by at most two triangles.
//
//          Case 1 (s dot n>0): the edge flux is outgoing from the triangle.
//          Case 2 (s dot n<0):
//              Case 2.1 (internal triangle): the edge flux is incoming from the adjacent triangle.
//              Case 2.2 (boundary triangle): the edge flux is incoming from the boundary source.
//
//          For each angular direction "s" on each mesh "slevel", the edge integrals can be assembled from the following:
//
//          bd[nt][9]:  in Case 2.1, it saves the adjacent triangle of the current triangle, the local order of shared nodes
//                      in the adjacent triangle in the order of edges: index "0" to "2" for edge "1" in the current triangle,
//                      index "3" to "5" for edge "2" and index "6" to "8" for edge "3".
//                      Example:    assume if the upwind flux for edge "2" of triangle "50" comes from triangle "10" with
//                                  the local order "1" and "3", then bd[50][3]=10, bd[50][4]=1 and bd[50][5]=3.
//                      in Case 2.2, the corresponding boundary source is found through the function "boundary".
//
//          bd2[nt][3]: the values of "s dot n * L (length of the edge)" by the order of edges of the current triangle.
//
// Note:    we define "23" as the 1st edge, "31" as the 2nd and "12" as the 3rd.

{
    int i,j;
    double a=theta[0],b=theta[1]; // a=cos(theta), b=sin(theta)
    double x1=0,y1=0,x2=0,y2=0,x3=0,y3=0,dx,dy,tol=1e-6;
    for(i=0;i<nt;i++)
    {
        x1=p[t[i][0]][0];y1=p[t[i][0]][1];
        x2=p[t[i][1]][0];y2=p[t[i][1]][1];
        x3=p[t[i][2]][0];y3=p[t[i][2]][1];

        // 1st edge "23"
        dx=y3-y2;dy=x2-x3;
        if(dx*(x1-x3)+dy*(y1-y3)>0) // to make sure the normal of edge "23" is outgoing w.r.t. the triangle.
        {dx=-dx;dy=-dy;}

        bd2[i][0]=a*dx+b*dy;// s dot n * L
        if(abs2(bd2[i][0])<tol)
        {bd2[i][0]=0;}
        if(bd2[i][0]<0&&so2[i][0]==-1)
        // "bd2[i][n]<0 or s dot n <0" means that this edge has upwind flux from the other adjacent triangle or the boundary source.
        // "so2[i][0]==-1" indicates that this triangle is not at the boundary.
        {   bd[i][0]=findtri(i,p2[t[i][1]],p2[t[i][2]]);// to find the adjacent triangle
            if(bd[i][0]>-1)
            {   for(j=0;j<3;j++)
                {   if(t[bd[i][0]][j]==t[i][1])// to find the local order in the adjacent triangle of node "2" of the current triangle
                    {bd[i][1]=j;goto stop3;}
                }
                stop3:;
                for(j=0;j<3;j++)
                {   if(t[bd[i][0]][j]==t[i][2])// to find the local order in the adjacent triangle of node "3" of the current triangle
                    {bd[i][2]=j;goto stop4;}
                }
                stop4:;
            }
        }

        // 2nd edge "31"
        dx=y1-y3;dy=x3-x1;
        if(dx*(x2-x3)+dy*(y2-y3)>0)
        {dx=-dx;dy=-dy;}

        bd2[i][1]=a*dx+b*dy;
        if(abs2(bd2[i][1])<tol)
        {bd2[i][1]=0;}
        if(bd2[i][1]<0&&so2[i][1]==-1)
        {   bd[i][3]=findtri(i,p2[t[i][0]],p2[t[i][2]]);
            if(bd[i][3]>-1)
            {   for(j=0;j<3;j++)
                {   if(t[bd[i][3]][j]==t[i][2])
                    {bd[i][4]=j;goto stop5;}
                }
                stop5:;
                for(j=0;j<3;j++)
                {   if(t[bd[i][3]][j]==t[i][0])
                    {bd[i][5]=j;goto stop6;}
                }
                stop6:;
            }
        }

        // 3rd edge "12"
        dx=y2-y1;dy=x1-x2;
        if(dx*(x3-x1)+dy*(y3-y1)>0)
        {dx=-dx;dy=-dy;}

        bd2[i][2]=a*dx+b*dy;
        if(abs2(bd2[i][2])<tol)
        {bd2[i][2]=0;}
        if(bd2[i][2]<0&&so2[i][2]==-1)
        {   bd[i][6]=findtri(i,p2[t[i][0]],p2[t[i][1]]);
            if(bd[i][6]>-1)
            {   for(j=0;j<3;j++)
                {   if(t[bd[i][6]][j]==t[i][0])
                    {bd[i][7]=j;goto stop1;}
                }
                stop1:;
                for(j=0;j<3;j++)
                {   if(t[bd[i][6]][j]==t[i][1])
                    {bd[i][8]=j;goto stop2;}
                }
                stop2:;
            }
        }
    }
}
