// #include "mex.h"
// struct detector det_fp(struct variables_fp fp,int SOR,const mxArray *para_sd)
struct detector det_fp(struct variables_fp fp, struct matlab_variables vars, int SOR)

// Purpose: this function is to write measurements to six ".txt" files given the input file "det.txt".
//          Input:
//              "det.txt": "B", "A", "n", and "coord".

{
    struct detector det;
    int vacuum,alevel,slevel,level;
    int i,j,k;
    int ne,nt,ns,**e,**t,**e2,w_a[2],tri=-1,edge;
    double **p,**n,**ec,**c,*a,**theta,**refraction,ratio_refraction,*distance,w_s[3],w_a2[2],l[2],area[3];
    double index_i,index_o,x,y,x1,x2,y1,y2,x3,y3,sn,dtheta,temp,**dxy;

   	alevel=fp.alevel;slevel=fp.slevel;level=fp.level;
   	refraction=fp.b[level].ratio;

   	vacuum=fp.vacuum;index_i=fp.index_i;index_o=fp.index_o;
   	ns=fp.amesh[alevel].ns;theta=fp.amesh[alevel].a;
   	ne=fp.smesh[slevel].ne;nt=fp.smesh[slevel].nt;n=fp.smesh[slevel].n;p=fp.smesh[slevel].p;
   	e=fp.smesh[slevel].e;e2=fp.smesh[slevel].e2;ec=fp.smesh[slevel].ec;t=fp.smesh[slevel].t;c=fp.smesh[slevel].c;a=fp.smesh[slevel].a;

   	dtheta=2*Pi/ns;

    // load "det.txt"
   	//FILE *f = fopen("det.txt", "r");
    //int Bsor, Asor, Nsor, Bdet, Adet, Ndet;
    //fscanf(f, "%d", &Bsor);
    //fscanf(f, "%d", &Asor);
    //fscanf(f, "%d", &Nsor);
    //fscanf(f, "%d", &Bdet);
    //fscanf(f, "%d", &Adet);
    //fscanf(f, "%d", &Ndet);
    //fclose(f);
    if(SOR==1)
    {   
        det.B = vars.Bsor;
        det.A = vars.Asor;
        det.n = vars.Nsor;
        // printf("slevel = %d\n", slevel);
        dxy = generatePoints(det.n, vars.R, 1, slevel);
        // printf("Source: B = %d, A = %d, n = %d. Points:\n", det.B, det.A, det.n);
        // for (int i = 0; i < det.n; ++i) 
        // {
        //     printf("%f\t%f\t%f\n", dxy[i][0], dxy[i][1], dxy[i][2]); 
        // }
    }
    else
    {   
        det.B = vars.Bdet;
        det.A = vars.Adet;
        det.n = vars.Ndet;
        dxy = generatePoints(det.n, vars.R, 0, slevel);
        // printf("Detector: B = %d, A = %d, n = %d. Points:\n", det.B, det.A, det.n);
        // for (int i = 0; i < det.n; ++i) 
        // {
        //     printf("%f\t%f\t%f\n", dxy[i][0], dxy[i][1], dxy[i][2]); 
        // }
    }

    det.t=malloc(det.n*sizeof(int));
    det.i=malloc(det.n*sizeof(double **));
    for(i=0;i<det.n;i++)
    {   det.i[i]=malloc(ns*sizeof(double *));
        for(j=0;j<ns;j++)
        {   det.i[i][j]=malloc(3*sizeof(double));
            for(k=0;k<3;k++)
            {det.i[i][j][k]=0;}
        }
    }

    // given "dxy", compute and save "output"

    if(det.B==1)
    {
    // boundary measurement
    distance=malloc(ne*sizeof(double));
    for(i=0;i<det.n;i++)
    {   // nodal interpolation at the boundary
        for(j=0;j<ne;j++)
        {   distance[j]=0;
            temp=0;
            for(k=0;k<2;k++)
            {   temp=dxy[i][k]-ec[j][k];
                distance[j]+=temp*temp;
            }
        }
        x=dxy[i][0];y=dxy[i][1];
        edge=findmin(ne,distance);
        // printf("%d\t%f\t", edge, distance[edge]);
        x1=p[e[edge][1]][0];y1=p[e[edge][1]][1];
        x2=p[e[edge][2]][0];y2=p[e[edge][2]][1];
        l[0]=LENGTH(x,y,x2,y2);
        l[1]=LENGTH(x1,y1,x,y);
        // printf("%f\t%f\t%e\t", l[0], l[1], l[0] / (l[1]+1.e-12) - 1);
        if(fabs(l[0] / (l[1] + 1.e-12) - 1) < 1.e-4) { // added by fyw
            temp = (l[0] + l[1]) * 0.5;
            l[0] = l[1] = temp;
            // printf("haha\t");
        }
        // printf("%f\t%f\t%e\t", l[0], l[1], l[0] / (l[1]+1.e-12) - 1);
        // printf("\n");
        temp=l[0]+l[1];
        for (k=0;k<2;k++)
        {w_s[k]=l[k]/temp;}

        tri=e[edge][0];
        det.t[i]=tri;

        //for(j=0;j<ns;j++){
        //    printf("theta\t%f\t%f\n", theta[j][0], theta[j][1]);
        //}
        // printf("n\t%.12e\t%.12e\n", n[edge][0], n[edge][1]);
        if(det.A==1)
        {   
            for(j=0;j<ns;j++)// angular integration with weights "s dot n" on the angular set with " sn>0"
            {   
                sn=theta[j][0]*n[edge][0]+theta[j][1]*n[edge][1];
                if(SOR==1&&sn<=0)
                {   //if(sn<-0.8)
                    {   for(k=0;k<2;k++)// boundary interpolation
                        {det.i[i][j][e2[edge][k]]+=(-sn)*w_s[k];}
//                      {det.i[i][j][e2[edge][k]]+=w_s[k];}
                    }
                }
                if(SOR==0&&sn>0)
                {   if(vacuum==1)
                    //if(1)
                    {   for(k=0;k<2;k++)// boundary interpolation
                        {det.i[i][j][e2[edge][k]]+=sn*w_s[k]*dtheta;}
                    }
                    else
                    {   for(k=0;k<2;k++)// boundary interpolation
                        {det.i[i][j][e2[edge][k]]+=sn*w_s[k]*dtheta*refraction[edge][j];}
                    }
                }
            }
        }
        else
        {   sn=cos(dxy[i][2])*n[edge][0]+sin(dxy[i][2])*n[edge][1];
            if(SOR==1&&sn<=0)
            {   intepolation_a(dxy[i][2],dtheta,ns,w_a,w_a2,1.0);
                for(k=0;k<2;k++)// boundary interpolation
                {   for(j=0;j<2;j++)// angular interpolation
                    {det.i[i][w_a[j]][e2[edge][k]]+=w_a2[j]*w_s[k];}
                }
            }
            if(SOR==0&&sn>0)
            {   intepolation_a(dxy[i][2],dtheta,ns,w_a,w_a2,1.0);
                if(vacuum==1)
                {   for(k=0;k<2;k++)// boundary interpolation
                    {   for(j=0;j<2;j++)// angular interpolation
                        {det.i[i][w_a[j]][e2[edge][k]]+=w_a2[j]*w_s[k];}
                    }
                }
                else
                {   ratio_refraction=1-reflection(acos(sn),index_i,index_o);
                    for(k=0;k<2;k++)// boundary interpolation
                    {   for(j=0;j<2;j++)// angular interpolation
                        {det.i[i][w_a[j]][e2[edge][k]]+=w_a2[j]*w_s[k]*ratio_refraction;}
                    }
                }
            }
        }
    }
    free(distance);
    }
    else
    {   distance=malloc(nt*sizeof(double));
        for(i=0;i<det.n;i++)
        {   // spatial interpolation
            for(j=0;j<nt;j++)
            {   distance[j]=0;
                temp=0;
                for(k=0;k<2;k++)
                {   temp=dxy[i][k]-c[j][k];
                    distance[j]+=temp*temp;
                }
            }
            x=dxy[i][0];y=dxy[i][1];
            for(j=0;j<nt;j++)
            {   tri=findmin(nt,distance);
                distance[tri]=1e10;
                x1=p[t[tri][0]][0];y1=p[t[tri][0]][1];
                x2=p[t[tri][1]][0];y2=p[t[tri][1]][1];
                x3=p[t[tri][2]][0];y3=p[t[tri][2]][1];
                area[0]=AREA(x,y,x2,y2,x3,y3);
                area[1]=AREA(x1,y1,x,y,x3,y3);
                area[2]=AREA(x1,y1,x2,y2,x,y);
                temp=area[0]+area[1]+area[2];
                if(abs2(temp-a[tri])/a[tri]<1e-2)
                {
                    for(k=0;k<3;k++)
                    {w_s[k]=area[k]/temp;}
                    goto stop;
                }
            }
            stop:;

            det.t[i]=tri;
            for(j=0;j<ns;j++)// angular integration with weights "s dot n" on the angular set with " sn>0"
            {   for(k=0;k<3;k++)// boundary interpolation
                {   if(SOR==1)
                    {det.i[i][j][k]+=w_s[k];}
                    else
                    {det.i[i][j][k]+=w_s[k]*dtheta;}
                }
            }
        }
        free(distance);
    }


    for(i=0;i<det.n;i++)
    {free(dxy[i]);}
    free(dxy);

    return(det);
}
