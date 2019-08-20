
double mod2pi(double x);

void bc_reflection(int ns,double **theta,struct spatialmesh smesh,double index_i,double index_o,struct boundarycoupling b)

// Purpose: this fucntion is to find the coupling relation between directions on the boundary
//          due to reflection and refraction in the presence of refraction index mismatch at the boundary.
//          For the data structure of "b", see "struct boundarycoupling" in "solver".

{
    int i,j,k,ne,**e,*ori,w[2],angle;
    double **n,**p,dx,dy,sn,dtheta=2*Pi/ns,w2[2];
    double theta_i,theta_reflection,ratio_reflection;
//    double theta_r,theta_refraction,ratio_refraction;
    ne=smesh.ne;p=smesh.p;n=smesh.n;e=smesh.e;ori=smesh.ori;

    for(i=0;i<ne;i++)
    {   dx=p[e[i][2]][0]-p[e[i][1]][0];
        dy=p[e[i][2]][1]-p[e[i][1]][1];
        if (ori[i]==0)
        {dx=-dx;dy=-dy;}// to make sure that (dx,dy) goes clockwisely.

        for(j=0;j<ns;j++)
        {b.ratio[i][j]=0;}
        for(j=0;j<ns;j++)
        {   b.nr[i][j]=0;
            for(k=0;k<BSIZE;k++)
            {b.r[i][j][k]=-1;b.r2[i][j][k]=0;}
        }

/*        for(j=0;j<ns;j++)
        {   for(k=0;k<2;k++)
            {b.s[i][j][k]=-1;b.s2[i][j][k]=0;}
        }*/

        for(j=0;j<ns;j++)
        {   sn=theta[j][0]*n[i][0]+theta[j][1]*n[i][1];// "s" dot "n" for the angle "s"
            if(sn>0)
            {   // contribution to incoming flux through internal reflection
                theta_i=acos(sn);
                ratio_reflection=reflection(theta_i,index_i,index_o);// reflection energy ratio of the total energy
                b.ratio[i][j]=1.0-ratio_reflection;

                if(theta[j][0]*dx+theta[j][1]*dy>0)// the ONLY place for clockwise (dx,dy)
                {theta_reflection=mod2pi(theta[j][2]+(Pi+2*theta_i));}
                else
                {theta_reflection=mod2pi(theta[j][2]+(Pi-2*theta_i));}

                intepolation_a(theta_reflection,dtheta,ns,w,w2,ratio_reflection);

                for(k=0;k<2;k++)
                {   if(w2[k]>1e-6)
                    {   angle=w[k];
                        b.r[i][angle][b.nr[i][angle]]=j;
                        b.r2[i][angle][b.nr[i][angle]]=w2[k];
                        b.nr[i][angle]++;
                    }
                }
                /*
                // contribution to outgoing flux after refraction: tracking-back method
                theta_r=theta_i;
                temp=sin(theta_r)*index_o/index_i;
                if(temp<1.0)
                {   theta_i=asin(temp);
                    ratio_refraction=1-reflection(theta_i,index_i,index_o);

                    if(theta[j][0]*dx+theta[j][1]*dy>0)// the ONLY place for clockwise (dx,dy)
                    {theta_refraction=mod2pi(theta[j][2]+(theta_r-theta_i));}
                    else
                    {theta_refraction=mod2pi(theta[j][2]-(theta_r-theta_i));}

                    intepolation_a(theta_refraction,dtheta,ns,b.s[i][j],b.s2[i][j],ratio_refraction);
                }*/
            }
        }

        for(j=0;j<ns;j++)
        {   if(b.nr[i][j]>BSIZE)
            {printf("BSIZE is too small\n");}
        }
    }
}

double reflection(double theta_i,double ni,double no)
// Purpose: this function is to find the reflection energy ratio for the incident angle "theta_i"
{
    double r,theta_t,temp1,temp2;
    if(abs2(theta_i)<1e-3)
    {
        temp1=(ni-no)/(ni+no);
        r=temp1*temp1;
    }
    else
    {
        temp1=sin(theta_i)*ni/no;
        if(temp1<1)
        {
            theta_t=asin(temp1);
            temp1=sin(theta_i-theta_t)/sin(theta_i+theta_t);
            temp2=tan(theta_i-theta_t)/tan(theta_i+theta_t);
            r=0.5*(temp1*temp1+temp2*temp2);
        }
        else
        {
            r=1;
        }
    }
    return r;
}

void intepolation_a(double theta_m,double dtheta,int ns,int *b,double *b2,double constant)
// Purpose: this function is to find two linearly intepolated angles "b" and weights "b2" for the angle "theta_m"
{
    int theta1,theta2;
    double w1,w2;

    theta1=(int)floor(theta_m/dtheta)+1;
    w2=(theta_m-(theta1-1)*dtheta)/dtheta;
    w1=1.0-w2;
    if(theta1==ns)
    {theta2=1;}
    else
    {theta2=theta1+1;}
    b[0]=theta1-1;b[1]=theta2-1;
    b2[0]=w1*constant;b2[1]=w2*constant;
}


double mod2pi(double x)
// Purpose: this function is to transfer angle "x" into [0 2*Pi)
{   double y;
    if (x<2*Pi)
    {   if(x<0)
        {y=x+2*Pi;}
        else
        {y=x;}
    }
    else
    {y=x-2*Pi;}
    return y;
}
