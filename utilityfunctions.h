struct c_num{
    double r;//real
    double i;//imag
};

#define Pi 3.141592653589793                                                    // PI
#define AREA(x1,y1,x2,y2,x3,y3) 0.5*abs2((y2-y1)*(x3-x1)-(y3-y1)*(x2-x1))       // area of the triangle
#define LENGTH(x1,y1,x2,y2)     sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1))           // the edge length
#define BSIZE 4

#define AMG 1                                                                   // indexes for multgrid methods
#define SMG 2
#define MG1 3
#define MG2 4
#define MG3 5
#define MG4_a 6
#define MG4_s 7

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double abs2(double x);
double min(double x,double y);
double max(double x,double y);
int max2(int x,int y);
int findmin(int n,double *d);
int mod(int n, int p);
int computelevel(int alevel,int alevel0,int slevel,int slevel0,int whichmg);
double dotproduct(int N, double *x, double *y);
int findtri(int tri,int *pt1,int *pt2);
int findtri2(int *pt1,int *pt2);

double** generatePoints(int, double, int);

// other utility functions
int intpower(int x,int y)
// Purpose: this function is to compute x to the power y.
{
	int i,res;
	res=x;
	for(i=1;i<y;i++)
	{
		res=x*res;
	}
	return res;
}

double abs2(double x)
// Purpose: this function is to compute the absolute value of x.
{return x>=0?x:(-x);}

double min(double x,double y)
// Purpose: this function is to find the maximum of two numbers x and y.
{return x<y?x:y;}

double max(double x,double y)
// Purpose: this function is to find the maximum of two numbers x and y.
{return x>y?x:y;}

int max2(int x,int y)
// Purpose: this function is to find the maximum of two integers x and y.
{return x>y?x:y;}

double dotproduct(int N, double *x, double *y)
{
    int i=0;
    double res=0;

    for(i=0;i<N;i++)
    {res+=x[i]*y[i];}

    return res;
}

int mod(int n,int p)
// Purpose: this function is to compute mod(n,p).
{
    int temp=(int)(n/p);

    if(n==temp*p)
    {return 0;}
    else
    {return 1;}
}

int computelevel(int alevel,int alevel0,int slevel,int slevel0,int whichmg)
{
    int ds=slevel-slevel0,da=alevel-alevel0,level=-1;
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
	return level;
}

int findmin(int n,double *d)
// Purpose: this function is to find the minimum from the vector d with the size n.
{
    int i;
    double dmin;
    int nmin;
    nmin=0;dmin=d[nmin];
    for(i=1;i<n;i++)
    {   if(d[i]<dmin)
        {nmin=i;dmin=d[i];}
    }
    return nmin;
}

int findtri2(int *pt1,int *pt2)
// Purpose: this function is to find the triangle containing the boundary.
// Note: "pt1" and "pt2" have to be boundary nodes of the same boundary.
{
    int i,j,tri=-1;
    for(i=1;i<=pt1[0];i++)
    {   for(j=1;j<=pt2[0];j++)
        {   if (pt1[i]==pt2[j])
            {
                tri=pt1[i];
                goto stop;
            }
        }
    }
    stop:return tri;
}

int findtri(int tri,int *pt1,int *pt2)
// Purpose: this function is find the adjacent triangle of the current triangle "tri"
// Note: the edge can not be at the boundary of the domain.
{
    int i,j,tri2=-1,flag=0;
    for(i=1;i<=pt1[0];i++)
    {   for(j=1;j<=pt2[0];j++)
        { if (pt1[i]==pt2[j])
            {   if(pt1[i]!=tri)
                {
                    flag=1;
                    tri2=pt1[i];
                    goto stop;
                }
            }
        }
    }
    if(flag==0)
    {printf("nonconforming mesh!!!\n");}
    stop:return tri2;
}


double** generatePoints(int n, double R, int flag)
{
    int i, j;
    double angle = 2. * Pi / (double)n;
    double start = 0.;
    if(flag == 0)
    {
        start = angle / 2.;
    }
    // start += 2. * Pi / (192.*0.125);
    // start += angle / 16.;
    double **var;
    var=malloc(n*sizeof(double *));
    for(i = 0; i < n; i++)
    {   
        var[i]=malloc(3*sizeof(double));
        for(j = 0; j < 3; j++)
        {
            var[i][0] = R * cos(start + (double)i * angle);
            var[i][1] = R * sin(start + (double)i * angle);
            var[i][2] = -(start + (double)i * angle);
        }
    }
    return var;
}
