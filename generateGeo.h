#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 
#include <math.h>

#ifndef Pi
#define Pi 3.141592653589793
#endif 


double myRand();
double** generateGeo(int n, double len, int n_msh);
int isShapeAvailable(double** shapes, int n, int n_msh);
int isInEllipse(double x1, double x2, double* shape);
int isInCircle(double* shape, int n_msh);
int isIntersect(double* shape1, double* shape2, int n_msh);
double myRand()
{
    return (double)rand() / (double)RAND_MAX;
}
int count = 0;
double** generateGeo(int n, double len, int n_msh)
{
    double** shape;
    double r, theta;
    double len2;
    len2 = len * 2. / 3.;
    int i;
    shape = malloc(n*sizeof(double*));
    for(i = 0; i < n; ++i) 
    {
        shape[i] = malloc(5*sizeof(double));
    }
    while(1) {
        for(i = 0; i < n; ++i)
        {
            r = (0.9-1.2*len) * myRand();
            theta = 2. * Pi * myRand();
            shape[i][0] = r * cos(theta);
            shape[i][1] = r * sin(theta);
            shape[i][2] = len * myRand() + len;
            shape[i][3] = len2 * myRand() + len2/2.;
            shape[i][4] = 2. * Pi * myRand();
        }
        count++;
        if(isShapeAvailable(shape, n, n_msh)) 
        {
            //printf("%d\n", count);
            return shape;
        }
    }
}

int isShapeAvailable(double** shapes, int n, int n_msh)
{
    int i, j;
    int flag;
    double factor = 1.2;
    for(i = 0; i < n; ++i)
    {
        for(j = i+1; j < n; ++j)
        {
            if( (shapes[i][0]-shapes[j][0])*(shapes[i][0]-shapes[j][0]) + (shapes[i][1]-shapes[j][1])*(shapes[i][1]-shapes[j][1]) < 0.25 )
              return 0;
        }
    }
    flag = 1;
    for(i = 0; i < n; ++i)
    {
        shapes[i][2] *= factor;
        shapes[i][3] *= factor;
    }
    for(i = 0; i < n; ++i)
    {
        if(!isInCircle(shapes[i], n_msh)) { flag = 0; }
    }
    for(i = 0; i < n; ++i)
    {
        for(j = i+1; j < n; ++j)
        {
            if(flag && isIntersect(shapes[i], shapes[j], n_msh))
            {
                flag = 0;
            }
        }
    }
    for(i = 0; i < n; ++i)
    {
        shapes[i][2] /= factor;
        shapes[i][3] /= factor;
    }
    return flag;
}

int isInCircle(double* shape, int n_msh)
{
    double c1, c2, a, b, theta;
    double x1, x2, y1, y2, angle;
    int i;
    double r = 0.9;
    c1 = shape[0];
    c2 = shape[1];
    a = shape[2];
    b = shape[3];
    theta = shape[4];
    if(sqrt(c1*c1+c2*c2) + min(a, b) > r) return 0;

    angle = 2. * Pi / (double)n_msh;
    for(i = 0; i < n_msh; ++i)
    {
        x1 = a * cos(angle * i);
        x2 = b * sin(angle * i);
        y1 = x1 * cos(theta) - x2 * sin(theta) + c1;
        y2 = x1 * sin(theta) + x2 * cos(theta) + c2;
        if(sqrt(y1*y1+y2*y2) > r) return 0;
    }
    return 1;
}
int isIntersect(double* shape1, double* shape2, int n_msh)
{
    double x1, x2, y1, y2, angle;
    int i;
    angle = 2. * Pi / (double)n_msh;
    for(i = 0; i < n_msh; ++i)
    {
        x1 = shape1[2] * cos(angle * i);
        x2 = shape1[3] * sin(angle * i);
        y1 = x1 * cos(shape1[4]) - x2 * sin(shape1[4]) + shape1[0];
        y2 = x1 * sin(shape1[4]) + x2 * cos(shape1[4]) + shape1[1];
        if(isInEllipse(y1, y2, shape2)) return 1;
    }
    return 0;
}

int isInEllipse(double x1, double x2, double* shape)
{
    double y1, y2;
    y1 = x1 - shape[0];
    y2 = x2 - shape[1];
    x1 = y1 * cos(shape[4]) + y2 * sin(shape[4]);
    x2 = -y1 * sin(shape[4]) + y2 * cos(shape[4]);
    if(x1*x1/(shape[2]*shape[2]) + x2*x2/(shape[3]*shape[3]) < 1) return 1;
    return 0;
}
