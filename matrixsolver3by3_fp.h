void matrixsolver(double A[3][3],double B[3],double sol[3])
// Purpose: this function is to solve a 3 by 3 nonsingular linear system by Gaussian elimination with pivoting.
{
    int ind[3]={0,1,2},i,j,temp2;
    double temp;

    if(A[0][0]==0)
    {
        if(A[1][0]==0)
        {ind[0]=2;ind[2]=0;}
        else
        {ind[0]=1;ind[1]=0;}
    }

    temp=A[ind[0]][0];
    B[ind[0]]=B[ind[0]]/temp;
    for(j=0;j<3;j++)
    {A[ind[0]][j]=A[ind[0]][j]/temp;}

    for(i=1;i<3;i++)
    {   if(A[ind[i]][0]==0)
        {}
        else
        {   temp=A[ind[i]][0];
            B[ind[i]]=B[ind[i]]/temp-B[ind[0]];
            for(j=0;j<3;j++)
            {A[ind[i]][j]=A[ind[i]][j]/temp-A[ind[0]][j];}
        }
    }

    if(A[ind[1]][1]==0)
    {temp2=ind[1];ind[1]=ind[2];ind[2]=temp2;}
    temp=A[ind[1]][1];
    B[ind[1]]=B[ind[1]]/temp;
    for(j=1;j<3;j++)
    {A[ind[1]][j]=A[ind[1]][j]/temp;}

    if(A[ind[2]][1]==0)
    {}
    else
    {   temp=A[ind[2]][1];
        B[ind[2]]=B[ind[2]]/temp-B[ind[1]];
        for(j=1;j<3;j++)
        {A[ind[2]][j]=A[ind[2]][j]/temp-A[ind[1]][j];}
    }

    sol[2]=B[ind[2]]/A[ind[2]][2];
    sol[1]=(B[ind[1]]-sol[2]*A[ind[1]][2])/A[ind[1]][1];
    sol[0]=(B[ind[0]]-sol[2]*A[ind[0]][2]-sol[1]*A[ind[0]][1])/A[ind[0]][0];
}
