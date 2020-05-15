/* Program Gauss_Jacobi
   Solution of a system of linear equations by Gauss-Jacobi's iteration
   method. Testing of diagonal dominance is also incorporated */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>


#define N 3

void main()
{
	double a[N][N] = {{2,1,0},{1,2,1},{0,1,2}},
    b[N] = {3,4,3},
    x[N],xn[N],epp=0.00001,sum;
	int i,j,n,flag;
    double residual = 0.0,res; 

// 	printf("\nEnter number of variables: ");
// 	scanf("%d",&n);
    n = N;
// 	printf11("\nEnter the coefficients row-wise: ");
    
// 	for(i=0;i<n;i++)
// 	{
// 		for(j=0;j<n;j++)
// 		{
// 			scanf("%f",&a[i][j]);
// 		}
// 	}
// 	printf("\nEnter right hand vector: ");
// 	for(i=0;i<n;i++)
// 		scanf("%f",&b[i]);
// 	for(i=0;i<n;i++)
// 		x[i]=0; //initialize

// 	//checking for row dominance
// 	flag=0;
// 	for(i=0;i<n;i++)
// 	{
// 		sum=0;
// 		for(j=0;j<n;j++)
// 			if(i!=j)
// 				sum+=fabs(a[i][j]);
// 			if(sum>fabs(a[i][i]))
// 				flag=1;
// 	}
// 
// 	//checking for column dominance
// 	if(flag==1)
// 		flag=0;
// 	for(j=1;j<n;j++)
// 	{
// 		sum=0;
// 		for(i=1;i<n;i++)
// 			if(i!=j)
// 				sum+=fabs(a[i][j]);
// 			if(sum>fabs(a[j][j]))
// 				flag=1;
// 	}
// 
// 	if(flag==1)
// 	{
// 		printf("The coefficient matrix is not diagonally dominant \n");
// 		printf("The Gauss-Jacobi method does not converge surely");
// 		exit(0);
// 	}
// 
// 	for(i=0;i<n;i++)
// 		printf(" x[%d] ",i);
// 	printf("\n");

	do
	{
		for(i=0;i<n;i++)
		{
			sum=b[i];
			for(j=0;j<n;j++)
				if(j!=i)
					sum-=a[i][j]*x[j];
				xn[i]=sum/a[i][i];
		}
		for(i=0;i<n;i++)
			printf("%8.5f ",xn[i]);
		printf("\n");
		//double 
		residual = 0.0; 
		
		for(i=0;i<n;i++)
			if((res = fabs(x[i]-xn[i])) > residual) 
				residual = res;
		if(flag==1)
			for(i=1;i<n;i++)
				x[i]=xn[i]; //reset x[i]	

	}while(residual > 1e-3);

	printf("Solution is \n");
	for(i=0;i<n;i++)
		printf("%8.5f ",xn[i]);
}
