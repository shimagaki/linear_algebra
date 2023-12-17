#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;

extern "C" {
  void dsysv_ ( const char& UPLO,
		const int& N, const int& NRHS,
		double** A, const int& LDA, int* IPIV,
		double** B, const int& LDB,
		double* WORK, const int& LWORK,
		int& INFO, int UPLOlen );
};

// 実対称行列 A の線形方程式 A x = b を解く簡易関数
// Input: A[N][N], b[N]  Output: x[N]
// ただし、A[i][j<=i] の下三角しか参照されない
//
template <int N> int dsysv( double A[N][N], double x[N], double b[N] )
{
  int i, j, info;
  static int ipiv[N];
  static double work[4*N];
  static double U[N][N];

  for( j=0; j<N; j++ ){
    for( i=0; i<N; i++ ){
      U[i][j] = A[j][i];
    }
    x[j] = b[j];
  }

  dsysv_( 'L', N, 1, (double**)U, N, ipiv, (double**)x, N, work, 4*N, info, 1 );

  return info;
}


int main(void)
{
  const int N = 2500;

  double A[N][N];
  double x[N], b[N];

  for(int i=0; i<N; i++ ){
    for(int  j=i; j<N; j++ ){
      A[i][j] = 1+j;
      A[j][i] = 1+j;
    }//A[i][i]+=5;
    b[i] = 1.0;
  }

  int info = dsysv( A, x, b );

  //printf("# info=%d.\n", info );

  //printf("# Solution.\n");
  /*
  for(int i=0; i<N; i++ ){
    printf("%+f\n", x[i] );
  }
  */
double error=0.0;	
for(int i=0;i<N-1;++i){
error += pow(x[i],2);
}
error += pow((x[N-1]-1.0/N),2);

cout << "\n|x-x_0| = "<< sqrt(error/N) << endl;

/*
  printf("# Error.\n");
  for(int i=0; i<N; i++ ){
    double sum=0.0;
    for(int j=0; j<N; j++ ){
      sum += A[i][j]*x[j];
    }
    sum -= b[i];

    printf("%+f\n", sum );
  }
  */

  return 0;
}
