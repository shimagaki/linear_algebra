#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <vector>  //std::back_inserter
using namespace std;
using namespace Eigen;
int L = 5000;
int n_iter = 10000;
double epsilon_pAp = 1e-30;
double epsilon_rr = 1e-30;
/**** Ax = y ****/
MatrixXd A = MatrixXd::Zero(L,L);
VectorXd y,x0,x,p0,p,r0,r;
double a,b;

void set_A_y(){
	y = VectorXd::Ones(L);
	for(int i=0;i<L;++i){
	for(int j=i;j<L;++j){
		A(i,j) = j+1;	
		A(j,i) = j+1;	
	}}
}

void init_vector(){
	x0= VectorXd::Zero(L);
	r0 = y;
	p0 = y;
}

double norm_x(){
	double error=0;	
	for(int i=0;i<L-1;++i){
	error += pow(x[i],2);
	}
	error += pow((x[L-1]-1.0/L),2);
	return sqrt(error/L);
}

int main(int argc, char const* argv[]){
	set_A_y();
	//cout << "\nA = \n" << A << endl;
	//cout << "\ny = \n" << y << endl;
	init_vector();
	double r0r0, p0Ap0, rr;	
	cout<<endl;	
	
	for(int n=0;n<n_iter;++n){
		r0r0 = r0.dot(r0);
		p0Ap0 = p0.dot(A * p0);
		
		if(abs(p0Ap0)<epsilon_pAp){
			cout<< "A too small variables appears in p0.dot(A * p0) " <<endl;
			break;
		}	
		
		a = r0r0/ p0Ap0; 
		x = x0 + a * p0;
		r = r0 - a * (A*p0);
		rr = r.dot(r);
		
		cout<< n << " " << rr << endl;	
		if(rr < epsilon_rr){
			cout<< "The residual |Ax-b| is sufficiently small. " <<endl;
			break;
		}	
		
		b = rr / r0r0;
	       	p = r + b*p0;	
		// Update variables. 
		x0 = x; p0 = p; r0 = r;	
	}
	//cout << "\nx = \n" << x << endl;
	cout << "\n|x-x_0| = "<<norm_x() << endl;
	return 0;
}
