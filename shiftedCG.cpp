#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <vector>  
using namespace std;
using namespace Eigen;

int L = 100;
int n_iter = 1000;
int n_s_max = 5;	
double epsilon_rr = 1e-30;
/**** Ax = y ****/
MatrixXd A = MatrixXd::Zero(L,L);
VectorXd y(L),x(L),p(L),r(L),q(L);
double a0,a,b0,b,rh0,rh;

MatrixXd x0_s,x_s(L,n_s_max),p_s(L,n_s_max); //,r_s(L,n_s_max);
VectorXd a_s(n_s_max),b_s(n_s_max);
VectorXd c0_s(n_s_max),c1_s(n_s_max),c2_s(n_s_max); //,rh0_s(n_s_max),rh_s(n_s_max);

double s;	// Shift variable.

void set_A_y(){
	y = VectorXd::Ones(L);
	for(int i=0;i<L;++i){
	for(int j=i;j<L;++j){
		A(i,j) = j+1;	
		A(j,i) = j+1;	
	}}
}

void init(){
	x= VectorXd::Zero(L);
	p = VectorXd::Zero(L);
	r = y;
	a0 = 1.0;
	rh0 = 1.0;
	b0 = 0.0;
}

void init_shifted(){
	x_s= MatrixXd::Zero(L,n_s_max);
	p_s= MatrixXd::Zero(L,n_s_max);
	
	c0_s = VectorXd::Ones(n_s_max); 
	c1_s = VectorXd::Ones(n_s_max); 
}

int main(int argc, char const* argv[]){
	set_A_y();
	cout << "A = \n" << A << endl;
	cout << "y = \n" << y << endl;
	init();
	init_shifted();
	double rr,r0r0; 
	VectorXd Ar0;

	cout<<endl;
	for(int k=0;k<n_iter;++k){
		//Seed System.
		rh = r.dot(r); 
		b =  - rh / rh0; // Init; rh0 = 1.0
		p = r - b*p;
		q = A * p;
		a = rh / r.dot(q);
		x = x + a*p;
		
		for(int j=0; j<n_s_max; ++j){
			s = 1 * j+1;
			// Shifted equations.
			c2_s[j] = (1.0+a*s)*c1_s[j] + (a*b/a0)*(c0_s[j]-c1_s[j]);
			b_s[j] = pow((c0_s[j]/c1_s[j]),2)*b;
			a_s[j] = (c1_s[j]/c2_s[j])*a;
			p_s.col(j) = (1.0/c1_s[j])*r - b_s[j]*p_s.col(j);
			x_s.col(j) = x_s.col(j) + a_s[j]*p_s.col(j);
			//Update.
			c0_s[j] = c1_s[j]; c1_s[j] = c2_s[j]; 
				
		}
	
		cout<<k<< " " << rh << endl;
		if(abs(rh)<epsilon_rr){
			cout<< "The residual |Ax-b| is sufficiently small. " <<endl;
			break;
		}
		
		//Update.
		r = r - a*q;
		rh0 = rh; a0 = a; 	
	}
	
	cout << "x = \n" << x << endl;
	for(int j=0; j<n_s_max; ++j){
		cout << "x" << j << " = \n" << x_s.col(j) << endl;
	}
	return 0;
}
