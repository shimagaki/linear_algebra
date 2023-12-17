#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <vector>  
using namespace std;
using namespace Eigen;

int L = 5;
int n_iter = 10;
int n_s_max = 5;	
double epsilon_rr = 1e-10;
/**** Ax = y ****/
MatrixXd A = MatrixXd::Zero(L,L);
VectorXd y,x0,x,p0,p,r0,r;
double a0,a,b0,b;

MatrixXd x0_s,p0_s,x_s(L,n_s_max),p_s(L,n_s_max),r_s(L,n_s_max);
VectorXd a_s(n_s_max),b_s(n_s_max);
VectorXd c0_s(n_s_max),c1_s(n_s_max),c2_s(n_s_max);


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
	x0= VectorXd::Zero(L);
	p0 = y;
	r0 = y;
	a0 = 1.0;
	b0 = 0.0;
}

void init_shifted(){
	x0_s= MatrixXd::Zero(L,n_s_max);
	p0_s= MatrixXd::Zero(L,n_s_max);
	
	for(int n_s=0;n_s<n_s_max;++n_s){
		p0_s.col(n_s) = y;
		//r_s.col(n_s) = y;
		r_s.col(n_s) = VectorXd::Zero(n_s_max);
	}
	
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
	for(int n=0;n<n_iter;++n){
		r0r0 = r0.dot(r0); 
		a = r0r0 / p0.dot(A*p0);
		x = x0 + a * p0;
		r = r0 - a * A * p0; 
		rr = r.dot(r); 
		b = rr / r0r0;
		p = r + b * p0; 
		cout<<n<<" ";
		for(int n_s=0; n_s<n_s_max; ++n_s){
			s = 1 * n_s;
			///*	
			// Shifted equations.
			c2_s[n_s] = c1_s[n_s]*c0_s[n_s]*a0 / ( c0_s[n_s]*a0*( 1 + a*s ) + a*b0*(c0_s[n_s] - c1_s[n_s]) );
			a_s[n_s] = c2_s[n_s]/c1_s[n_s] * a;
			x_s.col(n_s) = x0_s.col(n_s) + a_s[n_s] * p0_s.col(n_s);
			b_s[n_s] = pow((c2_s[n_s]/c1_s[n_s]),2) * b;
			r_s.col(n_s) = c2_s[n_s] * r;
			p_s.col(n_s) = r_s.col(n_s) + b_s[n_s] * p0_s.col(n_s);

			cout << " " << r_s.col(n_s).dot(r_s.col(n_s)) << endl;
			// Update of variables..
			x0_s.col(n_s)=x_s.col(n_s); 
			p0_s.col(n_s)=p_s.col(n_s); 
			c0_s[n_s] = c1_s[n_s]; c1_s[n_s] = c2_s[n_s]; 
			//*/	
		}cout<<endl;
			
		// Update of variables..
		a0=a; b0=b; x0=x; p0=p; r0=r;
			
		rr = r.dot(r);
		cout<< n << " " << rr << endl;	
		if(abs(rr)<epsilon_rr){
			cout<< "The residual |Ax-b| is sufficiently small. " <<endl;
			break;
		}


	}
	cout << "x = \n" << x << endl;
	for(int n_s=0; n_s<n_s_max; ++n_s){
		cout << "x" << n_s << " = \n" << x_s.col(n_s) << endl;
	}
	return 0;
}
