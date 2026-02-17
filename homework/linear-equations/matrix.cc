#include"matrix.h"
#include<string>
#include<algorithm>
#include<cmath>
#include<iostream>
#include<cassert>
#include<stdio.h>
#define SELF (*this)
#define FORV(i,v) for(int i=0;i<v.size();i++)
#define FOR_COLS(i,A) for(int i=0;i<A.size2();i++)
namespace pp{

// utility comparison functions
bool approx(double x, double y, double acc, double eps){
	if(std::abs(x-y) < acc) return true;
	if(std::abs(x-y) < eps*std::max(std::abs(x),std::abs(y))) return true;
	return false;
}

bool approx(const vector& a, const vector& b, double acc, double eps){
	if(a.size() != b.size()) return false;
	for(int i=0;i<a.size();i++)
		if(!approx(a[i],b[i],acc,eps))return false;
	return true;
}

// implementation details follow


vector& vector::operator+=(const vector& other) {
	FORV(i,SELF) data[i]+=other.data[i];
	return SELF; }

vector& vector::operator-=(const vector& other) {
	FORV(i,SELF) data[i]-=other.data[i];
	return SELF; }

vector& vector::operator*=(double x) {
	FORV(i,SELF) data[i]*=x;
	return SELF; }

vector& vector::operator/=(double x) {
	FORV(i,SELF) data[i]/=x;
	return SELF; }

vector& vector::add(double x){
	data.push_back(x);
	return SELF;}

vector& vector::push_back(double x){
	data.push_back(x);
	return SELF;}

double vector::norm() const {
	double s=0;
	FORV(i,SELF)s+=SELF[i]*SELF[i];
	return std::sqrt(s);
	}

vector vector::map(std::function<double(double)> f) const{
	vector r=SELF;
	for(int i=0;i<r.size();i++)r[i]=f(r[i]);
	return r;
	}

void vector::print(std::string s) const {
	std::cout << s;
	FORV(i,SELF)printf("%9.3g ",(double)SELF[i]);
	printf("\n");
	}

vector operator/(const vector& v, double x){
	vector r=v;
	r/=x;
	return r; }

vector operator*(const vector& v, double x){
	vector r=v;
	r*=x;
	return r; }

vector operator*(double x,const vector& a){ return a*x; }

vector operator+(const vector& a, const vector& b){
	vector r=a;
	r+=b;
	return r; }

vector operator-(const vector& a){
	vector r=a;
	for(int i=0;i<r.size();i++)r[i]=-r[i];
	return r; }

vector operator-(const vector& a, const vector& b){
	vector r=a;
	r-=b;
	return r; }

void matrix::resize(int n, int m){
	cols.resize(m);
	for(int i=0;i<m;++i)cols[i].resize(n);
	}

matrix& matrix::operator+=(const matrix& other) {
	FOR_COLS(i,SELF) SELF[i]+=other[i];
	return SELF; }

matrix& matrix::operator-=(const matrix& other) {
	FOR_COLS(i,SELF) SELF[i]-=other[i];
	return SELF; }

matrix& matrix::operator*=(double x) {
	FOR_COLS(i,SELF) SELF[i]*=x;
	return SELF; }

matrix& matrix::operator/=(double x) {
	FOR_COLS(i,SELF) SELF[i]/=x;
	return SELF; }

matrix operator/(const matrix& A,double x){
	matrix R=A;
	R/=x;
	return R; }

matrix operator*(const matrix& A,double x){
	matrix R=A;
	R*=x;
	return R; }

matrix operator*(double x,const matrix& A){
	return A*x; }

matrix operator+(const matrix& A, const matrix& B){
	matrix R=A;
	R+=B;
	return R; }

matrix operator-(const matrix& A, const matrix& B){
	matrix R=A;
	R-=B;
	return R; }

vector operator*(const matrix& M, const vector& v){
	vector r; r.resize(M.size1());
	for(int i=0;i<r.size();i++){
		double sum=0;
		for(int j=0;j<v.size();j++)sum+=M[i,j]*v[j];
		r[i]=sum;
		}
	return r;
	}

matrix operator*(const matrix& A, const matrix& B){
	matrix R(A.size1(),B.size2());
	for(int k=0;k<A.size2();k++)
	for(int j=0;j<B.size2();j++)
		{
		for(int i=0;i<A.size1();i++)R[i,j]+=A[i,k]*B[k,j];
		}
	return R;
	}

void matrix::setid(){
	assert(size1()==size2());
	for(int i=0;i<size1();i++){
		SELF[i,i]=1;
		for(int j=i+1;j<size1();j++)SELF[i,j]=SELF[j,i]=0;
		}
	}

matrix matrix::transpose() const {
	matrix R(size2(),size1());
	for(int i=0;i<R.size1();i++)
		for(int j=0;j<R.size2();j++) R[i,j]=SELF[j,i];
	return R;
	}

matrix matrix::T() const {return SELF.transpose();}

void matrix::print(std::string s) const {
	std::cout << s << std::endl;
	for(int i=0;i<size1();i++){
		for(int j=0;j<size2();j++)printf("%9.3g ",(double)SELF[i,j]);
		printf("\n");
		}
	printf("\n");
	}
// --- QR decomposition implementation ------------------------------------------------

// perform modified Gram-Schmidt on A (n x m, n>=m). Q will have orthonormal cols.
QR::QR(const matrix& A) {
    int n = A.size1();
    int m = A.size2();
    Q = A;
    // prepare R as m x m zero matrix
    R.resize(m, m);
    for(int i=0;i<m;i++) for(int j=0;j<m;j++) R[i,j]=0;

    for(int k=0;k<m;k++){
        // compute norm of k-th column of Q
        double norm = 0;
        for(int i=0;i<n;i++) norm += Q[i,k]*Q[i,k];
        norm = std::sqrt(norm);
        R[k,k] = norm;
        if(norm <= 0) continue; // zero column, leave Q column as is
        // normalize
        for(int i=0;i<n;i++) Q[i,k] /= norm;
        // orthogonalize remaining columns
        for(int j=k+1;j<m;j++){
            double dot = 0;
            for(int i=0;i<n;i++) dot += Q[i,k] * Q[i,j];
            R[k,j] = dot;
            for(int i=0;i<n;i++) Q[i,j] -= dot * Q[i,k];
        }
    }
}

vector QR::solve(const vector& b) const {
    int n = Q.size1();
    int m = R.size1();
    assert((int)b.size() == n);
    vector y(m);
    // y = Q^T * b
    for(int i=0;i<m;i++){
        double sum = 0;
        for(int k=0;k<n;k++) sum += Q[k,i] * b[k];
        y[i] = sum;
    }
    // back substitution R x = y
    vector x(m);
    for(int i=m-1;i>=0;i--){
        double sum = y[i];
        for(int j=i+1;j<m;j++) sum -= R[i,j] * x[j];
        x[i] = sum / R[i,i];
    }
    return x;
}

double QR::det() const {
    int m = R.size1();
    assert(m == R.size2());
    double d = 1;
    for(int i=0;i<m;i++) d *= R[i,i];
    return d;
}

matrix QR::inverse() const {
    int m = R.size1();
    assert(m == R.size2());
    matrix inv(m,m);
    // construct identity columns
    for(int j=0;j<m;j++){
        vector e(m);
        for(int i=0;i<m;i++) e[i] = (i==j ? 1.0 : 0.0);
        vector col = solve(e);
        for(int i=0;i<m;i++) inv[i,j] = col[i];
    }
    return inv;
}
}
