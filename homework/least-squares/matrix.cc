#include"matrix.h"
#include<cassert>
#include<iostream>
#include<stdio.h>
#define SELF (*this)
#define FOR_COLS(i,A) for(int i=0;i<A.size2();i++)
namespace pp{



// check approximate equality between matrices
bool mat_approx(const pp::matrix& A, const pp::matrix& B, double tol){
	if(A.size1()!=B.size1() || A.size2()!=B.size2()) return false;
	for(int i=0;i<A.size1();i++) for(int j=0;j<A.size2();j++)
		if(!pp::approx(A[i,j], B[i,j], tol, tol)) return false;
	return true;
}


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
}
