#include"vector.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<stdio.h>

#define SELF (*this)
#define FORV(i,v) for(int i=0;i<v.size();i++)

namespace pp{

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

}
