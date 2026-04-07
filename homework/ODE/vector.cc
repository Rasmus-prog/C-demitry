#include"vector.h"
#include<cmath>

#define SELF (*this)
#define FORV(i,v) for(int i=0;i<v.size();i++)

namespace pp{

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

double vector::norm() const {
	double s=0;
	FORV(i,SELF)s+=SELF[i]*SELF[i];
	return std::sqrt(s);
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
