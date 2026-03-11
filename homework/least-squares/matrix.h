#pragma once
#include<string>
#include<vector>
#include"vector.h"
namespace pp{
struct matrix {
	std::vector<pp::vector> cols;
	matrix()=default;
	matrix(int n,int m) : cols(m, pp::vector(n)) {}
	matrix(const matrix& other)=default;
	matrix(matrix&& other)=default;
	matrix& operator=(const matrix& other)=default;
	matrix& operator=(matrix&& other)=default;
	int size1() const {return cols.empty() ? 0 : cols[0].size(); }
	int size2() const {return cols.size();}
	double& operator()(int i, int j){return cols[j][i];}
	double& operator[](int i, int j){return cols[j][i];}
	const double& operator()(int i, int j)const{return cols[j][i];}
	const double& operator[](int i, int j)const{return cols[j][i];}
	vector& operator[](int i){return cols[i];}
	const vector& operator[](int i) const {return cols[i];}
	void resize(int n, int m);
	void setid();
	matrix transpose() const;
	matrix T() const;
	
	double get (int i, int j) {return cols[j][i];}
	void set(int i, int j, double value){cols[j][i] = value;}

	matrix& operator+=(const matrix&);
	matrix& operator-=(const matrix&);
	matrix& operator*=(const matrix&);
	matrix& operator*=(const double);
	matrix& operator/=(const double);
	matrix  operator^(int);

	void print(std::string s="") const;
};

matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const matrix&, double);
matrix operator*(double, const matrix&);
matrix operator/(const matrix&, double);
vector operator*(const matrix&, const vector&);
bool mat_approx(const matrix& A, const matrix& B, double tol=1e-6);
}
