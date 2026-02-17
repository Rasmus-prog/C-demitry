#pragma once
#include<string>
#include<vector>
#include<initializer_list>
#include<functional>
#include<cmath>
#include<iostream>
namespace pp{
struct vector {
	std::vector<double> data;
	vector(int n) : data(n) {}
	vector(std::initializer_list<double> list) :
		data(list.begin(),list.end()) {}
	vector()			=default;
	vector(const vector&)		=default;
	vector(vector&&)		=default;
	vector& operator=(const vector&)=default;
	vector& operator=(vector&&)	=default;
	int size() const {return data.size();}
	void resize(int n) {data.resize(n);}
	double& operator[](int i) {return data[i];}
	const double& operator[](int i) const {return data[i];}

	vector& operator+=(const vector& other);

	vector& operator-=(const vector& other);

	vector& operator*=(double c);

	vector& operator/=(double c);

	// mutators also available in implementation
	vector& add(double x);
	vector& push_back(double x);

	double norm() const;

	void print(std::string s="") const;

	vector map(std::function<double(double)> f) const;

}; //vector

// forward declaration of matrix so QR can refer to it

// global vector operators - defined in matrix.cc
vector operator+(const vector& a, const vector& b);
vector operator-(const vector& a);
vector operator-(const vector& a, const vector& b);
vector operator*(const vector& a, double c);
vector operator*(double c, const vector& a);
vector operator/(const vector& a, double c);

bool approx(double x, double y, double acc=1e-6, double eps=1e-6);

bool approx(const vector& a, const vector& b, double acc=1e-6, double eps=1e-6);

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
//	void resize(int n, int m);
	void setid();
	matrix transpose() const;
	matrix T() const;
	
	double get (int i, int j) {return cols[j][i];}
	void set(int i, int j, double value){cols[j][i] = value;}
//	vector get_col(int j);
//	void set_col(int j,vector& cj);

	matrix& operator+=(const matrix&);
	matrix& operator-=(const matrix&);
	matrix& operator*=(const matrix&);
	matrix& operator*=(const double);
	matrix& operator/=(const double);
	matrix  operator^(int);

	void print(std::string s="") const;
};

// QR decomposition helper class using modified Gram-Schmidt orthogonalization.
// Provides methods to decompose a (tall) matrix A = Q*R, solve linear systems,
// compute the determinant, and obtain the inverse for square matrices.
struct QR {
    matrix Q; // orthonormal columns (n x m)
    matrix R; // upper triangular (m x m)

    QR() = default;
    explicit QR(const matrix& A);          // perform decomposition
    static QR decomp(const matrix& A) { return QR(A); }

    // Solve QR x = b (least-squares when n>=m). b must have size n.
    vector solve(const vector& b) const;

    // Determinant of the original matrix (must be square).
    double det() const;

    // Inverse of the original matrix (must be square and non-singular).
    matrix inverse() const;
};

matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const matrix&, double);
matrix operator*(double, const matrix&);
matrix operator/(const matrix&, double);
vector operator*(const matrix&, const vector&);
}
