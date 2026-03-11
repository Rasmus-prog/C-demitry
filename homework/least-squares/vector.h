#pragma once
#include<string>
#include<vector>
#include<initializer_list>
#include<functional>

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

	vector& add(double x);
	vector& push_back(double x);

	double norm() const;
	void print(std::string s="") const;
	vector map(std::function<double(double)> f) const;
};

vector operator+(const vector& a, const vector& b);
vector operator-(const vector& a);
vector operator-(const vector& a, const vector& b);
vector operator*(const vector& a, double c);
vector operator*(double c, const vector& a);
vector operator/(const vector& a, double c);

bool approx(double x, double y, double acc=1e-6, double eps=1e-6);
bool approx(const vector& a, const vector& b, double acc=1e-6, double eps=1e-6);
}
