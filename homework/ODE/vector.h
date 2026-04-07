#pragma once
#include <initializer_list>
#include <vector>

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

	double norm() const;
};

vector operator+(const vector& a, const vector& b);
vector operator-(const vector& a);
vector operator-(const vector& a, const vector& b);
vector operator*(const vector& a, double c);
vector operator*(double c, const vector& a);
vector operator/(const vector& a, double c);
}
