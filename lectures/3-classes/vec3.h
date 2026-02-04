#pragma once
#include <iostream>
#include <string>
namespace pp {
    struct vec3 {
        double x,y,z;
        //constructors
        vec3(double a, double b, double c){ //parmeterized constructor
            std::cout<<"vec3 constructor called"<<std::endl;
            x=a; y=b; z=c;
        }
        vec3() : vec3(0,0,0) //defaul constructor
        {
            std::cout<<"vec3 default constructor called"<<std::endl;
        }
        vec3(const vec3&)=default; //copy constructor
        vec3(vec3&&)=default; //move constructor

        //destructor
        ~vec3(){
            std::cout<<"vec3 destructor called"<<std::endl;
        }
        //assignment operators
        vec3& operator=(const vec3&); //copy assignment
        vec3& operator=(vec3&&); //move assignment
        
        //member operations
        vec3& operator+=(double);
        vec3& operator-=(double);
        vec3& operator*=(double);
        vec3& operator/=(double);

        void print(const std::string& s="");

        // stream output
        friend std::ostream& operator<<(std::ostream&, const vec3&);
    };
// non-member operations
    vec3 operator-(const vec3&);
    vec3 operator+(const vec3&, const vec3&);
    vec3 operator-(const vec3&, const vec3&);
    vec3 operator*(const vec3&, double);
    vec3 operator*(double, const vec3&);
    vec3 operator/(const vec3&, double);

}