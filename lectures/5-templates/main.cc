#include<iostream>
#include<string>
#include<vector>
#include<functional>

struct func{
    double operator()(double x) {
        return x * x * x;
    }
};


double square(double x) {
    return x * x;
}

double eval(std::function<double(double)> f, double x=2) {
    return f(x);
}


template<typename T>
void echo(const T& arg) {
    std::cout << "echo: arg= " << arg << std::endl;
}



template<typename T>
struct list{
    T* data;
    int size;
    T& operator[](int index) { return data[index]; }
    list(int n) : size(n) {
        data = new T[n];
    }
    ~list() {size = 0; delete[] data;}

};




template<typename T>
struct bar{
    T datum;
    void print(){
        std::cout << "bar: datum= " << datum << std::endl;
    }
};


int main() {
    echo(1);
    echo(1.234);
    echo("Hello");
    bar<int> intbar {1};
    bar<double> doublebar {1.234};
    bar<std::string> stringbar {"Hello"};
    intbar.print();
    doublebar.print();
    stringbar.print();
    int n = 5;
    list<double> doublelist(n);
    doublelist[0] = 1.23;


    std::cout << "square(2) = "<< square(2) << std::endl;
    std::function<double(double)> f = square;
    std::cout << "f(2) = "<< f(2) << std::endl;

    std::cout << "eval(f,2) = "<< eval(f,2) << std::endl;
    std::cout << "eval(square,2) = "<< eval(square,2) << std::endl;
    func cube;
    std::function<double(double)> g = cube;
    std::cout << "cube(2) = "<< cube(2) << std::endl;
    std::cout << "g(2) = "<< g(2) << std::endl;

    double a = 1;
    auto h = [a]() -> double { return a;};
    auto q = [&a]() -> double { return a;};
    std::cout << "h() = "<< h() << std::endl;
    std::cout << "q() = "<< q() << std::endl;
    a = 2;
    std::cout << "h() = "<< h() << std::endl;
    std::cout << "q() = "<< q() << std::endl;
    return 0;
}