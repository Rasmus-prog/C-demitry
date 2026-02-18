#include <iostream>
#include <thread>
#include <string>
#include <vector>

//Calculate the harmonic sum, sumi=a…b 1/i, sharing the computation load between several processors. 
struct datum {
    int start;
    int end;
    double sum;
};

void harm(datum& p) {
    int start = p.start;
    int end = p.end;
    double sum = 0.0;
    for (int i = start; i <= end; ++i) {
        sum += 1.0 / i;
    }
    p.sum = sum;
}





int main(int argc ,char** argv) {
    int nterm = (int)1e9;
    int nthreads = 1;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-terms" && i + 1 < argc) {
            nterm = (int)std::stod(argv[++i]);
        }
        if (arg == "-threads" && i + 1 < argc) {
            nthreads = std::stoi(argv[++i]);
        }
    }
    std::cout << "terms: " << nterm << ", threads: " << nthreads << std::endl;

    std::vector < std::thread > threads;
    threads.reserve(nthreads);
    std::vector < datum > data(nthreads);

    int chunk = nterm / nthreads;
    for(int i = 0; i < nthreads; ++i) {
        data[i].start = i * (chunk) + 1;
        data[i].end = (i + 1) * (chunk);
        threads.emplace_back(harm, std::ref(data[i]));
    }
    for(std::thread& t : threads) {
        t.join();
    }
    double total_sum = 0.0;
    for (const datum& d : data) {
        total_sum += d.sum;        
    }
std::cout << "Harmonic sum: " << total_sum << std::endl;

return 0;
}