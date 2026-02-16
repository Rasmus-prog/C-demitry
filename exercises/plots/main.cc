#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>

double erf(double x){
// single precision error function (Abramowitz and Stegun, from Wikipedia)
if(x<0) return -erf(-x);
std::vector<double> a {0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
double t=1/(1+0.3275911*x);
double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
return 1-sum*std::exp(-x*x);
} 

int main() {
    // tabulated values for erf(x) from Abramowitz & Stegun (approx)
    std::vector<double> xs {0.0, 0.5, 1.0, 1.5, 2.0};
    std::vector<double> tab {0.0, 0.5204998778, 0.8427007929, 0.9661051465, 0.9953222650};

    // write tabulated points to file
    std::ofstream tabfile("tab.dat");
    for (size_t i = 0; i < xs.size(); ++i) {
        tabfile << xs[i] << " " << tab[i] << "\n";
    }
    tabfile.close();

    // write smooth curve data
    std::ofstream curvefile("curve.dat");
    for (double x = 0.0; x <= 2.0; x += 0.01) {
        curvefile << x << " " << erf(x) << "\n";
    }
    curvefile.close();

    // invoke gnuplot to display the plot
    // the -persist flag keeps the window open after gnuplot exits
    // produce both an interactive window and a PDF file
    std::string gp =
        "gnuplot -persist -e \"set terminal pdfcairo; set output 'erf.pdf'; \
        plot 'curve.dat' with lines title 'erf', \
             'tab.dat' with points pt 7 ps 1.5 title 'table'; \
        set terminal qt; set output;\"";
    std::system(gp.c_str());

    return 0;
}
