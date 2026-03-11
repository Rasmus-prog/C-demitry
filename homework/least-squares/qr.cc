#include"qr.h"
#include<cmath>
#include<cassert>

namespace pp{

QR::QR(const matrix& A) {
    int n = A.size1();
    int m = A.size2();
    Q = A;
    R.resize(m, m);
    for(int i=0;i<m;i++) for(int j=0;j<m;j++) R[i,j]=0;

    for(int k=0;k<m;k++){
        double norm = 0;
        for(int i=0;i<n;i++) norm += Q[i,k]*Q[i,k];
        norm = std::sqrt(norm);
        R[k,k] = norm;
        if(norm <= 0) continue;
        for(int i=0;i<n;i++) Q[i,k] /= norm;
        for(int j=k+1;j<m;j++){
            double dot = 0;
            for(int i=0;i<n;i++) dot += Q[i,k] * Q[i,j];
            R[k,j] = dot;
            for(int i=0;i<n;i++) Q[i,j] -= dot * Q[i,k];
        }
    }
}

vector QR::solve(const vector& b) const {
    int n = Q.size1();
    int m = R.size1();
    assert((int)b.size() == n);
    vector y(m);
    for(int i=0;i<m;i++){
        double sum = 0;
        for(int k=0;k<n;k++) sum += Q[k,i] * b[k];
        y[i] = sum;
    }
    vector x(m);
    for(int i=m-1;i>=0;i--){
        double sum = y[i];
        for(int j=i+1;j<m;j++) sum -= R[i,j] * x[j];
        x[i] = sum / R[i,i];
    }
    return x;
}

double QR::det() const {
    int m = R.size1();
    assert(m == R.size2());
    double d = 1;
    for(int i=0;i<m;i++) d *= R[i,i];
    return d;
}

matrix QR::inverse() const {
    int m = R.size1();
    assert(m == R.size2());
    matrix inv(m,m);
    for(int j=0;j<m;j++){
        vector e(m);
        for(int i=0;i<m;i++) e[i] = (i==j ? 1.0 : 0.0);
        vector col = solve(e);
        for(int i=0;i<m;i++) inv[i,j] = col[i];
    }
    return inv;
}

}
