#pragma once
#include"matrix.h"

namespace pp{
struct QR {
    matrix Q;
    matrix R;

    QR() = default;
    explicit QR(const matrix& A);
    static QR decomp(const matrix& A) { return QR(A); }

    vector solve(const vector& b) const;
    double det() const;
    matrix inverse() const;
};
}
