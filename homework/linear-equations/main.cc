#include"matrix.h"
#include <cstdlib>
#include <ctime>
#include <iostream>

pp::matrix eye(int n){
	pp::matrix M(n,n);
	for(int i=0;i<n;i++)M[i,i]=1;
	return M;
}

// helper to create an n-by-m matrix with random values in [0,1)
pp::matrix randmat(int n, int m){
	pp::matrix A(n,m);
	for(int j=0;j<m;j++){
		for(int i=0;i<n;i++){
			A[i,j] = (double)std::rand()/RAND_MAX;
		}
	}
	return A;
}

// helper to create random vector of length n
pp::vector randvec(int n){
	pp::vector v(n);
	for(int i=0;i<n;i++) v[i] = (double)std::rand()/RAND_MAX;
	return v;
}

// check approximate equality between matrices
bool mat_approx(const pp::matrix& A, const pp::matrix& B, double tol=1e-6){
	if(A.size1()!=B.size1() || A.size2()!=B.size2()) return false;
	for(int i=0;i<A.size1();i++) for(int j=0;j<A.size2();j++)
		if(!pp::approx(A[i,j], B[i,j], tol, tol)) return false;
	return true;
}

int main(int argc, char** argv){
	std::srand(std::time(nullptr));

	// check command-line for timing mode
	int timingSize = 0;
	for(int i=1;i<argc;i++){
		if(std::sscanf(argv[i],"-size:%d",&timingSize)==1) break;
	}
	if(timingSize > 0){
		// perform one QR decomposition on an NxN matrix and exit
		pp::matrix A = randmat(timingSize, timingSize);
		pp::QR qr = pp::QR::decomp(A);
		(void)qr; // suppress unused warning
		return 0;
	}

	// --- QR decomposition tests ------------------------------------------------
	{
		int n1 = 6, m1 = 4; // tall matrix
		pp::matrix A = randmat(n1, m1);
		pp::QR qr = pp::QR::decomp(A);

		std::cout << "A:\n"; A.print();
		std::cout << "Q:\n"; qr.Q.print();
		std::cout << "R:\n"; qr.R.print();

		// check upper triangular R
		bool upper = true;
		for(int i=0;i<m1;i++) for(int j=0;j<i;j++)
			if(std::abs(qr.R[i,j]) > 1e-8) upper = false;
		std::cout << "R upper triangular? " << (upper?"yes":"no") << std::endl;

		// check Q^T Q = I
		pp::matrix QtQ = qr.Q.T() * qr.Q;
		std::cout << "Q^T Q:\n"; QtQ.print();
		std::cout << "Q^T Q approx I? " << (mat_approx(QtQ, eye(m1))?"yes":"no") << std::endl;

		// check QR = A
		pp::matrix QRprod = qr.Q * qr.R;
		std::cout << "A:\n"; A.print();
		std::cout << "Q*R:\n"; QRprod.print();
		std::cout << "Q*R approx A? " << (mat_approx(QRprod, A)?"yes":"no") << std::endl;
	}

	// - -- solve test ------------------------------------------------------------
	{
		int n2 = 5;
		pp::matrix A = randmat(n2, n2);
		pp::vector b = randvec(n2);
		pp::QR qr = pp::QR::decomp(A);
		pp::vector x = qr.solve(b);
		pp::vector Ax = A * x;
		std::cout << "Ax:\n"; Ax.print();
		std::cout << "b:\n"; b.print();
		std::cout << "Ax ~= b? " << (pp::approx(Ax, b)?"yes":"no") << std::endl;

		// determinant and inverse (used later)
		double detA = qr.det();
		std::cout << "det(A) = " << detA << std::endl;
	}

	// --- inverse-only test ----------------------------------------------------
	{
		int n3 = 5;
		pp::matrix A = randmat(n3, n3);
		pp::QR qr = pp::QR::decomp(A);
		pp::matrix B = qr.inverse();
		pp::matrix AB = A * B;
		std::cout << "AB:\n"; AB.print();
		std::cout << "AB approx I? " << (mat_approx(AB, eye(n3))?"yes":"no") << std::endl;
	}
	return 0;
}
