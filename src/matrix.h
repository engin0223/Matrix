#pragma once
#include<numeric>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <string>
#include <sstream> 
#include <string>
#include <windows.h>
#include <random>


// Class for matrix.
// To create a matrix choose the type of entries->T, create a vector of vector(s)->myVec
// then write matrix<T> matrix_name(myVec).
//(E.g "matrix<int> A({{1,2,3}})" creates a 1x3 matrix(i.e row vector) which is [1,2,3]))
template <typename T>
class matrix {
public:
	std::vector<std::vector<T>> X;
	int m;
	int n;
	matrix<T>() {};
	matrix<T>(std::vector<std::vector<T>> B) {
		m = int(B.size());
		n = int(B[0].size());

		X = B;
	}

	T get(int i, int j) {
		return X[i - 1][j - 1];
	}

	void rep(int i, int j, T n) {
		X[i - 1][j - 1] = n;
	}

	void change(std::vector<std::vector<T>> N) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				rep(i, j, N.at(i, j));
			}
		}
	}
};


template <class C>
#define name(C) #C  // gets the variable name of class C object.


void s() {}


matrix<double> matrix_rng(int i, int j) {
	std::vector<std::vector<double>> B;
	std::vector<double> row(j, 0);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distr(-100, 100);
	for (int k = 0; k < i; ++k) {
		for (int l = 0; l < j; ++l) {
			row.at(l) = distr(gen);
		}
		B.push_back(row);
	}
	return matrix<double>(B);
}

// Creates a ixj matrix whose entries are "type(T) varible_name(n)".(E.g. "double number n = 5.0").
template <typename T>
matrix<T> matrix_n(int i, int j, T n) {
	std::vector<std::vector<T>> B;
	std::vector<T> row(j, (T) n);
	for (int k = 0; k < i; ++k) {
		B.push_back(row);
	}
	return matrix<T>(B);
}

// returns transpose of matrix.
template <typename K>
matrix<K> T(matrix<K> B) {
	matrix<K> G = B;
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.m; ++j) {
			G.rep(i, j, B.get(j, i));
		}
	}
	return G;
}

// returns trace of matrix.
template <typename T>
T Tr(matrix<T> B) {
	T sum{};
	for (int i = 1; i <= B.m; ++i) {
		sum += B.get(i, i);
	}
	return sum;
}


// Creates a NxN matrix Identity matrix.(i.e matrix[i][i] = 1 for all i that is 1<=i<=N and all other entries are 0)
template <typename T>
matrix<T> I(T N) {
	std::vector<std::vector<T>> B;
	std::vector<T> row(N, 0);
	for (int k = 0; k < N; ++k) {
		B.push_back(row);
	}
	for (int i = 0; i < N; ++i) {
		B[i][i] = 1;
	}
	return matrix<T>(B);
}


// returns the (m,n) minor of the matrix.
template <typename T>
matrix<T> M(int m, int n, matrix<T> B) {
	std::vector<std::vector<T>> G = B.X;
	int k = 0;
	for (int i = 1; i <= B.m; ++i) {
		if (i != m) {
			G[i - 1 - k].erase(G[i - 1 - k].begin() + n - 1);
		}
		else {
			G.erase(G.begin() + i - 1);
			k = 1;
		}
	}
	matrix<T> C(G);
	return C;
}


// return determinant of matrix.
template <typename T>
T det(matrix<T> B) {
	T sum{};
	if (B.n == 2) {
		return B.get(1, 1) * B.get(2, 2) - B.get(1, 2) * B.get(2, 1);
	}
	else if (B.n == 1) {
		return B.get(1, 1);
	}
	for (int j = 1; j <= B.m; ++j) {
		sum += (-2 * ((1 + j) % 2) + 1) * B.get(1, j) * det(M(1, j, B));
	}
	return sum;
}


// returns adjugate of matrix.
template <typename T>
matrix<T> adj(matrix<T> B) {
	matrix<T> G = B;
	if (B.m == 2) {
		G.rep(1, 1, B.get(2, 2));
		G.rep(1, 2, -B.get(1, 2));
		G.rep(2, 1, -B.get(2, 1));
		G.rep(2, 2, B.get(1, 1));
		return G;
	}
	else if (B.m == 1) {
		G.rep(1, 1, B.get(1, 1));
		return G;
	}
	else {
		for (int i = 1; i <= B.m; ++i) {
			for (int j = 1; j <= B.n; ++j) {
				G.rep(i, j, (-2 * ((i + j) % 2) + 1) * det(M(j, i, B)));
			}
		}
		return G;
	}
}


// returns invert of matrix.
template <typename T>
matrix<T> inv(matrix<T> B) {
	try {
		if (det(B) == 0) {
			throw (0);
		}
		else {
			return adj(B) * (T)(pow(det(B), -1.0));
		}
	}
	catch (int n) {
		std::cout << "determinant of the matrix is zero therefore the matrix is singular(not-invertible)";
		return matrix<T>();
	}
}


// Takes LU decomposition of matrix and returns a vector with entries are L and U.(i.e {L, U})
template <typename T>
std::vector<matrix<T>> LU(matrix<T> B) {
	int N = B.m;
	matrix<T> A_n(B.X);
	matrix<T> L(I<double>(N).X);
	matrix<T> L_n(I<double>(N).X);
	for (int n = 1; n <= N - 1; ++n) {
		L_n = I<double>(N);
		for (int i = n + 1; i <= N; ++i) {
			L_n.rep(i, n, -(A_n.get(i, n)) / (A_n.get(n, n)));
		}
		for (int i = n + 1; i <= N; ++i) {
			for (int j = 1; j <= N; ++j) {
				A_n.rep(i, j, A_n.get(i, j) + A_n.get(n, j) * L_n.get(i, n));
			}
		}
		L = L * inv(L_n);
	}
	matrix<T> U((inv(L) * B).X);
	return { L, U };

}


template <typename T>
double eigenfunc(double x, matrix<T> B) {
	double v = det(B - x * I<T>(B.m));
	return det(B - x * I<T>(B.m));
}


template <typename T>
double detderive(double x, matrix<T> B) {
	return Tr<T>((T)-1 * adj(B - x * I<T>(B.m)));
}

template <typename T>
double product_p(double x, std::vector<T> A) {
	double prod = 0;
	for (int i = 0; i < A.size(); ++i) {
		double val = 1;
		for (int j = 0; j < A.size(); ++j) {
			if (i != j) {
				val *= x - A.at(j);
			}
		}
		prod += val;
	}
	return prod;
}

template <typename T>
double product(double x, std::vector<T> A) {
	double prod = 1;
	for (int i = 0; i < A.size(); ++i) {
		prod *= x - A.at(i);
	}
	return prod;
}


template <typename T>
std::vector<T> eigenvals(matrix<T> B, double x = 0) {
	std::vector<T> vals = {};
	T val = 0;
	while (!isinf(val)) {
		while (abs(eigenfunc(x, B)) > std::pow(10, -10)) {
			if (!vals.empty()) {
				double g = product(x, vals);
				double g_p = product_p(x, vals);
				double v = eigenfunc(x, B) * g / (detderive(x, B) * g - eigenfunc(x, B) * g_p);
				x -= eigenfunc(x, B) * g / (detderive(x, B) * g - eigenfunc(x, B) * g_p);
				val = x;
 			}
			else {
				x -= eigenfunc(x, B) / detderive(x, B);
				val = x;
			}
			if (isnan(val)) { break; }
		}
		if (isnan(val)) { break; }
		if (!isinf(val)) {
			vals.push_back(val);
		}
		x = 0;
	}
	return vals;
}


template <typename T>
matrix<T> eigenvecs(matrix<T> B) {
	matrix<T> C({ eigenvals(B) });
	std::vector<std::vector<T>> X(B.m);
	for (int i = 0; i < B.m; ++i) {
		X[i].resize(B.m);
	}
	for (int i = 1; i <= B.m; ++i) {
		X[i-1] = solve_row_echelon(gaussian_elimination(B - C.get(1, i) * I((T)B.m))).X[0];
	}
	matrix<T> A(X);
	return T(A);
}


template <typename T>
matrix<T> gaussian_elimination(matrix<T> B) {
	matrix<T> X = B;
	for (int j = 1; j <= X.n; ++j) {
		for (int i = j + 1; i <= X.m; ++i) {
			if (j != X.n) {
				double scale = X.get(i, j) / X.get(j, j);
				for (int l = j; l <= X.n; ++l) {
					X.rep(i, l, X.get(i, l) - X.get(j, l) * scale);
				}
			}
		}
	}
	return X;
}


template <typename T>
matrix<T> solve_row_echelon(matrix<T> B) {
	matrix<T> X = matrix_n(1, B.n, 0.0);
	int first_nonzero_idx = 0;
	while (B.X[B.m - (B.get(B.m, B.n) > 0.000001 || B.get(B.m, B.n) < -0.000001) != 0 ? 1 : 2][first_nonzero_idx] > 0.000001 || B.X[B.m - (B.get(B.m, B.n) > 0.000001 || B.get(B.m, B.n) < -0.000001) != 0 ? 1 : 2][first_nonzero_idx] < -0.000001) {
		first_nonzero_idx += 1;
	}
	if (first_nonzero_idx != B.n) {
		!(B.get(B.m, B.n) > 0.000001 || B.get(B.m, B.n) < -0.000001) ? X.rep(1, B.m, 1) : X.rep(1, B.m, B.get(B.m, B.n));
		for (int i = 1; i < B.m; ++i) {
			double sum = 0;
			for (int j = 0; j < i; ++j) {
				sum += B.get(B.m - i, B.n - j) * X.get(1, B.n - j);
			}
			X.rep(1, B.m - i, -sum / B.get(B.m - i, B.m - i));
		}
	}
	return X;
}


// Eigendecomposition of matrix
template <typename T>
std::vector<matrix<T>> ED(matrix<T> B) {
	matrix<T> E = eigenvecs(B);
	matrix<T> E_p = inv(E);
	std::vector<std::vector<T>> V(B.m);
	for (int i = 0; i < B.m; ++i) {
		V[i].resize(B.m);
	}
	for (int i = 0; i < B.m; ++i) {
		V[i][i] = eigenvals(B)[i];
	}
	matrix<T> D({ V });
	std::vector<matrix<T>> X{ E, D, E_p };
	return X;
}


template <typename T>
std::vector<T> proj(std::vector<T> u, std::vector<T> a) {
	std::vector<T> x = u;
	double v = inner_product(u.begin(), u.end(), u.begin(), 0);
	double scalar = inner_product(u.begin(), u.end(), a.begin(), 0) / (double)inner_product(u.begin(), u.end(), u.begin(), 0);
	for (int i = 0; i < u.size(); ++i) {
		x.at(i) *= scalar;
	}

	return x;
}


template <typename K>
std::vector<matrix<K>> QR(matrix<K> B) {
	std::vector<matrix<K>> X;
	std::vector<std::vector<K>> a = T(B).X;
	std::vector<std::vector<K>> U;
	std::vector<K> u = a[0];
	U.push_back(u);
	for (int i = 1; i < B.n; ++i) {
		u = a[i];
		for (int j = 0; j < i; ++j) {
			for (int l = 0; l < u.size(); ++l) { u[l] -= proj(U[j], a[i])[l];}
		}
		U.push_back(u);
	}
	std::vector<std::vector<K>> E;
	for (int i = 0; i < U.size(); ++i) {
		E.push_back(U[i]);
		double v = inner_product(U[i].begin(), U[i].end(), U[i].begin(), 0);
		for (int l = 0; l < E[i].size(); ++l) { E[i][l] /= sqrt(inner_product(U[i].begin(), U[i].end(), U[i].begin(), 0.0)); }
	}

	matrix<K> Q_T(E);
	matrix<K> R = Q_T * B;

	X.push_back(T(Q_T));
	X.push_back(R);

	return X;
}


int fact(int n) {
	return (n == 0) || (n == 1) ? 1 : n * fact(n - 1);
}


template <typename T>
matrix<T> PowerOfMatrix(matrix<T> B, int n) {
	matrix<T> D = ED(B).at(1);
	for (int i = 1; i <= B.m; ++i) {
		D.rep(i, i, std::pow(D.get(i, i), n));
	}
	matrix<T> X = ED(B).at(0) * D * ED(B).at(2);
	return X;
}


// Takes the exponential of matrix
template <typename T>
matrix<T> ExponentialOfMatrix(matrix<T> B) {
	matrix<T> D = matrix_n(B.m, B.n, 0.0);

	for (int i = 1; i <= B.m; ++i) {
		D.rep(i, i, std::exp(ED(B).at(1).get(i, i)));
	}
	
	return ED(B).at(0) * D * ED(B).at(2);
}


// prints a matrix's entries with a precision(default is 5).
template <typename T>
void show(matrix<T> B, int precision = 5) {
	int current_length = 0;
	std::vector<T> values(B.n * B.m);
	for (int j = 1; j <= B.n; ++j) {
		for (int k = 1; k <= B.m; ++k) {
			if (B.get(k, j) == 0) { values[(j - 1) * B.m + k - 1] = 2 + precision; }
			else if (B.get(k, j) > 0) { values[(j - 1) * B.m + k - 1] = trunc(log10(abs(B.get(k, j)))) + 2 + precision; }
			else if (B.get(k, j) < -1) { values[(j - 1) * B.m + k - 1] = trunc(log10(abs(B.get(k, j)))) + 3 + precision; }
			else { values[(j - 1) * B.m + k - 1] = 3 + precision; }
		}
	}
	int max_length = *std::max_element(values.begin(), values.end());

	std::cout << "|";
	for (int i = 1; i < max_length * B.n + 5 * B.n; ++i) { std::cout << "-";}
	std::cout << "|";

	std::cout << "\n|";
	for (int i = 1; i < max_length * B.n + 5 * B.n; ++i) { std::cout << " "; }
	std::cout << "|\n";
 
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.n; ++j) {
			if (j == 1) { std::cout << "| "; }
			if (precision != 0) {
				current_length = 0;
				if (0 <= B.get(i, j) && B.get(i, j) < 10) {
					current_length = 7;
				}
				else if (-10 < B.get(i, j) && B.get(i, j) < 0) {
					current_length = 8;
				}
				else {
					if (B.get(i, j) > 0) { current_length = trunc(log10(abs(B.get(i, j)))) + 2 + precision; }
					else if (B.get(i, j) < 0) { current_length = trunc(log10(abs(B.get(i, j)))) + 3 + precision; }
					else { current_length = 2 + precision; }
				}
				if (j == 1) { std::cout << "|"; }
				for (int l = 1; l <= max_length - current_length; ++l) {
					std::cout << " ";
				}
				std::cout << std::setw(precision + 2) << std::setfill(' ') << std::fixed << std::setprecision(precision) << B.get(i, j);
			}
			else {
				current_length = 0;
				if (0 <= B.get(i, j) && B.get(i, j) < 10) {
					current_length = 7;
				}
				else if (-10 < B.get(i, j) && B.get(i, j) < 0) {
					current_length = 8;
				}
				else {
					if (B.get(i, j) > 0) { current_length = trunc(log10(abs(B.get(i, j)))) + 2 + precision; }
					else if (B.get(i, j) < 0) { current_length = trunc(log10(abs(B.get(i, j)))) + 3 + precision; }
					else { current_length = 2 + precision; }
				}
				std::cout << std::setw(1) << std::setfill(' ') << std::fixed << std::setprecision(0) << B.get(i, j);
			}
			if (j != B.n) { 
				std::cout << "|   |";
			}
			else { std::cout << "| |"; }
		}
		std::cout << "\n|";
		for (int i = 1; i < max_length * B.n + 5 * B.n; ++i) { std::cout << " "; }
		std::cout << "|\n";
	}
	std::cout << "|";
	for (int i = 1; i < max_length * B.n + 5 * B.n; ++i) { std::cout << "-"; }
	std::cout << "|";
	std::cout << "\n\n";
}

// multiplies the matrix with itself n times and returns that matrix.
template <typename T>
matrix<T> operator^(matrix<T> B, int n) {
	matrix<T> G = B;
	for (int i = 1; i < n; ++i) {
		G = G * B;
	}
	return G;
}


// returns the product of matrices. (be aware that order matters (i.e A*B!=B*A))
template <typename T>
matrix<T> operator*(matrix<T> a, matrix<T> b) {
	if (a.n == b.m) {
		T sum{};
		matrix<T> X = matrix_n(a.m, b.n, a.get(1, 1));
		for (int i = 1; i <= a.m; ++i) {
			for (int j = 1; j <= b.n; ++j) {
				sum = 0;
				for (int k = 1; k <= b.m; ++k) {
					sum += a.get(i, k) * b.get(k, j);
				}
				X.rep(i, j, sum);
			}
		}
		return X;
	}
}

// returns the sum of matrices.
template <typename T>
matrix<T> operator+(matrix<T> a, matrix<T> b) {
	if (a.m == b.m && a.n == b.n) {
		matrix<T> X = a;
		for (int i = 1; i <= a.m; ++i) {
			for (int j = 1; j <= b.n; ++j) {
				X.rep(i, j, a.get(i, j) + b.get(i, j));
			}
		}
		return X;
	}
}

// returns the difference of matrices.
template <typename T>
matrix<T> operator-(matrix<T> a, matrix<T> b) {
	if (a.m == b.m && a.n == b.n) {
		matrix<T> X = a;
		for (int i = 1; i <= a.m; ++i) {
			for (int j = 1; j <= b.n; ++j) {
				X.rep(i, j, a.get(i, j) - b.get(i, j));
			}
		}
		return X;
	}
}


template <typename T>
matrix<T> operator*(T a, matrix<T> b) {
 	for (int i = 1; i <= b.m; ++i) {
		for (int j = 1; j <= b.n; ++j) {
			b.rep(i, j, (T)(a * b.get(i, j)));
		}
	}
	return b;
}

// above and below operators are the same. Reason for 2 seemingly same operator is order matters for operators so you can't do k*A and A*k at the same time.
// returns the scalar multiplication of matrix. 

template <typename T>
matrix<T> operator*(matrix<T> b, T a) {
	for (int i = 1; i <= b.m; ++i) {
		for (int j = 1; j <= b.n; ++j) {
			b.rep(i, j, (T)(a * b.get(i, j)));
		}
	}
	return b;
}

template <typename T>
matrix<T> operator/(matrix<T> b, T a) {
	for (int i = 1; i <= b.m; ++i) {
		for (int j = 1; j <= b.n; ++j) {
			b.rep(i, j, (T)(b.get(i, j) / a));
		}
	}
	return b;
}


template <typename T>
class cplx {
public:
	T R;
	T I;
	cplx<T>(T reel, T img) {
		R = reel;
		I = img;
	}
};

template <typename T>
class matrix_c {
public:
	std::vector<std::vector<T>> Reel;
	std::vector<std::vector<T>> Img;
	int m;
	int n;
	matrix<T> R;
	matrix<T> I;
	matrix_c<T>(std::vector<std::vector<T>> X, std::vector<std::vector<T>> Y) {
		m = X.size();
		n = X[0].size();
		Reel = X;
		Img = Y;
		R = matrix<T>(X);
		I = matrix<T>(Y);
	}

	T get_R(int i, int j) {
		return Reel[i - 1][j - 1];
	}

	void rep_R(int i, int j, T n) {
		Reel[i - 1][j - 1] = n;
	}

	T get_I(int i, int j) {
		return Img[i - 1][j - 1];
	}

	void rep_I(int i, int j, T n) {
		Img[i - 1][j - 1] = n;
	}

	cplx<T> get(int i, int j) {
		cplx<T> A(get_R(i, j), get_I(i, j));
		return A;
	}

	void rep(int i, int j, cplx<T> G) {
		rep_R(i, j, G.R);
		R.rep(i, j, G.R);
		rep_I(i, j, G.I);
		I.rep(i, j, G.I);
	}
};

template <typename T>
cplx<T> operator+(cplx<T> A, cplx<T> B) {
	cplx<T> S(A.R + B.R, A.I + B.I);
	return S;
}

template <typename T>
cplx<T> operator-(cplx<T> A, cplx<T> B) {
	cplx<T> S(A.R - B.R, A.I - B.I);
	return S;
}

template <typename T>
cplx<T> operator*(cplx<T> A, cplx<T> B) {
	cplx<T> S(A.R * B.R - A.I * B.I, A.R * B.I + A.I * B.R);
	return S;
}

template <typename T>
cplx<T> operator/(cplx<T> A, cplx<T> B) {
	cplx<T> S(A.R * B.R + A.I * B.I, -A.R * B.I + A.I * B.R);
	S = S / (B.R * B.R + B.I * B.I);
	return S;
}

template <typename T>
cplx<T> operator*(cplx<T> A, T k) {
	cplx<T> S(A.R * k, A.I * k);
	return S;
}

template <typename T>
cplx<T> operator*(T k, cplx<T> A) {
	cplx<T> S(A.R * k, A.I * k);
	return S;
}

template <typename T>
cplx<T> operator/(cplx<T> A, T k) {
	cplx<T> S(A.R / k, A.I / k);
	return S;
}

template <typename T>
cplx<T> operator/(T k, cplx<T> A) {
	cplx<T> S(A.R / k, A.I / k);
	return S;
}

template <typename T>
cplx<T> operator-(cplx<T> A) {
	cplx<T> S(-A.R, -A.I);
	return S;
}

template <typename T>
cplx<T> operator^(cplx<T> A, T n) {
	if (n == (T)1) {
		return A;
	}
	else if (n == (T)2) {
		return A * A;
	}
	else {
		return A * (A ^ (n - 1));
	}
}

template <typename K>
matrix_c<K> T(matrix_c<K> B) {
	matrix<K> R = B.R;
	matrix<K> I = B.I;
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.m; ++j) {
			R.rep(i, j, B.R.get(j, i));
		}
	}
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.m; ++j) {
			I.rep(i, j, B.I.get(j, i));
		}
	}
	matrix_c<K> G(R.X, I.X);
	return G;
}

template <typename T>
cplx<T> Tr(matrix_c<T> B) {
	cplx<T> sum(0, 0);
	for (int i = 1; i <= B.m; ++i) {
		sum = sum + cplx<T> (B.get_R(i, i), B.get_I(i, i));
	}
	return sum;
}

template <typename T>
matrix_c<T> M(int m, int n, matrix_c<T> B) {
	std::vector<std::vector<T>> R = B.Reel;
	std::vector<std::vector<T>> I = B.Img;
	int k = 0;
	for (int i = 1; i <= B.m; ++i) {
		if (i != m) {
			R[i - 1 - k].erase(R[i - 1 - k].begin() + n - 1);
			I[i - 1 - k].erase(I[i - 1 - k].begin() + n - 1);
		}
		else {
			R.erase(R.begin() + i - 1);
			I.erase(I.begin() + i - 1);
			k = 1;
		}
	}
	matrix_c<T> C(R, I);
	return C;
}

template <typename T>
cplx<T> det(matrix_c<T> B) {
	cplx<T> sum(0, 0);
	if (B.n == 2) {
		return B.get(1, 1) * B.get(2, 2) - B.get(1, 2) * B.get(2, 1);
	}
	else if (B.n == 1) {
		return B.get(1, 1);
	}
	for (int j = 1; j <= B.m; ++j) {
		sum = sum + (T)(-2 * ((1 + j) % 2) + 1) * B.get(1, j) * det(M(1, j, B));
	}
	return sum;
}

template <typename T>
matrix_c<T> adj(matrix_c<T> B) {
	matrix_c<T> G = B;
	if (B.m == 2) {
		G.rep(1, 1, B.get(2, 2));
		G.rep(1, 2, -B.get(1, 2));
		G.rep(2, 1, -B.get(2, 1));
		G.rep(2, 2, B.get(1, 1));
		return G;
	}
	else if (B.m == 1) {
		G.rep(1, 1, B.get(1, 1));
		return G;
	}
	else {
		for (int i = 1; i <= B.m; ++i) {
			for (int j = 1; j <= B.n; ++j) {
				G.rep(i, j, (T)(-2 * ((i + j) % 2) + 1) * det(M(j, i, B)));
			}
		}
		return G;
	}
}

template <typename T>
matrix_c<T> inv(matrix_c<T> B) {
	try {
		if (det(B).R == 0 && det(B).I == 0) {
			throw (0);
		}
		else {
			return adj(B) / det(B);
		}
	}
	catch (int n) {
		std::cout << "determinant of the matrix is zero therefore the matrix is singular(not-invertible)";
	}
}

template <typename T>
matrix_c<T> operator+(matrix_c<T> X, matrix_c<T> Y) {
	matrix<T> R;
	matrix<T> I;
	R = operator+(X.R, Y.R);
	I = operator+(X.I, Y.I);

	matrix_c<T> G(R.X, I.X);
	return G;
}

template <typename T>
matrix_c<T> operator-(matrix_c<T> X, matrix_c<T> Y) {
	matrix<T> R;
	matrix<T> I;
	R = operator-(X.R, Y.R);
	I = operator-(X.I, Y.I);

	matrix_c<T> G(R.X, I.X);
	return G;
}

template <typename T>
matrix_c<T> operator*(matrix_c<T> X, matrix_c<T> Y) {
	matrix<T> R;
	matrix<T> I;
	R = X.R * Y.R - X.I * Y.I;
	I = X.R * Y.I + X.I * Y.R;

	matrix_c<T> G(R.X, I.X);
	return G;
}

template <typename T>
matrix_c<T> operator*(matrix_c<T> X, T k) {
	matrix<T> R;
	matrix<T> I;
	R = operator*(X.R, k);
	I = operator*(X.I, k);

	matrix_c<T> G(R.X, I.X);
	return G;
}

template <typename T>
matrix_c<T> operator*(T k, matrix_c<T> X) {
	matrix<T> R;
	matrix<T> I;
	R = operator*(k, X.R);
	I = operator*(k, X.I);

	matrix_c<T> G(R.X, I.X);
	return G;
}

template <typename T>
matrix_c<T> operator*(matrix_c<T> X, cplx<T> k) {
	matrix<T> A = k.R * X.R - k.I * X.I;
	matrix<T> B = k.R * X.I + k.I * X.R;

	matrix_c<T> G(A.X, B.X);

	return G;
}

template <typename T>
matrix_c<T> operator/(matrix_c<T> X, cplx<T> k) {
	matrix<T> A = k.R * X.R + k.I * X.I;
	matrix<T> B = k.R * X.I - k.I * X.R;

	A = A / (k.R * k.R + k.I * k.I);
	B = B / (k.R * k.R + k.I * k.I);

	matrix_c<T> G(A.X, B.X);

	return G;
}

template <typename T>
void show(matrix_c<T> X, int precision = 5) {
	std::vector<std::vector<T>> A;
	std::vector<T> row;
	for (int i = 1; i <= X.m; ++i) {
		for (int j = 1; j <= X.n; ++j) {
			row.push_back(X.R.get(i, j));
		}
		for (int j = 1; j <= X.n; ++j) {
			row.push_back(X.I.get(i, j));
		}
		A.push_back(row);
		row.clear();
	}

	matrix<T> B(A);

	int dec = 0;
	int max_length;
	int max_length_i;
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.n / 2; ++j) {
			std::vector<T> values(B.n);
			std::vector<T> values_i(B.n);
			std::vector<T> values_m(B.m);
			for (int k = 1; k <= B.m; ++k) {
				if (B.get(k, j) > 1) { values[k - 1] = trunc(log10(abs(B.get(k, j)))) + 2 + precision; }
				else if (B.get(k, j) < -1) { values[k - 1] = trunc(log10(abs(B.get(k, j)))) + 3 + precision; }
				else if (B.get(k, j) > 0 && B.get(k, j) == 0) { values[k - 1] = 2 + precision; }
				else { values[k - 1] = 3 + precision; }
			}
			max_length = *std::max_element(values.begin(), values.end());

			for (int k = 1; k <= B.m; ++k) {
				if (abs(B.get(k, j + B.n / 2)) > 1) { values_i[k - 1] = trunc(log10(abs(B.get(k, j + B.n / 2)))) + 2 + precision; }
				else { values_i[k - 1] = 2 + precision; }
			}
			max_length_i = *std::max_element(values_i.begin(), values_i.end());


			if (precision != 0) {
				if ((0 < B.get(i, j) && B.get(i, j) < 10) || B.get(i, j) == 0) {
					dec = 7;
				}
				else if (-10 < B.get(i, j) && B.get(i, j) < 0) {
					dec = 8;
				}
				else {
					dec = trunc(log10(abs(B.get(i, j)))) + 2 + precision;
					if (B.get(i, j) < 0) { ++dec; }
				}
				for (int l = 1; l <= max_length - dec; ++l) {
					std::cout << " ";
				}

				std::cout << std::setw(precision + 2) << std::setfill(' ') << std::fixed << std::setprecision(precision) << B.get(i, j);

				if (B.get(i, j + B.n / 2) >= 0) { std::cout << " + "; }
				if (B.get(i, j + B.n / 2) < 0) { std::cout << " - "; }

				if (abs(B.get(i, j + B.n / 2)) > 1) {
					dec = trunc(log10(abs(B.get(i, j + B.n / 2)))) + 2 + precision;
				}
				else {
					dec = 7;
				}
				for (int l = 1; l <= max_length_i - dec; ++l) {
					std::cout << " ";
				}

				std::cout << std::setw(precision + 2) << std::setfill(' ') << std::fixed << std::setprecision(precision) << abs(B.get(i, j + B.n / 2));
				std::cout << "i";
			}
			else {
				if (0 <= B.get(i, j) && B.get(i, j) < 10) {
					dec = 7;
				}
				else if (-10 < B.get(i, j) && B.get(i, j) < 0) {
					dec = 8;
				}
				else {
					dec = trunc(log10(abs(B.get(i, j)))) + 2 + precision;
					if (B.get(i, j) < 0) { ++dec; }
				}
				for (int l = 1; l <= max_length - dec; ++l) {
					std::cout << " ";
				}
				std::cout << std::setw(1) << std::setfill(' ') << std::fixed << std::setprecision(0) << B.get(i, j);

				if (B.get(i, j + B.n / 2) > 0) { std::cout << " + "; }
				if (B.get(i, j + B.n / 2) < 0) { std::cout << " - "; }

				if (-10 < B.get(i, j + B.n / 2) && B.get(i, j + B.n / 2) < 10) {
					dec = 7;
				}
				else {
					dec = trunc(log10(abs(B.get(i, j + B.n / 2)))) + 2 + precision;
					if (B.get(i, j) < 0) { ++dec; }
				}
				for (int l = 1; l <= max_length_i - dec; ++l) {
					std::cout << " ";
				}

				std::cout << std::setw(1) << std::setfill(' ') << std::fixed << std::setprecision(0) << abs(B.get(i, j + B.n / 2));
				std::cout << "i";
			}
			if (j != B.n / 2) {
				std::cout << "   ";
			}
		}
		std::cout << "\n\n";
	}
}

template <typename T>
void show(cplx<T> X, int precision = 5) {
	std::cout << std::setw(precision + 2) << std::setfill(' ') << std::fixed << std::setprecision(precision) << X.R;
	if (X.I > 0) { std::cout << " + "; }
	else { std::cout << " - "; }
	std::cout << abs(X.I) << "i" << "\n\n";
}
