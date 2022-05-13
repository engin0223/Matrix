#pragma once
#include <vector>
#include<algorithm>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <string>
#include <cmath>
#include <sstream> 
#include <string> 


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

};


template <class C>
#define name(C) #C  // gets the variable name of class C object.


void s() {}




// Creates a ixj matrix whose entries are "type(T) varible_name(n)".(E.g. "double number").
template <typename T>
matrix<T> matrix_n(int i, int j, T n) {
	std::vector<std::vector<T>> B;
	std::vector<T> row(j, n);
	for (int k = 0; k < i; ++k) {
		B.push_back(row);
	}
	return matrix<T>(B);
}


// Creates a NxN matrix Identity matrix.(i.e matrix[i][i] = 1 for all i that is 1<=i<=N and all other entries are 0)
template <typename T>
matrix<T> I(int N) {
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
		sum += (-2 * ((1 + j) % 2) + 1) * B.get(1, j) * get_minor(1, j, B);
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
		G.rep(1, 1, 1);
		return G;
	}
	else {
		for (int i = 1; i <= B.m; ++i) {
			for (int j = 1; j <= B.n; ++j) {
				G.rep(i, j, (-2 * ((i + j) % 2) + 1) * get_minor(j, i, B));
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
			return adj(B) * (pow(det(B), -1.0)) * 1.0;
		}
	}
	catch (int n) {
		std::cout << "determinant of the matrix is zero therefore the matrix is singular(not-invertible)";
	}
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


// prints a matrix's entries with a precision(default is 5).
template <typename T>
void show(matrix<T> B, int precision = 5) {
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.n; ++j) {
			std::vector<T> values(B.n);
			for (int k = 1; k <= B.n; ++k) {
				values[k - 1] = B.get(k, j);
			}
			int max_length = 7;
			int max_num = *std::max_element(values.begin(), values.end());
			if (max_num == 0) {
				for (int k = 1; k <= B.n; ++k) {
					if (values[k - 1] < 0) { max_length = 8; break; }
				}
			}
			else {
				if (0 <= max_num < 10) {
					max_length = 7;
				}
				else if (-10 < max_num < 0) {
					max_length = 8;
				}
				else {
					max_length = trunc(log10(*std::max_element(values.begin(), values.end(),
						[](const int& a, const int& b) {return abs(a) < abs(b); }))) + 2 + precision;
				}
			}
			if (precision != 0) {
				std::cout << std::setw(precision + 2) << std::setfill(' ') << std::fixed << std::setprecision(precision) << B.get(i, j);
			}
			else {
				std::cout << std::setw(1) << std::setfill(' ') << std::fixed << std::setprecision(0) << B.get(i, j);
			}
			if (j != B.n) { 
				int dec = 0;
				if (0 <= B.get(i, j) && B.get(i, j) < 10) {
					dec = 7;
				}
				else if (-10 < B.get(i, j) && B.get(i, j) < 0) {
					dec = 8;
				}
				else {
					dec = trunc(log10(abs(B.get(i, j)))) + 2 + precision;
				}
				for (int l = 1; l <= max_length - dec; ++l) {
					std::cout << " ";
				}
				std::cout << " - ";
			}
		}
		std::cout << std::endl;
	}
}

// multiplies the matrix with itself n times and returns that matrix.
template <typename T>
matrix<T> mul_n(matrix<T> B, int n) {
	matrix<T> G = B;
	for (int i = 1; i < n; ++i) {
		G = G * B;
	}
	return G;
}


// returns the (m,n) minor of the matrix.
template <typename T>
T get_minor(int m, int n, matrix<T> B) {
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
	T M = det(C);
	return M;
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
matrix<double> operator*(double a, matrix<T> b) {
 	for (int i = 1; i <= b.m; ++i) {
		for (int j = 1; j <= b.n; ++j) {
			b.rep(i, j, a * b.get(i, j) * 1.0);
		}
	}
	return b;
}

// above and below operators are the same. Reason for 2 seemingly same operator is order matters for operators so you can't do k*A and A*k at the same time.
// returns the scalar multiplication of matrix. 

template <typename T>
matrix<double> operator*(matrix<T> b, double a) {
	for (int i = 1; i <= b.m; ++i) {
		for (int j = 1; j <= b.n; ++j) {
			b.rep(i, j, a * b.get(i, j) * 1.0);
		}
	}
	return b;
}
