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

template <typename T>
class matrix {
public:
	std::vector<std::vector<T>> X;
	int m;
	int n;
	matrix<T>(std::vector<std::vector<T>> B) {
		m = int(B.size());
		if (not std::is_same<std::type_info(B[0]), std::vector<T>>) {
			n = m;
			m = 1;
		}
		else {
			n = int(B[0].size());
		}
		X = B;
	};

	T get(int i, int j) {
		if (m != 1){ return X[i - 1][j - 1]; }
		else { return X[j - 1]; }
	}

	void rep(int i, int j, T n) {
		if (m != 1) { X[i - 1][j - 1] = n; }
		else { X[j - 1] = n; }
	}

};


template <typename T>
matrix<T> matrix_n(int i, int j, T n) {
	std::vector<std::vector<T>> B;
	std::vector<int> row(j, n);
	for (int k = 0; k < i; ++k) {
		B.push_back(row);
	}
	return matrix<T>(B);
}

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
		sum += (-2 * ((1 + j) % 2) + 1) * B.get(1, j) * det(get_minor(1, j, B));
	}
	return sum;
}


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
				G.rep(i, j, (-2 * ((i + j) % 2) + 1) * det(get_minor(j, i, B)));
			}
		}
		return G;
	}
}

template <typename T>
matrix<T> inv(matrix<T> B) {
	return adj(B) * (pow(det(B), -1.0)) * 1.0;
}


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

template <typename T>
T Tr(matrix<T> B) {
	T sum{};
	for (int i = 1; i <= B.m; ++i) {
		sum += B.get(i, i);
	}
	return sum;
}


template <typename T>
void show(matrix<T> B, int precision=5) {
	for (int i = 1; i <= B.m; ++i) {
		for (int j = 1; j <= B.n; ++j) {
			if (precision != 0) {
				std::cout << std::setw(precision + 2) << std::setfill(' ') << std::fixed << std::setprecision(precision) << B.get(i, j);
			}
			else {
				std::cout << std::setw(1) << std::setfill(' ') << std::fixed << std::setprecision(0) << B.get(i, j);
			}
			if (j != B.n) { std::cout << " - "; }
		}
		std::cout << std::endl;
	}
}

template <typename T>
matrix<T> mul_n(matrix<T> B, int n) {
	matrix<T> G = B;
	for (int i = 1; i < n; ++i) {
		G = G * B;
	}
	return G;
}


template <typename T>
matrix<T> get_minor(int m, int n, matrix<T> B) {
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


template <typename T>
matrix<double> operator*(matrix<T> b, double a) {
	for (int i = 1; i <= b.m; ++i) {
		for (int j = 1; j <= b.n; ++j) {
			b.rep(i, j, a * b.get(i, j) * 1.0);
		}
	}
	return b;
}



class complex {
public:
	double reel;
	double img;
	complex(double a, double b) {
		reel = a;
		img = b;
	}
};



class matrix_cx {
public:
	std::vector<std::vector<complex>> X;
	std::vector<std::vector<double>> A;
	std::vector<std::vector<double>> B;
	int m;
	int n;
	matrix_cx(std::vector<std::vector<complex>> D) {
		m = int(D.size());
		X = D;
		if (typeid(D[0]).name() == typeid(std::vector<complex>).name()) {
			n = m;
			m = 1;
		}
		else {
			n = int(D[0].size());
		}
		std::vector<double> R;
		std::vector<double> I;
		for (int i = 1; i <= m; ++i) {
			for (int j = 1; j <= n; ++j) {
				R.push_back(X[i - 1][j - 1].reel);
				I.push_back(X[i - 1][j - 1].img);
			}
			A.push_back(R);
			B.push_back(I);
			R = {};
			I = {};
		}
	};

	complex get(int i, int j) {
		return X[i - 1][j - 1];
	}

	void rep(int i, int j, complex n) {
		X[i - 1][j - 1] = n;
	}

};
