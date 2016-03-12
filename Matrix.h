#ifndef MAXTRIX_H
#define MAXTRIX_H

class Matrix{
public:
	int m, n;
	double** matrix;
	Matrix(double** matrix, int m_, int n_);
	Matrix* operator*(Matrix that) const;
	Matrix* operator-(Matrix that);
	double*& operator[](int i) const;
	~Matrix();
};

#endif