#include "Matrix.h"
#include <stdio.h>
#include <cstdlib>

Matrix::Matrix(double**matrix_, int m_, int n_):matrix(matrix_),m(m_),n(n_){
}

Matrix* Matrix::operator*(Matrix that) const{
	if(that.m!= this->n){
		perror("inconsistent matrix multiplication!");
		exit(-1);
	}
	double** newMatrix = new double*[this->m];
	for(int i = 0; i < this->m; i++){
		newMatrix[i] = new double[that.n];
	}
	for(int i = 0; i < this->m; i++)
		for(int j = 0; j < that.n; j++){
			double x = 0;
			for(int k = 0; k < this->n; k++){
				x+=this->matrix[i][k]*that.matrix[k][j];
			}
			newMatrix[i][j] = x;
		}
	return new Matrix(newMatrix, this->m, that.n);
}

Matrix* Matrix::operator-(Matrix that){
	if(this->m!= that.m || this->n != that.n){
		perror("inconsistent matrix subtraction!");
		exit(-1);
	}
	for(int i = 0; i < this->m; i++)
		for(int j =0; j < this->n; j++){
			this->matrix[i][j]-= that.matrix[i][j];
		}
	return this;
}

double*& Matrix::operator[](int i)const{
	if(i >= this->m){
		perror("illegal subscript!");
		exit(-1);
	}
	return this->matrix[i];
}

Matrix::~Matrix(){
	 // for(int i = 0; i < m; i++){
	 // 	delete this->matrix[i];
	 // }
	 // delete this->matrix;
}
