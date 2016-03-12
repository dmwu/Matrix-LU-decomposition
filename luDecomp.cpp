#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include "Matrix.h"
struct Pair{ double val; int pos; };
#pragma omp declare reduction (max: struct Pair: omp_out = omp_out.val > omp_in.val? omp_out:omp_in)
void print(double** data, int n,const char*);
double** formulatep(int* p, int n);
double** formulatel(double** data, int n);
double** formulateu(double** data, int n);
double l21normAndDestruct(double**,double**,double**,double**,int);
void freeMatrix(double**data, int n);
int PARALLEL_Limit=100;
int main (int argc, char *argv[]){
    int n, nthreads, seed,tid,maxRow;
    struct drand48_data drand_buf;
    double t0, t1,pivot;
    n  = atoi(argv[1]);
    nthreads = atoi(argv[2]);

    int* P = new int[n];
    for(int i = 0; i < n; i++)
    	P[i] = i;
    omp_set_num_threads (nthreads);
    double** data = new double*[n];
    double** dataCopy = new double*[n];
    t0 = omp_get_wtime();

    /*generate matrix in parallel*/
	#pragma omp parallel private(tid, seed, drand_buf) \
	  	  firstprivate(n,nthreads) \
	  	  shared(data,dataCopy)

	{
		  tid = omp_get_thread_num();	 
		  seed = tid*111237;
		  srand48_r(seed, &drand_buf);
		  double xx =0;
		  #pragma omp parallel for
		  for (int i= 0; i < n; i++){
		  	data[i] = new double[n];
		  	dataCopy[i] = new double[n];
		  	for(int j = 0; j < n; j++){
		      drand48_r(&drand_buf, &xx);
		      data[i][j] = xx;
		      dataCopy[i][j] = xx;
		  	}
		  }
	}

	
    /*start LU decomposition*/
    for(int k = 0; k < n; k++){

	 	struct Pair pair;
	 	pair.val = data[k][k];
	 	pair.pos = k;
	 	/*a parallel region that pivoting*/
	 	#pragma omp parallel if(n-k >= PARALLEL_Limit) firstprivate(k,n) shared(data) 
	 	{
	 	 	#pragma omp parallel for reduction(max:pair) 
		 	for(int i = k; i < n; i++){
		   		if( pair.val < std::abs(data[i][k]) ){
		   			pair.val = std::abs(data[i][k]);
		   			pair.pos = i;
		   		}
		   		printf("i'm ID %d i find max row %d, %f\n", omp_get_thread_num(), pair.pos, pair.val);
		 	}

		}
		
		if(pair.val <=0){
			perror("singular matrix");
			exit(-1);
		}
		maxRow = pair.pos;
		printf("max row %d, val %f\n", pair.pos,pair.val);
			/*
			 * swap permutation matrix P pivoting data matrix data	
			 */
		print(data,n,"data before swap");
		if (maxRow != k){
			int temp = P[k];
			P[k] = P[maxRow];
			P[maxRow] = temp;
			#pragma omp parallel for schedule(static)
			for(int i = 0; i < n; i++){
				double tempp = data[k][i];
				data[k][i] = data[maxRow][i];
				data[maxRow][i] = tempp;
			}			
		}

		print(data,n,"data after swap");
		
		pivot = data[k][k];
		#pragma omp parallel if(n-k >= PARALLEL_Limit) firstprivate(k,n,pivot) shared(data)
		{
			
			#pragma omp for
			for(int i = k+1; i < n; i++){
				data[i][k] = data[i][k]/pivot;
			}
		}

		#pragma omp parallel if(n-k >= PARALLEL_Limit) firstprivate(k,n,pivot) shared(data)
		{ 	
			#pragma omp for
			for(int i = k+1; i < n; i++)
				for(int j = k+1; j < n; j++){
					data[i][j] -= data[k][j]*data[i][k];
			} 
		}

 	}

 	t1 = omp_get_wtime();
  	printf("Time: %7.2f\n", t1-t0);
 	double** pmatrix = formulatep(P,n);
 	double** lmatrix = formulatel(data,n);
 	double** umatrix = formulateu(data,n);

 	printf("l21norm: %f\n", l21normAndDestruct(pmatrix, dataCopy,lmatrix,umatrix,n));
 	freeMatrix(data,n);
 	freeMatrix(pmatrix,n);
 	freeMatrix(lmatrix,n);
 	freeMatrix(umatrix,n);
 	freeMatrix(dataCopy,n);
  	return 0;
};

void freeMatrix(double**data, int n){
	for (int i = 0; i < n; i++){
		delete[] data[i];
	}
	delete[] data;
}
void print(double** data, int n,const char* str){
	printf("===============================\n");
	printf("matrix for %s\n",str);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			printf("%f ",data[i][j]);
	    }
	    printf("\n");
	}
}

double** formulatep(int* p, int n){
	double** pmatrix = new double*[n];
	for(int i = 0; i < n; i++){
		pmatrix[i] = new double[n]();
		pmatrix[i][p[i]]=1;
	}
	return pmatrix;
}

double** formulatel(double** data, int n){
	double** lmatrix = new double*[n];
	for(int i = 0; i < n; i++){
		lmatrix[i] = new double[n];
		for(int j = 0; j < n; j++){
			if(j < i) lmatrix[i][j] = data[i][j];
			else if(j == i) lmatrix[i][j] = 1;
			else lmatrix[i][j] = 0;
		}
	}
	return lmatrix;
}


double** formulateu(double** data, int n){
	double** umatrix = new double*[n];
	for(int i = 0; i < n; i++){
		umatrix[i] = new double[n];
		for(int j = 0; j < n; j++){
			if(j < i) umatrix[i][j] = 0;
			else umatrix[i][j] = data[i][j];
		}
	}
	return umatrix;
}

double l21normAndDestruct(double** p, double** a, double** l, double** u, int n){
	Matrix* pp = new Matrix(p,n,n);
	Matrix* aa = new Matrix(a,n,n);
	Matrix* ll = new Matrix(l,n,n);
	Matrix* uu = new Matrix(u,n,n);
	Matrix* res1 = (*pp)*(*aa);
	Matrix* res2 = (*ll)*(*uu);
	res1 = (*res1)-(*res2);
	print(res1->matrix,n,"pa");
	print(res2->matrix,n,"lu");
	double ret = 0;
	for (int j = 0; j < n; j++){
		double col = 0;
		for(int i = 0; i < n; i++){
			col += (*res1)[i][j] * (*res1)[i][j];
		}
		ret+= sqrt(col);
	}
	return ret;
}

