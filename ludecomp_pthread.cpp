/*
 * ludecomp.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: dingming
 */
#include <stdio.h>
#include <cstdlib>
#include <pthread.h>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include "Matrix.h"
void print(double** data, int n, const char*);
double** formulatep(int* p, int n);
double** formulatel(double** data, int n);
double** formulateu(double** data, int n);
double l21normAndDestruct(double**,double**,double**,double**,int);
void freeMatrix(double**data, int n);
void printToFile(double**, int, const char*);
struct Thread_data{
	int k;
	int n;
	int id;
	int lowIndex;
	int highIndex;
	double** shareData;
	void* optional;
	Thread_data(int k_, int n_, int id_, int low_, int high_, double** s_, void* opt=NULL):\
			k(k_),n(n_),id(id_), lowIndex(low_),highIndex(high_),shareData(s_),optional(opt){

	}
};

struct Pair{
	double val; int pos;
	Pair(int v, int p):val(v),pos(p){}
	Pair():val(0.0),pos(0){}
};

int PARALLEL_LIMIT=50;

void* initData_routine(void*);
void* searchMax_routine(void*);
void* updateL_routine(void*);
void* updateU_routine(void*);
void generateData(int n, int nthreads, double** data, double** dataCopy, pthread_t* threads);
int searchMax(int k, int n, int nthreads, double** data,pthread_t* threads);
void updateL(int k, int n, int nthreads,  double**data, pthread_t*threads);
void updateU(int k, int n, int nthreads,double** data, pthread_t*threads);

int main(int argc, char*argv[]){
    double t0, t1;
    struct Pair pair;
    int k,n,nthreads;
    n  = atoi(argv[1]);
    nthreads = atoi(argv[2]);
    struct timeval tv1,tv2;
    struct timezone tz;
    gettimeofday(&tv1,&tz);
    double** data = new double*[n];
    double** dataCopy = new double*[n];
    int* P = new int[n];
    for(int i = 0; i < n; i++)
    	P[i] = i;

    pthread_t* threads = new pthread_t[nthreads];
    generateData(n,nthreads,data,dataCopy,threads);
#ifdef VERBOSE
    print(data,n,"data");
    print(dataCopy,n,"dataCopy");
#endif
    for(k = 0; k < n; k++){
    	int maxRow = searchMax(k,n,nthreads,data,threads);
#ifdef VERBOSE
    	print(data,n,"data before swap");
    	printf("found maxRow: %d",maxRow);
#endif
    	if(maxRow!=k){
			int temp = P[k];
			P[k] = P[maxRow];
			P[maxRow] = temp;
			double* tempp = data[k];
			data[k] = data[maxRow] ;
			data[maxRow] = tempp;
    	}
#ifdef VERBOSE
    	print(data,n,"data after swap");
#endif

    	updateL(k,n,nthreads,data, threads);
#ifdef VERBOSE
    	print(data,n, "after updating L");
#endif
    	updateU(k,n,nthreads,data,threads);
#ifdef VERBOSE
    	print(data,n,"after updating U");
#endif
    }
    gettimeofday(&tv2,&tz);
    printf("time elapsed %7.2f\n",((tv2.tv_sec-tv1.tv_sec)*1000+(tv2.tv_usec-tv1.tv_usec)/1000)/1000.0);
    delete[] threads;
    double** pmatrix = formulatep(P,n);
 	double** lmatrix = formulatel(data,n);
 	double** umatrix = formulateu(data,n);
 	printf("l21norm: %f\n", l21normAndDestruct(pmatrix, dataCopy,lmatrix,umatrix,n));
 	freeMatrix(data,n);
 	freeMatrix(dataCopy,n);
    return 0;
}

void generateData(int n, int nthreads, double** data, double** dataCopy, pthread_t* threads){
		int lowIndex,highIndex;
	 	int avgDivision = n/nthreads;
	    pthread_attr_t attr;
	    pthread_attr_init(&attr);
	    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	    int rc;
	    void* status;
	    Thread_data** ttdata = new Thread_data*[nthreads];
	    /* generate random matrix in parallel,
	     * each thread is in charge of several rows
	    */
	    for(int i =0; i < nthreads; i++){
	    	lowIndex = avgDivision*i;
	    	highIndex = lowIndex>=(nthreads-1)*avgDivision ? n:lowIndex+avgDivision;
	    	ttdata[i] = new Thread_data(0,n,i,lowIndex,highIndex,data, dataCopy);
	    	rc = pthread_create(&threads[i], &attr ,initData_routine,(void*)ttdata[i]);
	        if (rc) {
	           printf("return code from pthread_create() for initData is %d\n", rc);
	           exit(-1);
	        }
	    }

	    /*  Wait for the other threads and free attribute*/

	    for(int i = 0; i < nthreads; i++) {
	       rc = pthread_join(threads[i], &status);
	       if (rc) {
	          printf("ERROR; return code from pthread_join() for initData is %d\n", rc);
	          exit(-1);
	          }
	       	  delete ttdata[i];
	    }
	    delete ttdata;
	    pthread_attr_destroy(&attr);

}

int searchMax(int k, int n, int nthreads, double** data,pthread_t* threads){
	   /*
	    * search for the row index of maximum abs value
	    * only parallelzes when #rows is larger than PARALLEL_LIMIT*/
	   struct Pair pair;
	   pair.val = 0;
	   pair.pos = k;
	   int avgDivision = (n-k)/nthreads;
	   int rc, lowIndex,highIndex;
	   if (avgDivision <= PARALLEL_LIMIT){
		  for(int i = k; i < n; i++){
			if( pair.val < std::fabs(data[i][k])){
				pair.val = std::fabs(data[i][k]);
				pair.pos = i;
			}
		  }
	   }
	   else{
		   /*
		    * every thread will store there local maximum to pairs[i]
			* and the master thread then find the global maximum via a for loop
			* */
		    Thread_data** ttdata = new Thread_data*[nthreads];
		    struct Pair* pairs = new Pair[nthreads]();
		    pthread_attr_t attr;
		    void* status;
		    pthread_attr_init(&attr);
		    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		    for(int i = 0; i < nthreads; i++){
			   lowIndex = k+i*avgDivision;
			   highIndex = lowIndex>=(nthreads-1)*avgDivision?n:lowIndex+avgDivision;
			   /*i is passed to the routine so that the routine can figure out
			    * the corresponding position it can modify */
			   ttdata[i] = new Thread_data(k,n,i,lowIndex, highIndex,data,pairs);
			   rc = pthread_create(&threads[i], &attr, searchMax_routine, (void*)ttdata[i]);
		       if (rc) {
		           printf("return code from pthread_create() for searchMax is %d\n", rc);
		           exit(-1);
		       }
		    }

		    for(int i = 0; i < nthreads; i++) {
		       rc = pthread_join(threads[i], &status);
		       if (rc) {
		          printf("ERROR; return code from pthread_join() for initData is %d\n", rc);
		          exit(-1);
		          }
		       delete ttdata[i];
		    }

		    for(int i = 0; i < nthreads; i++){
				if(pair.val < pairs[i].val){
					pair.val = pairs[i].val;
					pair.pos = pairs[i].pos;
				}
			 }
		     delete ttdata;
		     pthread_attr_destroy(&attr);

	    }

	   if(pair.val <=0){
		   perror("sigular matrix!");
		   exit(-1);
	   }
	   return pair.pos;
}

void updateL(int k, int n, int nthreads, double**data, pthread_t*threads){
	int lowIndex,highIndex;
 	int avgDivision = (n-k-1)/nthreads;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int rc;
    void* status;
    Thread_data** ttdata = new Thread_data*[nthreads];
    double pivot = data[k][k];
    if(avgDivision <= PARALLEL_LIMIT){
    	for(int i = k+1; i< n; i++){
    		data[i][k] /=pivot;
    	}
    }
    else{
		for(int i =0; i < nthreads; i++){
			lowIndex = avgDivision*i+1+k;
			highIndex = lowIndex>=(nthreads-1)*avgDivision +k+1? n:lowIndex+avgDivision;
			//printf("thread %d low %d high %d\n",i,lowIndex, highIndex);
			ttdata[i] = new Thread_data(k,n,i,lowIndex,highIndex,data, (void*)&pivot);
			rc = pthread_create(&threads[i], &attr ,updateL_routine,(void*)ttdata[i]);
			if (rc) {
			   printf("return code from pthread_create() for updaetL is %d\n", rc);
			   exit(-1);
			}
		}
    /*  Wait for the other threads and free attribute*/
		for(int i = 0; i < nthreads; i++) {
		   rc = pthread_join(threads[i], &status);
		   if (rc) {
			  printf("ERROR; return code from pthread_join() for updateL is %d\n", rc);
			  exit(-1);
			  }
			  delete ttdata[i];
		 }
    }
    delete[] ttdata;
    pthread_attr_destroy(&attr);

}

void updateU(int k, int n, int nthreads,double** data, pthread_t*threads){
	int lowIndex,highIndex;
 	int avgDivision = (n-k-1)/nthreads;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int rc;
    void* status;
    Thread_data** ttdata = new Thread_data*[nthreads];
    if(avgDivision <= PARALLEL_LIMIT/10){
    	for(int i = k+1; i < n; i++){
    		for(int j = k+1; j < n; j++){
    			data[i][j] -= data[k][j]*data[i][k];
    		}
    	}
    }
    else{
		for(int i =0; i < nthreads; i++){
			lowIndex = avgDivision*i+1+k;
			highIndex = lowIndex>=(nthreads-1)*avgDivision+k+1 ? n:lowIndex+avgDivision;
		//	printf("updateU thread %d low %d high %d\n",i,lowIndex, highIndex);
			ttdata[i] = new Thread_data(k,n,i,lowIndex,highIndex,data);
			rc = pthread_create(&threads[i], &attr ,updateU_routine,(void*)ttdata[i]);
			if (rc) {
			   printf("return code from pthread_create() for updaetU is %d\n", rc);
			   exit(-1);
			}
		}
    /*  Wait for the other threads and free attribute*/
		for(int i = 0; i < nthreads; i++) {
		   rc = pthread_join(threads[i], &status);
		   if (rc) {
			  printf("ERROR; return code from pthread_join() for updateU is %d\n", rc);
			  exit(-1);
			  }
			  delete ttdata[i];
		}
    }
    delete[] ttdata;
    pthread_attr_destroy(&attr);
}


void* initData_routine(void* threadData){
	Thread_data* tdata = (Thread_data*) threadData;
	int n = tdata->n;
	int low = tdata->lowIndex;
	int high = tdata->highIndex;
	int seed = tdata->id*11371;
	struct drand48_data drand_buf;
	srand48_r(seed, &drand_buf);
	double xx =0;
	for(int i = low; i < high; i++)
	{
	  	tdata->shareData[i] = new double[n];
	  	((double**)tdata->optional)[i] = new double[n];
		for(int j = 0; j < n; j++){
	      drand48_r(&drand_buf, &xx);
	      tdata->shareData[i][j] = xx;
	      ((double**)tdata->optional)[i][j] = xx;
		}
	}
	pthread_exit(NULL);
}


void* searchMax_routine(void* threadData){
	Thread_data* tdata = (Thread_data*)threadData;
	int low = tdata->lowIndex;
	int high = tdata->highIndex;
	int k = tdata->k;
	int i = tdata->id;
	Pair* pairs = (Pair*)tdata->optional;
	Pair pair = Pair(0,0);

	for(int j = low; j < high; j++){
		if(pair.val < std::fabs(tdata->shareData[j][k])){
			pair.val = std::fabs(tdata->shareData[j][k]);
			pair.pos = j;
		}
	}
	pairs[i].val = pair.val;
	pairs[i].pos = pair.pos;
	pthread_exit(NULL);
}

void* updateL_routine(void*threadData){
	Thread_data* tdata = (Thread_data*)threadData;
	int low = tdata->lowIndex;
	int high = tdata->highIndex;
	int k = tdata->k;
	double* pivot = (double*)(tdata->optional);
	//printf("updateL, k:%d, thread:%d, low:%d, high:%d, pivot:%f\n",k,tdata->id,low,high,(*pivot));
	for(int i = low; i< high; i++){
		tdata->shareData[i][k] /=(*pivot);
	}
	pthread_exit(NULL);
}

void* updateU_routine(void* threadData){
	Thread_data* tdata = (Thread_data*)threadData;
	int n= tdata->n;
	int low = tdata->lowIndex;
	int high = tdata->highIndex;
	int k = tdata->k;
	for(int i = low; i < high; i++){
		for(int j = k+1; j < n; j++){
			tdata->shareData[i][j] -= tdata->shareData[k][j]*tdata->shareData[i][k];
		}
	}
	pthread_exit(NULL);
}

void freeMatrix(double**data, int n){
	for (int i = 0; i < n; i++){
		delete[] data[i];
	}
	delete[] data;
}

void print(  double** data, int n, const char* str){
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
	Matrix* res1 = (*pp)*aa;
	Matrix* res2 = (*ll)*uu;
#ifdef VERBOSE
	print(pp->matrix,n,"p");
	print(aa->matrix,n,"a");
	print(ll->matrix,n,"l");
	print(uu->matrix,n,"u");
	print(res1->matrix,n,"pa");
	print(res2->matrix,n,"lu");
#endif

	res1 = (*res1)-res2;
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

void printToFile(double** data, int n, const char* file){
	FILE* fp = fopen(file,"w+");
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			fprintf(fp,"%f\n",data[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

}
