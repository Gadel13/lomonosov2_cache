#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <random>
#include <chrono>


using namespace std;

 
//all sizes in KB
#define L1size 32
#define L2size 256
#define L3size 35840

#define NUM_EVENTS 3
#define datasize 100000000

#define fiMin 0
#define fiMax datasize-1

#define fMin -1000000.0
#define fMax 1000000.0


void cache_test(double* &A, int minute) {	
	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);


	double sum = 0;

	int num_sum = 110000000;

	int* tmpi = new int[num_sum];
	for (int i = 0; i < num_sum; i++){
		tmpi[i] = f(gen);
	}


	auto begin = std::chrono::high_resolution_clock::now();
	for(int j = 0; j < minute; j++){
		for(int i = 0; i < num_sum; i++) {
			sum += A[tmpi[i]];
		}
	}
	auto end = std::chrono::high_resolution_clock::now();

	delete[] tmpi;


	auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

	ofstream REZ("REZ_without_papi_double_1sum.csv", ios::in|ios::app);
	REZ << "double; "  << sum << ';' << datasize << ';' << num_sum << ';' << time/(1000) << endl ;
	REZ.close();
}

void trash_to_cache(double* &trash) {

	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);

	clock_t start,stop;

	double sum = 0;

	start = clock();
	for(int i = 0; i < datasize; i++) {
		sum = sum + trash[f(gen)] - trash[f(gen)] + trash[f(gen)] - trash[f(gen)];
	}
	stop = clock();

	long double time = (long double)(stop-start)/CLOCKS_PER_SEC;

	ofstream REZ("REZ_without_papi_trash.csv", ios::in|ios::app);
	REZ << sum << ";" << datasize << ';' << datasize << ';' << time << endl ;
	REZ.close();
}

int main (int argc, char** argv) // time min
{
	mt19937 gen(time(0)); 
    uniform_real_distribution<> f(fMin, fMax); 

    mt19937 geni(time(0)); 
    uniform_int_distribution<int> fi(fiMin, fiMax);

	double *data = new double[datasize];
	double *trash;


	for(int i = 0; i < datasize; i++) {
		data[i] = f(gen);
	}

	int minute = atoi(argv[1]);

	for(double k = 8; k <= 8; k *= 2) {
		trash = new double[datasize];
		for(int k = 0; k < datasize; k++) {
			trash[k] = f(gen);
		}
		cout << data[fi(gen)] + trash[fi(gen)];

		cache_test(data, minute);

		for(int j = 0; j < 5; j++) {
			trash_to_cache(trash);
		}

		cout << data[fi(gen)] + trash[fi(gen)] << endl;
		delete[] trash;
	}

	delete[] data;

	return 0;
}
