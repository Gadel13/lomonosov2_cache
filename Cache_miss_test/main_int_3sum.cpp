#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <papi.h>
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


void cache_test(int* &A, int num_sum) {	
	int Events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM};
	long long values[NUM_EVENTS];
	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);


	int sum = 0;

	int* tmpi = new int[3*num_sum];
	for (int i = 0; i < 3*num_sum; ++i)
		tmpi[i] = f(gen);

	if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
		cout << "PAPI ERROR";


	auto begin = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < num_sum; i++) {
		sum += A[tmpi[3*i]] + A[tmpi[3*i + 1]] + A[tmpi[3*i + 2]];
		//sum = sum + A[f(gen)];// + A[f(gen)] - A[f(gen)] - A[f(gen)];
	}
	auto end = std::chrono::high_resolution_clock::now();

	if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
		cout << "PAPI ERROR";

	delete[] tmpi;


	auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

	ofstream REZ("REZ_int_3sum.csv", ios::in|ios::app);
	REZ << "int; "   << sum << ';' << datasize << ';' << 3*num_sum << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
	REZ.close();
}

void trash_to_cache(double* &trash) {
	int Events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM};
	long long values_trash[NUM_EVENTS];

	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);

	clock_t start,stop;

	double sum = 0;

	if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
		cout << "PAPI ERROR";
	start = clock();
	for(int i = 0; i < datasize; i++) {
		sum = sum + trash[f(gen)] - trash[f(gen)] + trash[f(gen)] - trash[f(gen)];
	}
	stop = clock();
	if (PAPI_stop_counters(values_trash, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";
	long double time = (long double)(stop-start)/CLOCKS_PER_SEC;

	ofstream REZ("REZ_trash.csv", ios::in|ios::app);
	REZ << sum << ";" << datasize << ';' << datasize << ';' << time << ';'  << values_trash[0] << ';' << values_trash[1] << ';' <<  values_trash[2] << ';' << values_trash[2]/time << endl ;
	REZ.close();
}

int main (int argc, char** argv) // N - N times call cache_test
{
	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fMin, fMax); 

    mt19937 geni(time(0)); 
    uniform_int_distribution<int> fi(fiMin, fiMax);

	int *data = new int[datasize];
	double *trash;


	for(int i = 0; i < datasize; i++) {
		data[i] = f(gen);
	}

	for(double k = 1.0/32; k <= 8; k *= 2)
		for (int i = 0; i < atoi(argv[1]); ++i) {
			trash = new double[datasize];
			for(int k = 0; k < datasize; k++) {
				trash[k] = f(gen);
			}
			cout << data[fi(gen)] + trash[fi(gen)];

			cache_test(data, (L3size/4) * 1024 * k);

			for(int j = 0; j < 5; j++) {
				trash_to_cache(trash);
			}


			cout << data[fi(gen)] + trash[fi(gen)] << endl;
			delete[] trash;
		}




	delete[] data;

	return 0;
}
