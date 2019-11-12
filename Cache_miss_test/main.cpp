#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <papi.h>
#include <random>


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


void cache_test(double* &A, int num_sum) {	
	int Events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM};
	long long values[NUM_EVENTS];
	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);

	clock_t start,stop;

	double sum = 0;

	if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
		cout << "PAPI ERROR";

	start = clock();
	for(int i = 0; i < num_sum; i++) {
		sum = *((long long*) &sum) | *((long long*) &A[f(gen)]);
		//sum = sum + A[f(gen)];// + A[f(gen)] - A[f(gen)] - A[f(gen)];
	}
	stop = clock();
	if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
		cout << "PAPI ERROR";


	long double time = (long double)(stop-start)/CLOCKS_PER_SEC;

	ofstream REZ("REZ.csv", ios::in|ios::app);
	REZ << sum << ';' << datasize << ';' << num_sum << ';' << time << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << values[2]/time << endl ;
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
    uniform_real_distribution<> f(fMin, fMax); 

    mt19937 geni(time(0)); 
    uniform_int_distribution<int> fi(fiMin, fiMax);

	double *data = new double[datasize];
	double *trash;


	for(int i = 0; i < datasize; i++) {
		data[i] = f(gen);
	}


	for (int i = 0; i < atoi(argv[1]); ++i) {
		trash = new double[datasize];
		for(int k = 0; k < datasize; k++) {
			trash[k] = f(gen);
		}
		cout << data[fi(gen)] + trash[fi(gen)];

		cache_test(data, (L3size/8) * 1024);

		for(int j = 0; j < 5; j++) {
			trash_to_cache(trash);
		}


		cout << data[fi(gen)] + trash[fi(gen)] << endl;
		delete[] trash;
	}




	delete[] data;

	return 0;
}
