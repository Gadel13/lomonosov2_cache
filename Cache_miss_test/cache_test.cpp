#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <papi.h>
#include <random>
#include <chrono>


using namespace std;

 
//all sizes in Bytes
#define L1size 32000
#define L2size 256000
#define L3size 35000000

#define number_trash_to_cache_calls 5
#define datasize 1e8
#define NUM_EVENTS 3

#define fiMin 0
#define fiMax datasize-1

#define fMin_double -1e-308
#define fMax_double 1e308

#define fMin_int -2e-9
#define fMax_int 2e9

template<typename data_type>
void cache_test(data_type* &A, char d_type, int num_sum_per_oper, int num_sum) {	
	int Events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM};
	long long values[NUM_EVENTS];

	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);


	data_type sum = 0;
	int imax = num_sum_per_oper*num_sum;

	int* tmpi = new int[imax];
	for (int i = 0; i < imax; ++i)
		tmpi[i] = f(gen);


	if(num_sum_per_oper == 1) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 2) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 3) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 4) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 5) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3] + A[tmpi[i] + 4];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 6) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3] + A[tmpi[i] + 4] + A[tmpi[i] + 5];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 7) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3] + A[tmpi[i] + 4] + A[tmpi[i] + 5] + A[tmpi[i] + 6];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 8) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3] + A[tmpi[i] + 4] + A[tmpi[i] + 5] + A[tmpi[i] + 6] + A[tmpi[i] + 7];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 9) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3] + A[tmpi[i] + 4] + A[tmpi[i] + 5] + A[tmpi[i] + 6] + A[tmpi[i] + 7] + A[tmpi[i] + 8];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 10) {
		if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
			cout << "PAPI ERROR";

		auto begin = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < imax; i += num_sum_per_oper) {
			sum += A[tmpi[i]] + A[tmpi[i] + 1] + A[tmpi[i] + 2] + A[tmpi[i] + 3] + A[tmpi[i] + 4] + A[tmpi[i] + 5] + A[tmpi[i] + 6] + A[tmpi[i] + 7] + A[tmpi[i] + 8] + A[tmpi[i] + 9];
		}
		auto end = std::chrono::high_resolution_clock::now();

		if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";

		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		ofstream REZ("Result.csv", ios::in|ios::app);
		REZ << d_type << "; "  << sum << ';' << datasize << ';' << num_sum << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << values[0] << ';' << values[1] << ';' <<  values[2] << ';' << (long double)values[2]/(time/1000) << endl ;
		REZ.close();
	}

	delete[] tmpi;
}

template<typename data_type>
void trash_to_cache(data_type* &trash) {
	int Events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM};
	long long values_trash[NUM_EVENTS];

	mt19937 gen(time(0)); 
    uniform_int_distribution<int> f(fiMin, fiMax);

	clock_t start,stop;

	data_type sum = 0;

	if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
		cout << "PAPI ERROR";

	start = clock();
	for(int i = 0; i < datasize; i++) {
		sum += trash[f(gen)] + trash[f(gen)] + trash[f(gen)] + trash[f(gen)];
	}
	stop = clock();

	if (PAPI_stop_counters(values_trash, NUM_EVENTS) != PAPI_OK)
			cout << "PAPI ERROR";
	long double time = (long double)(stop-start)/CLOCKS_PER_SEC;

	ofstream REZ("Trash.csv", ios::in|ios::app);
	REZ << sum << ";" << datasize << ';' << datasize << ';' << time << ';'  << values_trash[0] << ';' << values_trash[1] << ';' <<  values_trash[2] << ';' << values_trash[2]/time << endl ;
	REZ.close();
}


int main (int argc, char** argv) { // <Lx x - Cache level> <d/i/c - double/int/char> <num_sum - number of + operations> <num_iter>

	if(argv[2][0] == 'd') {
		mt19937 gen(time(0)); 
		uniform_real_distribution<> f(fMin_double, fMax_double); 

		mt19937 geni(time(0)); 
		uniform_int_distribution<int> fi(fiMin, fiMax);

		double *data = new double[datasize];
		double *trash;

		for(int i = 0; i < datasize; i++) {
			data[i] = f(gen);
		}

		for(double k = 1.0/32; k <= 8; k *= 2) { 
			for (int i = 0; i < atoi(argv[4]); ++i) {
				trash = new double[datasize];
				for(int k = 0; k < datasize; k++) {
					trash[k] = f(gen);
				}
				cout << data[fi(gen)] + trash[fi(gen)];

				cache_test<double>(data, argv[2][0], atoi(argv[3]), (L3size/sizeof(double) + 1) * k);

				for(int j = 0; j < number_trash_to_cache_calls; j++) {
					trash_to_cache<double>(trash);
				}

				cout << data[fi(gen)] + trash[fi(gen)] << endl;
				delete[] trash;
			}
		}
		delete[] data;
	} else if(argv[2][0] == 'i') {
		mt19937 gen(time(0)); 
	    uniform_int_distribution<int> f(fMin_int, fMax_int); 

	    mt19937 geni(time(0)); 
	    uniform_int_distribution<int> fi(fiMin, fiMax);

		int *data = new int[datasize];
		int *trash;


		for(int i = 0; i < datasize; i++) {
			data[i] = f(gen);
		}

		for(double k = 1.0/32; k <= 8; k *= 2)
			for (int i = 0; i < atoi(argv[4]); ++i) {
				trash = new int[datasize];
				for(int k = 0; k < datasize; k++) {
					trash[k] = f(gen);
				}
				cout << data[fi(gen)] + trash[fi(gen)];

				cache_test<int>(data, argv[2][0], atoi(argv[3]), (L3size/sizeof(int) + 1) * k);

				for(int j = 0; j < number_trash_to_cache_calls; j++) {
					trash_to_cache<int>(trash);
				}


				cout << data[fi(gen)] + trash[fi(gen)] << endl;
				delete[] trash;
			}

		delete[] data;
	}

	return 0;
}
