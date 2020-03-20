#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <papi.h>
#include <random>
#include <chrono>
#include <limits>
#include <omp.h>


using namespace std;


#define thread_num 2

 
//all sizes in Bytes
#define L1size 32000
#define L2size 256000
#define L3size 6000000

#define RAM_size 15000000000
#define ch_num 2

#define number_trash_to_cache_calls 3

long long datasize;

#define NUM_EVENTS 3

#define fiMin 0
#define fiMax datasize-1

#define fMin_double std::numeric_limits<double>::min()
#define fMax_double std::numeric_limits<double>::max()

#define fMin_int std::numeric_limits<int>::min()
#define fMax_int std::numeric_limits<int>::max()

mt19937 gen_cache_test_ind(time(0)); 

template<typename data_type>
void cache_test(data_type** &A, char d_type, int num_sum_per_oper, int num_sum) {
	PAPI_library_init(PAPI_VER_CURRENT);
	int Events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM};
	long long **values;
	values = new long long *[thread_num];


	omp_set_num_threads(thread_num);
	

	uniform_int_distribution<int> f_cache_test_ind(0, datasize/thread_num - 1);



	data_type* sum = new data_type[thread_num];
	int imax = (num_sum_per_oper*num_sum)/thread_num;

	int** tmpi = new int*[thread_num];
	for(int th = 0; th < thread_num; th++) {
		values[th] = new long long [NUM_EVENTS];
		sum[th] = 0;
		tmpi[th] = new int [imax];

		for (int i = 0; i < imax; i++){
			tmpi[th][i] = f_cache_test_ind(gen_cache_test_ind);
		}
	}

	if(num_sum_per_oper == 1) {

		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	}if(num_sum_per_oper == 2) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 3) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 4) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 5) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 6) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 7) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 8) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 9) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	} if(num_sum_per_oper == 10) {
		auto begin = std::chrono::high_resolution_clock::now();
		
		int cur_thr;
		#pragma omp parallel private(cur_thr)
		{
			cur_thr = omp_get_thread_num();

			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK) 
				cout << "PAPI ERROR";

			for(int i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]] + A[cur_thr][tmpi[cur_thr][i+9]];
			}

			if (PAPI_stop_counters(values[cur_thr], NUM_EVENTS) != PAPI_OK)
				cout << "PAPI ERROR";
		}

		auto end = std::chrono::high_resolution_clock::now();


		auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

		data_type total_sum = 0;
		long long tot_L1 = 0;
		long long tot_L2;
		long long tot_L3;
		for(int th = 0; th < thread_num; th++){
			total_sum += sum[th];
			tot_L1 += values[th][0];
			tot_L2 += values[th][1];
			tot_L3 += values[th][2];
		}

		ofstream REZ("ResultOMP_ver2_2thr.csv", ios::in|ios::app);
		REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum << ";" << thread_num << ";" << num_sum_per_oper << ';' << time/(1000) << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L3/(time/1000) << endl ;
		REZ.close();
	}

	for(int th = 0; th < thread_num; th++){
		delete[] tmpi[th];
		delete[] values[th];
	}

	delete[] tmpi;
	delete[] values;
	delete[] sum;
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
	#pragma omp parallel for
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
		// datasize = 2*RAM_size/(5*sizeof(double)); 
		datasize = RAM_size/(sizeof(double));

		mt19937 gen(time(0)); 
		uniform_real_distribution<> f(fMin_double, fMax_double); 

		mt19937 geni(time(0)); 
		uniform_int_distribution<int> fi(fiMin, fiMax);


		double **data = new double* [thread_num];
		for(int th = 0; th < thread_num; th++){
			data[th] = new double [datasize/thread_num];
		}
		// double *trash;

		for(int th = 0; th < thread_num; th++) {
			for (int i = 0; i < datasize/thread_num; i++){
				data[th][i] = f(gen);
			}
		}

		for(double k = 1.0/32; k <= 8; k *= 2) { 
			for (int i = 0; i < atoi(argv[4]); ++i) {
				// trash = new double[datasize];
				// for(int k = 0; k < datasize; k++) {
				// 	trash[k] = f(gen);
				// }

				double tmp = 0;
				mt19937 tmp_geni(time(0));
				uniform_int_distribution<int> tmp_fi(0, datasize/thread_num - 1);
				for (int th = 0; th < thread_num; th++) {
					tmp += data[th][tmp_fi(tmp_geni)];
				}
				cout << tmp;// + trash[fi(geni)];

				cache_test<double>(data, argv[2][0], atoi(argv[3]), (L3size/sizeof(double) + 1) * k);

				// for(int j = 0; j < number_trash_to_cache_calls; j++) {
				// 	trash_to_cache<double>(trash);
				// }

				int tmp2 = 0;
				mt19937 tmp2_geni(time(0));
				uniform_int_distribution<int> tmp2_fi(0, datasize/thread_num - 1);
				for (int th = 0; th < thread_num; th++) {
					tmp2 += data[th][tmp2_fi(tmp2_geni)];
				}
				cout << tmp2 << endl;// + trash[fi(gen)] << endl;
				// delete[] trash;
			}
		}
		for(int th = 0; th < thread_num; th++){
			delete[] data[th];
		}

		delete[] data;

	} else if(argv[2][0] == 'i') {
	    // datasize = 2*RAM_size/(5*sizeof(int));

	    datasize = RAM_size/(sizeof(int));

		mt19937 gen(time(0)); 
	    uniform_int_distribution<int> f(fMin_int, fMax_int); 

	    mt19937 geni(time(0)); 
	    uniform_int_distribution<int> fi(fiMin, fiMax);


		int **data = new int* [thread_num];
		for(int th = 0; th < thread_num; th++){
			data[th] = new int [datasize/thread_num];
		}
		// int *trash;


		for(int th = 0; th < thread_num; th++) {
			for (int i = 0; i < datasize/thread_num; i++){
				data[th][i] = f(gen);
			}
		}

		mt19937 tmp_geni(time(0));
		uniform_int_distribution<int> tmp_fi(0, datasize/thread_num - 1);
		for(double k = 1.0/32; k <= 8; k *= 2)
			for (int i = 0; i < atoi(argv[4]); ++i) {
				// trash = new int[datasize];
				// for(int k = 0; k < datasize; k++) {
					// trash[k] = f(gen);
				// }

				int tmp1 = 0;
				for (int th = 0; th < thread_num; th++) {
					tmp1 += data[th][tmp_fi(tmp_geni)];
				}
				cout << tmp1;// + trash[fi(geni)];

				cache_test<int>(data, argv[2][0], atoi(argv[3]), (L3size/sizeof(int) + 1) * k);

				// for(int j = 0; j < number_trash_to_cache_calls; j++) {
				// 	trash_to_cache<int>(trash);
				// }


				int tmp2 = 0;
				for (int th = 0; th < thread_num; th++) {
					tmp2 += data[th][tmp_fi(tmp_geni)];
				}
				cout << tmp2 << endl;// + trash[fi(gen)] << endl;
				// delete[] trash;
			}

		for(int th = 0; th < thread_num; th++){
			delete[] data[th];
		}

		delete[] data;
	}

	return 0;
}
