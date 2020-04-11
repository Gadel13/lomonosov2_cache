#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <papi.h>
#include <random>
#include <chrono>
#include <limits>
#include <omp.h>

#include <vector>

#include <sys/time.h>


using namespace std;


#define thread_num 4

#define ch_num 2

#define number_trash_to_cache_calls 3
 

//all sizes in Bytes
long long L1size = 32*1024;
long long L2size = 256*1024;
long long L3size = 6*1024*1024;

long long RAM_size = (long long)8*1024*1024*1024;
long long datasize;
long long datasize_per_thread;

#define NUM_EVENTS 3

#define fiMin 0
#define fiMax datasize-1

#define fMin_double std::numeric_limits<double>::min()
#define fMax_double std::numeric_limits<double>::max()

#define fMin_int std::numeric_limits<int>::min()
#define fMax_int std::numeric_limits<int>::max()

mt19937 gen_cache_test_ind(time(0));

template<typename data_type>
void cache_test(vector<vector<data_type>>& A, char d_type, int num_sum_per_oper, long long num_sum, long long segment) {
	omp_set_num_threads(thread_num);

	PAPI_library_init(PAPI_VER_CURRENT);

	if (PAPI_thread_init(reinterpret_cast<long unsigned int (*)()>(omp_get_thread_num) ) != PAPI_OK)
		cout << "thread PAPI ERROR" << endl;


	uniform_int_distribution<long long> f_offset(0, datasize_per_thread/segment - 1);
	long long offset = f_offset(gen_cache_test_ind);

	uniform_int_distribution<long long> f_cache_test_ind(segment*offset, segment*(offset+1));


	vector<data_type> sum(thread_num);
	long long imax = num_sum_per_oper*(num_sum/thread_num);

	vector<vector<long long>> tmpi(thread_num, vector<long long>(imax));
	for(int th = 0; th < thread_num; th++) {
		sum[th] = 0;
		for (long long i = 0; i < imax; i++){
			tmpi[th][i] = f_cache_test_ind(gen_cache_test_ind);;
		}
	}

	long double begin[thread_num], end[thread_num];
	long long tot_L1 = 0;
	long long tot_L2 = 0;
	long long tot_L3 = 0;

	// auto begin = std::chrono::high_resolution_clock::now();
	if(num_sum_per_oper == 1) {
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}


			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}
		}


	} else if(num_sum_per_oper == 2) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 3) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 4) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 5) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 6) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 7) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 8) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 9) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	} else if(num_sum_per_oper == 10) {
		
		#pragma omp parallel
		{
			int cur_thr;
			int papi_status;
			int EventSet = PAPI_NULL;
			long long values[NUM_EVENTS];

			#pragma omp critical
			{
				if(PAPI_create_eventset(&EventSet) != PAPI_OK)
					cout << "EventSet PAPI ERROR" << endl;

				if (PAPI_add_event(EventSet, PAPI_L1_LDM) != PAPI_OK)
			    	cout << "1add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L2_LDM) != PAPI_OK)
			    	cout << "2add event PAPI ERROR" << endl;

			    if (PAPI_add_event(EventSet, PAPI_L3_LDM) != PAPI_OK)
			    	cout << "3add event PAPI ERROR" << endl;
				cur_thr = omp_get_thread_num();
				papi_status = PAPI_start(EventSet);

			}

			if (papi_status != PAPI_OK)
				cout << "start PAPI ERROR" << endl;

			#pragma omp barrier

			begin[cur_thr] = omp_get_wtime() * 1000000;
			for(long long i = 0; i < imax; i += num_sum_per_oper) {
				sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]] + A[cur_thr][tmpi[cur_thr][i+9]];
			}
			end[cur_thr] = omp_get_wtime() * 1000000;


			#pragma omp critical
			{
				papi_status = PAPI_stop(EventSet, values);

				if (papi_status != PAPI_OK)
					cout << "stop PAPI ERROR" << endl;
				tot_L1 += values[0];
				tot_L2 += values[1];
				tot_L3 += values[2];
			}

		}

	}
	// auto end = std::chrono::high_resolution_clock::now();


	// auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

	long double time = 0;

	data_type total_sum = 0;
	for(int th = 0; th < thread_num; th++){

		// long double tmp_time = (end[th].tv_sec  - begin[th].tv_sec) * 1000000 + end[th].tv_usec - begin[th].tv_usec;
		long double tmp_time = end[th] - begin[th]; 
		if(time < tmp_time){
			time = tmp_time;
		}
		total_sum += sum[th]/100000;
	}

	ofstream REZ("ResultOMP_ver2_4thr_gettimeofday.csv", ios::in|ios::app);
	REZ << d_type << "; "  << total_sum << ';' << datasize << ';' << num_sum*num_sum_per_oper << ";" << thread_num << ";" << num_sum << ";" << num_sum_per_oper << ";" << segment << ';' << time << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << "L1, L2, L3, Cache miss per second is;" << (long double)tot_L1/(time) << ';' << (long double)tot_L2/(time) << ';' << (long double)tot_L3/(time) << endl;
	REZ.close();

	return;

}


int main (int argc, char** argv) { // <Lx x - Cache level> <d/i/c - double/int/char> <num_sum - number of + operations> <num_iter>

	if(argv[2][0] == 'd') {
    	datasize = RAM_size/(sizeof(double));
    	datasize_per_thread = datasize/thread_num;

		mt19937 gen(time(0)); 
		uniform_real_distribution<> f(fMin_double, fMax_double); 

		mt19937 geni(time(0)); 
		uniform_int_distribution<int> fi(fiMin, fiMax);


		vector<vector<double>> data(thread_num, vector<double>(datasize_per_thread));
		for(int th = 0; th < thread_num; th++){
			for (long long i = 0; i < datasize_per_thread; i++){
				data[th][i] = f(gen);
			}
		}
		// double *trash;

		mt19937 tmp_geni(time(0));
		uniform_int_distribution<long long> tmp_fi(0, datasize_per_thread-1);
		for (int i = 0; i < atoi(argv[4]); ++i) {
			for(double k = 16; k >= 1.0/256; k /= 2) { 

				double tmp1 = 0;
				for (int th = 0; th < thread_num; th++) {
					tmp1 += data[th][tmp_fi(tmp_geni)]/1000000;
				}
				cout << tmp1 << "  " << endl;// + trash[fi(geni)];

				cache_test<double>(data, argv[2][0], atoi(argv[3]), k * L1size/sizeof(double), L2size/(64*thread_num));
				cache_test<double>(data, argv[2][0], atoi(argv[3]), k * L2size/sizeof(double), L3size/(64*thread_num));
				cache_test<double>(data, argv[2][0], atoi(argv[3]), k * L3size/sizeof(double), datasize_per_thread);


				double tmp2 = 0;
				for (int th = 0; th < thread_num; th++) {
					tmp2 += data[th][tmp_fi(tmp_geni)]/1000000;
				}
				cout << tmp2 << endl;// + trash[fi(gen)] << endl;
			}
		}

	} else if(argv[2][0] == 'i') {
	    datasize = RAM_size/(sizeof(int));
    	datasize_per_thread = datasize/thread_num;



		mt19937 gen(time(0)); 
	    uniform_int_distribution<int> f(fMin_int, fMax_int); 

	    mt19937 geni(time(0)); 
	    uniform_int_distribution<int> fi(fiMin, fiMax);


		vector<vector<int>> data(thread_num, vector<int>(datasize_per_thread));

		for(int th = 0; th < thread_num; th++){
			for (long long i = 0; i < datasize_per_thread; i++){
				data[th][i] = f(gen);
			}
		}
		// int *trash;


		mt19937 tmp_geni(time(0));
		uniform_int_distribution<int> tmp_fi(0, datasize_per_thread - 1);
		for (int i = 0; i < atoi(argv[4]); ++i) {
			for(double k = 16; k >= 1.0/256; k /= 2) {
				int tmp1 = 0;
				for (int th = 0; th < thread_num; th++) {
					tmp1 += data[th][tmp_fi(tmp_geni)]/1000000;
				}
				cout << tmp1 << "  " << endl;// + trash[fi(geni)];

				cache_test<int>(data, argv[2][0], atoi(argv[3]), k * L1size/sizeof(int), L2size/(64*thread_num));
				cache_test<int>(data, argv[2][0], atoi(argv[3]), k * L2size/sizeof(int), L3size/(64*thread_num));
				cache_test<int>(data, argv[2][0], atoi(argv[3]), k * L3size/sizeof(int), datasize_per_thread);

				int tmp2 = 0;
				for (int th = 0; th < thread_num; th++) {
					tmp2 += data[th][tmp_fi(tmp_geni)]/1000000;
				}
				cout << tmp2 << endl;// + trash[fi(gen)] << endl;
			}
		}
	}

	return 0;
}
