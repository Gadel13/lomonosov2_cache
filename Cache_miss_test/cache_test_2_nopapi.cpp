#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <limits>
#include <omp.h>
#include <string>
#include <pthread.h>


#include <vector>

#include <sys/time.h>


using namespace std;


#define thread_num omp_get_max_threads()

#define ch_num 2

#define number_trash_to_cache_calls 3
 

//all sizes in Bytes
long long L1size = 32*1024;
long long L2size = 256*1024;
long long L3size = 35*1024*1024;

long long RAM_size = (long long)thread_num*1024*1024*1024;
long long datasize;
long long datasize_per_thread;

int cache_lv; //1 2 3
char op; //s - store l - load

#define NUM_EVENTS 3

#define fiMin 0
#define fiMax datasize-1

#define fMin_double std::numeric_limits<double>::min()
#define fMax_double std::numeric_limits<double>::max()

#define fMin_int std::numeric_limits<int>::min()
#define fMax_int std::numeric_limits<int>::max()

mt19937 gen_cache_test_ind(time(0));

template<typename data_type>
void cache_test(data_type** A, char d_type, int num_sum_per_oper, long long num_sum, long long segment) {
	// omp_set_num_threads(thread_num);

	uniform_int_distribution<long long> f_offset(0, datasize_per_thread/segment - 1);
	long long offset = f_offset(gen_cache_test_ind);

	uniform_int_distribution<long long> f_cache_test_ind(segment*offset, segment*(offset+1));


	vector<data_type> sum(thread_num);
	long long imax = num_sum_per_oper*(num_sum/thread_num);
	if (op == 't') {
		imax += (num_sum/thread_num);
	}


	if (op == 's') {
		data_type tmp_sum = 0;
		#pragma omp parallel
		{
			data_type th_tmp_sum = 0;
			int th = omp_get_thread_num();
			for(int i = segment*offset; i <= segment*(offset+1); i++){
				th_tmp_sum += A[th][i];
			}

			#pragma omp critical
			{
				tmp_sum += th_tmp_sum;
			}
		}
		cout << "tmp sum for store " << tmp_sum << endl;
	}

	long long ** tmpi = new long long* [thread_num];
	#pragma omp parallel
	{
		int th = omp_get_thread_num();
		tmpi[th] = new long long [imax];
		sum[th] = 0;
		for (long long i = 0; i < imax; i++){
			tmpi[th][i] = f_cache_test_ind(gen_cache_test_ind);;
		}
	}

	long double begin[thread_num], end[thread_num];
	long long tot_L1 = 0;
	long long tot_L2 = 0;
	long long tot_L3 = 0;
	int cur_thr;

	// auto begin = std::chrono::high_resolution_clock::now();
	cout << "start counting" << endl;
	#pragma omp parallel private (cur_thr)
	{
		cur_thr = omp_get_thread_num();
		long long values[NUM_EVENTS] = {0, 0, 0};

		#pragma omp barrier

		for(int rep = 0; rep < 5000000; rep++) {
			if(num_sum_per_oper == 1) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 2) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 3) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l')  {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 4) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l')  {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 5) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
						A[cur_thr][tmpi[cur_thr][i+4]] = i+4;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l')  {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
						A[cur_thr][tmpi[cur_thr][i+12]] = A[cur_thr][tmpi[cur_thr][i+13]] + A[cur_thr][tmpi[cur_thr][i+14]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 6) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
						A[cur_thr][tmpi[cur_thr][i+4]] = i+4;
						A[cur_thr][tmpi[cur_thr][i+5]] = i+5;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
						A[cur_thr][tmpi[cur_thr][i+12]] = A[cur_thr][tmpi[cur_thr][i+13]] + A[cur_thr][tmpi[cur_thr][i+14]];
						A[cur_thr][tmpi[cur_thr][i+15]] = A[cur_thr][tmpi[cur_thr][i+16]] + A[cur_thr][tmpi[cur_thr][i+17]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 7) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
						A[cur_thr][tmpi[cur_thr][i+4]] = i+4;
						A[cur_thr][tmpi[cur_thr][i+5]] = i+5;
						A[cur_thr][tmpi[cur_thr][i+6]] = i+6;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
						A[cur_thr][tmpi[cur_thr][i+12]] = A[cur_thr][tmpi[cur_thr][i+13]] + A[cur_thr][tmpi[cur_thr][i+14]];
						A[cur_thr][tmpi[cur_thr][i+15]] = A[cur_thr][tmpi[cur_thr][i+16]] + A[cur_thr][tmpi[cur_thr][i+17]];
						A[cur_thr][tmpi[cur_thr][i+18]] = A[cur_thr][tmpi[cur_thr][i+19]] + A[cur_thr][tmpi[cur_thr][i+20]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 8) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
						A[cur_thr][tmpi[cur_thr][i+4]] = i+4;
						A[cur_thr][tmpi[cur_thr][i+5]] = i+5;
						A[cur_thr][tmpi[cur_thr][i+6]] = i+6;
						A[cur_thr][tmpi[cur_thr][i+7]] = i+7;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
						A[cur_thr][tmpi[cur_thr][i+12]] = A[cur_thr][tmpi[cur_thr][i+13]] + A[cur_thr][tmpi[cur_thr][i+14]];
						A[cur_thr][tmpi[cur_thr][i+15]] = A[cur_thr][tmpi[cur_thr][i+16]] + A[cur_thr][tmpi[cur_thr][i+17]];
						A[cur_thr][tmpi[cur_thr][i+18]] = A[cur_thr][tmpi[cur_thr][i+19]] + A[cur_thr][tmpi[cur_thr][i+20]];
						A[cur_thr][tmpi[cur_thr][i+21]] = A[cur_thr][tmpi[cur_thr][i+22]] + A[cur_thr][tmpi[cur_thr][i+23]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 9) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
						A[cur_thr][tmpi[cur_thr][i+4]] = i+4;
						A[cur_thr][tmpi[cur_thr][i+5]] = i+5;
						A[cur_thr][tmpi[cur_thr][i+6]] = i+6;
						A[cur_thr][tmpi[cur_thr][i+7]] = i+7;
						A[cur_thr][tmpi[cur_thr][i+8]] = i+8;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
						A[cur_thr][tmpi[cur_thr][i+12]] = A[cur_thr][tmpi[cur_thr][i+13]] + A[cur_thr][tmpi[cur_thr][i+14]];
						A[cur_thr][tmpi[cur_thr][i+15]] = A[cur_thr][tmpi[cur_thr][i+16]] + A[cur_thr][tmpi[cur_thr][i+17]];
						A[cur_thr][tmpi[cur_thr][i+18]] = A[cur_thr][tmpi[cur_thr][i+19]] + A[cur_thr][tmpi[cur_thr][i+20]];
						A[cur_thr][tmpi[cur_thr][i+21]] = A[cur_thr][tmpi[cur_thr][i+22]] + A[cur_thr][tmpi[cur_thr][i+23]];
						A[cur_thr][tmpi[cur_thr][i+24]] = A[cur_thr][tmpi[cur_thr][i+25]] + A[cur_thr][tmpi[cur_thr][i+26]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}

			} else if(num_sum_per_oper == 10) {

				if(op == 's') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						A[cur_thr][tmpi[cur_thr][i]] = i;
						A[cur_thr][tmpi[cur_thr][i+1]] = i+1;
						A[cur_thr][tmpi[cur_thr][i+2]] = i+2;
						A[cur_thr][tmpi[cur_thr][i+3]] = i+3;
						A[cur_thr][tmpi[cur_thr][i+4]] = i+4;
						A[cur_thr][tmpi[cur_thr][i+5]] = i+5;
						A[cur_thr][tmpi[cur_thr][i+6]] = i+6;
						A[cur_thr][tmpi[cur_thr][i+7]] = i+7;
						A[cur_thr][tmpi[cur_thr][i+8]] = i+8;
						A[cur_thr][tmpi[cur_thr][i+9]] = i+9;
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} else if(op == 'l') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper) {
						sum[cur_thr] += A[cur_thr][tmpi[cur_thr][i]] + A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]] + A[cur_thr][tmpi[cur_thr][i+3]] + A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]] + A[cur_thr][tmpi[cur_thr][i+6]] + A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]] + A[cur_thr][tmpi[cur_thr][i+9]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				} if(op == 't') {
					begin[cur_thr] = omp_get_wtime() * 1000000;
					for(long long i = 0; i < imax; i += num_sum_per_oper*3) {
						A[cur_thr][tmpi[cur_thr][i]] = A[cur_thr][tmpi[cur_thr][i+1]] + A[cur_thr][tmpi[cur_thr][i+2]];
						A[cur_thr][tmpi[cur_thr][i+3]] = A[cur_thr][tmpi[cur_thr][i+4]] + A[cur_thr][tmpi[cur_thr][i+5]];
						A[cur_thr][tmpi[cur_thr][i+6]] = A[cur_thr][tmpi[cur_thr][i+7]] + A[cur_thr][tmpi[cur_thr][i+8]];
						A[cur_thr][tmpi[cur_thr][i+9]] = A[cur_thr][tmpi[cur_thr][i+10]] + A[cur_thr][tmpi[cur_thr][i+11]];
						A[cur_thr][tmpi[cur_thr][i+12]] = A[cur_thr][tmpi[cur_thr][i+13]] + A[cur_thr][tmpi[cur_thr][i+14]];
						A[cur_thr][tmpi[cur_thr][i+15]] = A[cur_thr][tmpi[cur_thr][i+16]] + A[cur_thr][tmpi[cur_thr][i+17]];
						A[cur_thr][tmpi[cur_thr][i+18]] = A[cur_thr][tmpi[cur_thr][i+19]] + A[cur_thr][tmpi[cur_thr][i+20]];
						A[cur_thr][tmpi[cur_thr][i+21]] = A[cur_thr][tmpi[cur_thr][i+22]] + A[cur_thr][tmpi[cur_thr][i+23]];
						A[cur_thr][tmpi[cur_thr][i+24]] = A[cur_thr][tmpi[cur_thr][i+25]] + A[cur_thr][tmpi[cur_thr][i+26]];
						A[cur_thr][tmpi[cur_thr][i+27]] = A[cur_thr][tmpi[cur_thr][i+28]] + A[cur_thr][tmpi[cur_thr][i+29]];
					}
					end[cur_thr] = omp_get_wtime() * 1000000;
				}
			}
		}

		#pragma omp critical
		{
			tot_L1 += values[0];
			tot_L2 += values[1];
			tot_L3 += values[2];
		}
	}
	cout <<"stop counting" << endl;
	// auto end = std::chrono::high_resolution_clock::now();


	// auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

	long double time = 0.0;

	data_type total_sum = 0;
	for(int th = 0; th < thread_num; th++){

		// long double tmp_time = (end[th].tv_sec  - begin[th].tv_sec) * 1000000 + end[th].tv_usec - begin[th].tv_usec;
		long double tmp_time = end[th] - begin[th]; 
		if(time < tmp_time){
			time = tmp_time;
		}
		total_sum += sum[th];
	}

	string filename;
	if (op == 's') {
		filename = "nopapi_store_lom2_ResultOMP_ver2_";
	} else if (op == 'l') {
		filename = "nopapi_load_lom2_ResultOMP_ver2_";
	} if (op == 't') {
		filename = "nopapi_total_lom2_ResultOMP_ver2_";
	}
	filename += to_string(cache_lv);
	filename += "lvl_";
	filename += to_string(thread_num);
	filename += "thr.csv";

	ofstream REZ(filename, ios::in|ios::app);
	REZ << cache_lv << ";" << d_type << ";"  << total_sum << ";" << datasize << ";" << num_sum*num_sum_per_oper << ";" << thread_num << ";" << num_sum << ";" << num_sum_per_oper << ";" << segment << ';' << time << ';'  << tot_L1 << ';' << tot_L2 << ';' <<  tot_L3 << ';' << (long double)tot_L1/(time) << ';' << (long double)tot_L2/(time) << ';' << (long double)tot_L3/(time) << endl;
	REZ.close();

	for(int th = 0; th < thread_num; th++) {
		delete[] tmpi[th];
	}
	delete[] tmpi;

	return;

}


int main (int argc, char** argv) { // <Lx x - Cache level> <s/l - store load test> <d/i/c - double/int/char> <num_sum - number of + operations> <num_iter>

	cache_lv = atoi(argv[1]);
	op = argv[2][0];

	if(argv[3][0] == 'd') {
    	datasize = RAM_size/(sizeof(double));
    	datasize_per_thread = datasize/thread_num;

		mt19937 gen(time(0)); 
		uniform_real_distribution<> f(fMin_double, fMax_double); 

		mt19937 geni(time(0)); 
		uniform_int_distribution<int> fi(fiMin, fiMax);


		double ** data = new double* [thread_num];
		for(int th = 0; th < thread_num; th++){
			data[th] = new double [datasize_per_thread];
			for (long long i = 0; i < datasize_per_thread; i++){
				data[th][i] = f(gen);
			}
		}

		mt19937 tmp_geni(time(0));
		uniform_int_distribution<long long> tmp_fi(0, datasize_per_thread-1);

		double tmp1 = 0;
		for (int th = 0; th < thread_num; th++) {
			tmp1 += data[th][tmp_fi(tmp_geni)]/1000000;
		}
		cout << tmp1 << "  " << endl;

		if(cache_lv == 1){
			cache_test<double>(data, argv[3][0], atoi(argv[4]), L2size*thread_num/(sizeof(double)), L2size/sizeof(double));
		} else if(cache_lv == 2){
			cache_test<double>(data, argv[3][0], atoi(argv[4]), L3size/(sizeof(double)), L3size/sizeof(double));
		} else {
			cache_test<double>(data, argv[3][0], atoi(argv[4]), L3size*thread_num/(sizeof(double)), datasize_per_thread);
		}

		double tmp2 = 0;
		for (int th = 0; th < thread_num; th++) {
			tmp2 += data[th][tmp_fi(tmp_geni)]/1000000;
		}
		cout << tmp2 << endl;

		for(int th = 0; th < thread_num; th++){
			delete[] data[th];
		}
		delete[] data;

	} else if(argv[3][0] == 'i') {
	    datasize = RAM_size/(sizeof(int));
    	datasize_per_thread = datasize/thread_num;



		mt19937 gen(time(0)); 
	    uniform_int_distribution<int> f(fMin_int, fMax_int); 

	    mt19937 geni(time(0)); 
	    uniform_int_distribution<int> fi(fiMin, fiMax);


		int** data = new int* [thread_num];

		for(int th = 0; th < thread_num; th++){
			data[th] = new int [datasize_per_thread];
			for (long long i = 0; i < datasize_per_thread; i++){
				data[th][i] = f(gen);
			}
		}

		mt19937 tmp_geni(time(0));
		uniform_int_distribution<int> tmp_fi(0, datasize_per_thread - 1);

		int tmp1 = 0;
		for (int th = 0; th < thread_num; th++) {
			tmp1 += data[th][tmp_fi(tmp_geni)]/1000000;
		}
		cout << tmp1 << "  " << endl;

		if(cache_lv == 1){
			cache_test<int>(data, argv[3][0], atoi(argv[4]), L2size*thread_num/(sizeof(int)), L2size/sizeof(int));
		} else if(cache_lv == 2){
			cache_test<int>(data, argv[3][0], atoi(argv[4]), L3size/(sizeof(int)), L3size/sizeof(int));
		} else {
			cache_test<int>(data, argv[3][0], atoi(argv[4]), L3size*thread_num/(sizeof(int)), datasize_per_thread);
		}

		int tmp2 = 0;
		for (int th = 0; th < thread_num; th++) {
			tmp2 += data[th][tmp_fi(tmp_geni)]/1000000;
		}
		cout << tmp2 << endl;

		for(int th = 0; th < thread_num; th++){
			delete[] data[th];
		}
		delete[] data;
	}

	return 0;
}
