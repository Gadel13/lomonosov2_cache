
all: main


main: main.cpp 
	g++ -O3 -o cache_test cache_test.cpp -L /usr/local/lib/ -lpapi -lm -std=c++11 -fopenmp -lpthread
	g++ -O3 -o main main.cpp -L papiso/ -lpapi -lm -lpapi -std=c++11
	g++ -O3 -o main_double_2sum main_double_2sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_double_3sum main_double_3sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_double_4sum main_double_4sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_double_5sum main_double_5sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_double_6sum main_double_6sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_double_7sum main_double_7sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_1sum main_int_1sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_2sum main_int_2sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_3sum main_int_3sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_4sum main_int_4sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_5sum main_int_5sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_6sum main_int_6sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_int_7sum main_int_7sum.cpp -L papiso/ -lpapi -lm -std=c++11
	g++ -O3 -o main_without_papi main_without_papi.cpp -L papiso/ -lpapi -lm -std=c++11

multi_test:
	./cache_test 3 i 1 3
	./cache_test 3 i 2 3
	./cache_test 3 i 3 3
	./cache_test 3 i 4 3
	./cache_test 3 i 5 3
	./cache_test 3 i 6 3
	./cache_test 3 i 7 3
	./cache_test 3 i 8 3
	./cache_test 3 i 9 3
	./cache_test 3 i 10 3

	./cache_test 3 d 1 3
	./cache_test 3 d 2 3
	./cache_test 3 d 3 3
	./cache_test 3 d 4 3
	./cache_test 3 d 5 3
	./cache_test 3 d 6 3
	./cache_test 3 d 7 3
	./cache_test 3 d 8 3
	./cache_test 3 d 9 3
	./cache_test 3 d 10 3

test_int: 
	./main_int_1sum 5;
	./main_int_2sum 5;
	./main_int_3sum 5;
	./main_int_4sum 5;
	./main_int_5sum 5;
	./main_int_6sum 5;

test_double:
	./main 5;
	./main_double_2sum 5;
	./main_double_3sum 5;
	./main_double_4sum 5;
	./main_double_5sum 5;
	./main_double_6sum 5;

test:
	./main_int_5sum 5;
	./main_double_5sum 5;
	./main_int_6sum 5;
	./main_double_6sum 5;


	./main_int_1sum 5;
	./main_int_2sum 5;
	./main_int_3sum 5;
	./main_int_4sum 5;
	
	./main 5;
	./main_double_2sum 5;
	./main_double_3sum 5;
	./main_double_4sum 5;

test1: 
	./main_int_1sum 5;
	./main_int_2sum 5;
	./main_int_3sum 5;
	./main_int_4sum 5;
	./main_int_5sum 5;
	./main_int_6sum 5;

	./main 5;
	./main_double_2sum 5;
	./main_double_3sum 5;
	./main_double_4sum 5;
	./main_double_5sum 5;
	./main_double_6sum 5;
	

clean:
	rm -f *.o main main_double_2sum main_double_3sum main_double_4sum \
		main_int_1sum main_int_2sum main_int_3sum main_int_4sum \
		main_int_5sum main_int_6sum main_int_7sum main_double_5sum \
		main_double_6sum main_double_7sum main_without_papi cache_test
	rm -f *.o commands
	rm -f *.txt
	rm -f *.png


