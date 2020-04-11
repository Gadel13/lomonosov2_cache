
all: main


main: main.cpp 
	g++ -O3 -o cache_test_2 cache_test_2.cpp -L papiso/ -lpapi -lm -std=c++11 -fopenmp -lpthread
	
qqq:
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
	./cache_test_2 3 i 1 3
	./cache_test_2 3 i 2 3
	./cache_test_2 3 i 3 3
	./cache_test_2 3 i 4 3
	./cache_test_2 3 i 5 3
	./cache_test_2 3 i 6 3
	./cache_test_2 3 i 7 3
	./cache_test_2 3 i 8 3
	./cache_test_2 3 i 9 3
	./cache_test_2 3 i 10 3

	./cache_test_2 3 d 1 3
	./cache_test_2 3 d 2 3
	./cache_test_2 3 d 3 3
	./cache_test_2 3 d 4 3
	./cache_test_2 3 d 5 3
	./cache_test_2 3 d 6 3
	./cache_test_2 3 d 7 3
	./cache_test_2 3 d 8 3
	./cache_test_2 3 d 9 3
	./cache_test_2 3 d 10 3


t1:
	./cache_test_2 3 i 1 1
	./cache_test_2 3 d 1 1

t2:
	./cache_test_2 3 i 2 1
	./cache_test_2 3 d 2 1

t3:
	./cache_test_2 3 i 3 1
	./cache_test_2 3 d 3 1

t4:
	./cache_test_2 3 i 4 1
	./cache_test_2 3 d 4 1

t5:
	./cache_test_2 3 i 5 1
	./cache_test_2 3 d 5 1

t6:
	./cache_test_2 3 i 6 1
	./cache_test_2 3 d 6 1

t7:
	./cache_test_2 3 i 7 1
	./cache_test_2 3 d 7 1

t8:
	./cache_test_2 3 i 8 1
	./cache_test_2 3 d 8 1

t9:
	./cache_test_2 3 i 9 1
	./cache_test_2 3 d 9 1

t10:
	./cache_test_2 3 i 10 1
	./cache_test_2 3 d 10 1

multi_test_v2:
	make t10
	make t7
	make t8
	make t9
	make t5
	make t6
	make t1
	make t2
	make t3
	make t4


	make t1
	make t2
	make t3
	make t4
	make t5
	make t6
	make t7
	make t8
	make t9
	make t10



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
		main_double_6sum main_double_7sum main_without_papi cache_test cache_test_2
	rm -f *.o commands
	rm -f *.txt
	rm -f *.png


