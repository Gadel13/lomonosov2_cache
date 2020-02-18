
all: main


main: main.cpp 
	g++ -O3 -o main main.cpp -L papiso/ -lpapi -lm -std=c++11
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
		main_int_1sum main_int_2sum main_int_3sum main_int_4sum
	rm -f *.o commands
	rm -f *.txt
	rm -f *.png


