
all: main


main: main.cpp 
	g++ -O3 -o main main.cpp -L"./papiso" -lpapi -lm -std=c++11



test: 
	./main 5
	

clean:
	rm -f *.o main gen print eql \
		A.dat B.dat C.dat testC.dat
	rm -f *.o commands
	rm -f *.txt
	rm -f *.png


