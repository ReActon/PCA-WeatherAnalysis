CC=g++
CCFLAGS=-std=c++11
EFLAG=-I eigen

build:
	$(CC) $(EFLAG) Calcs.cpp -c $(CCFLAGS) 
	$(CC) Calcs.o -o pca $(CCFLAGS)
	
clean:
	rm -f pca
	rm -f *.o
	rm -f answers.txt
