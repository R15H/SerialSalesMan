#make file
 
CC=gcc  #compiler
TARGET=tsp #target file name
 
all:
	$(CC) -O2 -lm -fopenmp -g tsp.c  -o $(TARGET) -lm

clean:
	rm $(TARGET)
