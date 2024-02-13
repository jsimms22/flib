CC = g++
OPT = -O1
CFLAGS = -mavx2

TARGETS = test

all:  $(TARGETS)

test: test.o
	$(CC) -o $@ $(LIBS) test.o

test.o: test.cpp
	$(CC) -c $(CFLAGS) test.cpp

clean:
	rm -f *.o $(TARGETS)

