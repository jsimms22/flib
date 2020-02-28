CC = g++
OPT = -O1
CFLAGS = -mavx2  -DMKL_ILP64 -m64 -I${MKLROOT}/include $(OPT)

LDLIBS = -Wl,--start-group \
${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a \
${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group \
-lpthread -lm -ldl

#LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
#LDLIBS = -lrt  -I$(MKLROOT)/include -Wl,-L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

TARGETS = test

all:  $(TARGETS)

test: test.o
	$(CC) -o $@ $(LIBS) test.o

test.o: test.cpp
	$(CC) -c $(CFLAGS) test.cpp

clean:
	rm -f *.o $(TARGETS) *.stdouts *.txt

