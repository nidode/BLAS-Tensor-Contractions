include $(INCLUDE_FILE)

# General flags for compilation:
CFLAGS    += -fopenmp
LDFLAGS   += -fopenmp

# Dependencies:
TEST_CUBMAT = test_cubmat.o time_cubmat.o \
              case1.o case2.o case3.o case4.o case5.o \
              extract_flat_matrices.o \
              MyFLA_Utils.o

all : test_cubmat.x

test_cubmat.x : $(TEST_CUBMAT) Makefile
	$(LINKER) $(LDFLAGS) -o $@ $(TEST_CUBMAT) $(LDLIBS)

clean : 
	rm -f *.o *~ core *.x

